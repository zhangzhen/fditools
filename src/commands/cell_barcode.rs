use anyhow::Result;
use bio::io::fastq;
use bk_tree::{BKTree, Metric};
use clap::Parser;
use crossbeam_channel::{bounded, Receiver};
use fditools_lib::{
    par_chunks::{self, ChunkProcessor},
    utils,
};
use log::*;
use serde::Serialize;
use std::{
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str, thread,
};
use triple_accel::hamming;

const DISTANCE_THRESHOLD: u32 = 1;

#[derive(Serialize)]
struct ReadPair(fastq::Record, fastq::Record);

pub struct Hamming;

impl<K: AsRef<[u8]> + ?Sized> Metric<K> for Hamming {
    fn distance(&self, a: &K, b: &K) -> u32 {
        let a_bytes = a.as_ref();
        let b_bytes = b.as_ref();
        hamming(a_bytes, b_bytes)
    }

    fn threshold_distance(&self, a: &K, b: &K, threshold: u32) -> Option<u32> {
        let a_bytes = a.as_ref();
        let b_bytes = b.as_ref();
        let distance = hamming(a_bytes, b_bytes);
        if distance <= threshold {
            Some(distance)
        } else {
            None
        }
    }
}

#[derive(Parser, Debug)]
pub struct CellBarcode {
    /// Input fastq file R1
    #[arg(short = 'a', long)]
    in_reads1: PathBuf,

    /// Input fastq file R2
    #[arg(short = 'b', long)]
    in_reads2: PathBuf,

    /// Output fastq file R1
    #[arg(short = 'c', long)]
    out_reads1: PathBuf,

    /// Output fastq file R2
    #[arg(short = 'd', long)]
    out_reads2: PathBuf,

    #[arg(short = 'w', long)]
    bc_whitelist: PathBuf,

    /// Number of CPU cores to use
    #[arg(short = 't', long, default_value = utils::NUM_CPU.as_str())]
    threads: usize,

    #[arg(short = 'n', long, default_value = par_chunks::NUM_CHUNKS_STR.as_str())]
    num_chunks: usize,

    #[arg(short = 'C', long, default_value = par_chunks::CHUNKSIZE_STR.as_str())]
    chunksize: usize,

    /// barcode match length at the 5' end
    #[arg(short = 'm', long)]
    match_len: usize,

    /// skip the number of bases following the barcode at the 5' end
    #[arg(short = 's', long)]
    skip_len: usize,

    /// Optionally b-gzip the output.
    #[arg(short = 'Z', long)]
    gzip: bool,

    /// The number of threads to use for compressing output (specified by --gzip)
    #[arg(short = 'T', long, default_value = "4")]
    compression_threads: usize,

    /// The level to use for compressing output (specified by --gzip)
    #[arg(short = 'L', long, default_value = "6")]
    compression_level: u32,
}

impl CellBarcode {
    pub fn run(self) -> Result<()> {
        // Spawn a thread to read and chunk the FASTQ files
        let in_reads1 = self.in_reads1.clone();
        let in_reads2 = self.in_reads2.clone();

        let chunk_rx = load_read_pairs(
            in_reads1,
            in_reads2,
            self.chunksize,
            self.num_chunks,
            self.gzip,
        )?;

        let whitelist = load_barcodes(self.bc_whitelist)?;

        let processor = CellBarcodeProcessor::new(whitelist, self.match_len, self.skip_len);

        let par_chunks_runner = par_chunks::ParChunks::new(
            Some(std::cmp::max(self.threads.checked_sub(1).unwrap_or(0), 1)),
            Some(self.chunksize),
            chunk_rx,
            Some(self.num_chunks),
            processor,
        );

        let mut writer1 = utils::get_writer(
            &Some(self.out_reads1),
            self.gzip,
            self.compression_threads,
            self.compression_level,
        )?;

        let mut writer2 = utils::get_writer(
            &Some(self.out_reads2),
            self.gzip,
            self.compression_threads,
            self.compression_level,
        )?;

        let result_rx = par_chunks_runner.process()?;
        for ReadPair(read1, read2) in result_rx {
            writer1.write_record(&read1)?;
            writer2.write_record(&read2)?;
        }
        writer1.flush()?;
        writer2.flush()?;

        Ok(())
    }
}

fn load_read_pairs(
    reads1: PathBuf,
    reads2: PathBuf,
    chunksize: usize,
    num_chunks: usize,
    gzip: bool,
) -> Result<Receiver<Vec<ReadPair>>> {
    info!("Start loading read pairs");
    let (chunk_tx, chunk_rx) = bounded(num_chunks);

    thread::spawn(move || {
        let reader1 =
            utils::get_reader(&Some(reads1), gzip).expect("Failed to open read1 FASTQ file");
        let reader2 =
            utils::get_reader(&Some(reads2), gzip).expect("Failed to open read2 FASTQ file");

        // Zip the two readers to process paired reads
        let paired_reads = reader1.records().zip(reader2.records());

        let mut current_chunk = Vec::with_capacity(chunksize);
        for (read1_result, read2_result) in paired_reads {
            match (read1_result, read2_result) {
                (Ok(read1), Ok(read2)) => {
                    // Here you would need to convert the paired reads to R::P
                    // This is a placeholder - actual implementation depends on the R::P type
                    current_chunk.push(ReadPair(read1, read2));

                    // If the chunk is full, send it for processing
                    if current_chunk.len() >= chunksize {
                        chunk_tx
                            .send(current_chunk)
                            .expect("Sent a serializable to processor");
                        current_chunk = Vec::with_capacity(chunksize);
                    }
                }
                _ => {
                    error!("Error reading paired FASTQ records");
                    // Handle the error as appropriate for your application
                }
            }
        }

        // Send any remaining reads in the final chunk
        if !current_chunk.is_empty() {
            chunk_tx
                .send(current_chunk)
                .expect("Failed to send final chunk - channel closed");
        }

        info!("Finished loading read pairs");
    });

    Ok(chunk_rx)
}

struct CellBarcodeProcessor {
    whitelist: BKTree<Vec<u8>, Hamming>,
    match_len: usize,
    skip_len: usize,
}

impl CellBarcodeProcessor {
    pub fn new(whitelist: Vec<Vec<u8>>, match_len: usize, skip_len: usize) -> Self {
        let mut tree = BKTree::new(Hamming);
        for barcode in whitelist {
            tree.add(barcode);
        }
        Self {
            whitelist: tree,
            match_len,
            skip_len,
        }
    }
}

impl ChunkProcessor for CellBarcodeProcessor {
    type P = ReadPair;

    fn process_chunk(&self, chunk: Vec<ReadPair>) -> Vec<ReadPair> {
        trace!("Processing a chunk of size {}", chunk.len());
        let mut results = Vec::with_capacity(chunk.len());
        for ReadPair(read1, read2) in chunk {
            // Process each read pair as needed
            let found1 = self
                .whitelist
                .find(&read1.seq()[..self.match_len], DISTANCE_THRESHOLD)
                .collect::<Vec<_>>();
            let found2 = self
                .whitelist
                .find(&read2.seq()[..self.match_len], DISTANCE_THRESHOLD)
                .collect::<Vec<_>>();
            if found1.len() != 1 || found2.len() != 1 {
                continue;
            }

            let (_, bc1) = found1[0];
            let (_, bc2) = found2[0];
            let desc1 = format!(
                "{} CB:Z:{}{}",
                read1.desc().unwrap_or(""),
                str::from_utf8(bc1).unwrap(),
                str::from_utf8(bc2).unwrap()
            );
            let read1_new = fastq::Record::with_attrs(
                read1.id(),
                Some(&desc1),
                &read1.seq()[self.match_len + self.skip_len..],
                &read1.qual()[self.match_len + self.skip_len..],
            );
            let desc2 = format!(
                "{} CB:Z:{}{}",
                read2.desc().unwrap_or(""),
                str::from_utf8(bc1).unwrap(),
                str::from_utf8(bc2).unwrap()
            );
            let read2_new = fastq::Record::with_attrs(
                read2.id(),
                Some(&desc2),
                &read2.seq()[self.match_len + self.skip_len..],
                &read2.qual()[self.match_len + self.skip_len..],
            );
            results.push(ReadPair(read1_new, read2_new));
        }
        results
    }
}

fn load_barcodes(bc_file: impl AsRef<Path>) -> Result<Vec<Vec<u8>>> {
    let file = File::open(bc_file)?;
    let reader = BufReader::new(file);
    let mut barcodes = Vec::new();

    for line in reader.lines() {
        if let Ok(barcode) = line {
            let barcode = barcode.trim();
            if !barcode.is_empty() {
                barcodes.push(barcode.as_bytes().to_owned());
            }
        }
    }

    Ok(barcodes)
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use super::*;
    use fditools_lib::{
        par_chunks::{self, CHUNKSIZE_STR, NUM_CHUNKS_STR},
        utils,
    };
    use rstest::*;

    #[fixture]
    fn barcodes() -> Vec<Vec<u8>> {
        load_barcodes("test_data/barcodes.txt").unwrap()
    }

    #[fixture]
    fn read_pairs(barcodes: Vec<Vec<u8>>) -> Vec<ReadPair> {
        let chunksize = CHUNKSIZE_STR.parse().unwrap();
        let num_chunks = NUM_CHUNKS_STR.parse().unwrap();
        let chunk_rx = load_read_pairs(
            PathBuf::from("test_data/r1.fastq.gz"),
            PathBuf::from("test_data/r2.fastq.gz"),
            chunksize,
            num_chunks,
            true,
        )
        .unwrap();

        let processor = CellBarcodeProcessor::new(barcodes, 11, 19);

        let par_chunks_runner = par_chunks::ParChunks::new(
            Some(2),
            Some(chunksize),
            chunk_rx,
            Some(num_chunks),
            processor,
        );
        let mut records = vec![];
        par_chunks_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|record| {
                records.push(record);
            });
        records
    }

    #[rstest]
    fn test_read_pairs(read_pairs: Vec<ReadPair>) {
        assert_eq!(read_pairs.len(), 10);
        let rp = read_pairs.first().unwrap();

        let desc1 = rp.0.desc().unwrap();
        let desc1_parts: Vec<&str> = desc1.split(' ').collect();
        let cb_tag1 = desc1_parts[desc1_parts.len() - 1];
        assert_eq!(cb_tag1, "CB:Z:TGAAGTATTGTTATATTGTGAG");

        let desc2 = rp.1.desc().unwrap();
        let desc2_parts: Vec<&str> = desc2.split(' ').collect();
        let cb_tag2 = desc2_parts[desc2_parts.len() - 1];
        assert_eq!(cb_tag2, cb_tag1);
    }
}
