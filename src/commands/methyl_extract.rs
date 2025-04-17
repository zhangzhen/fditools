use anyhow::{bail, Context, Result};
use bio::alphabets::dna::revcomp;
use clap::Parser;
use fditools_lib::par_granges::{ParGranges, RegionProcessor};
use rayon::prelude::*;
use rust_htslib::{
    bam,
    bam::{
        record::{Aux, Cigar},
        Read,
    },
};
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct MethylExtractArgs {
    /// Output file path
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Input BAM/SAM alignment file
    #[arg(short = 'a', long)]
    alignment: PathBuf,

    /// Number of CPU cores to use
    #[arg(short = 'c', long, default_value = "1")]
    cpu: usize,

    /// Trimming for Tn5
    #[arg(short = 't', long, default_value = "15")]
    trimming_tn5: i32,

    /// Mapping quality filter threshold
    #[arg(short = 'q', long, default_value = "5")]
    map_quality_filter: u8,

    /// Whether this is ATAC-seq data
    #[arg(long, default_value = "Y")]
    atac: String,

    /// Whether to remove chimeric reads
    #[arg(long, default_value = "Y")]
    remove_chimera: String,
}

lazy_static::lazy_static! {
    static ref TAG_DICT: HashMap<&'static str, char> = {
        let mut m = HashMap::new();
        m.insert("ACA", 'A'); m.insert("ACT", 'B'); m.insert("ACC", 'C'); m.insert("ACG", 'D');
        m.insert("TCA", 'E'); m.insert("TCT", 'F'); m.insert("TCC", 'G'); m.insert("TCG", 'H');
        m.insert("CCA", 'I'); m.insert("CCT", 'J'); m.insert("CCC", 'K'); m.insert("CCG", 'L');
        m.insert("GCA", 'M'); m.insert("GCT", 'N'); m.insert("GCC", 'O'); m.insert("GCG", 'P');
        m
    };
}

#[derive(Debug, Clone)]
struct ReadPair {
    read1: bam::Record,
    read2: bam::Record,
}

impl ReadPair {
    fn new(read1: bam::Record, read2: bam::Record) -> Self {
        Self { read1, read2 }
    }
}

fn get_context(seq: &[u8], start: usize, end: usize) -> char {
    if start >= seq.len() || end >= seq.len() {
        return '.';
    }

    let context = std::str::from_utf8(&seq[start..=end])
        .unwrap_or("...")
        .to_uppercase();

    if context.contains('N') {
        return 'U';
    }

    let context = if context.as_bytes()[1] == b'G' {
        String::from_utf8(revcomp(context.as_bytes().to_vec())).unwrap_or_default()
    } else {
        context
    };

    TAG_DICT.get(context.as_str()).copied().unwrap_or('.')
}

fn call_methyl(
    read_seq: &[u8],
    ref_seq: &[u8],
    ref_pos: &[Option<i64>],
    ref_start: u32,
    xg: &str,
) -> (String, Vec<usize>, Vec<usize>) {
    let mut is_methylated = String::new();
    let mut ga_converted_ref_pos = Vec::new();
    let mut ct_converted_ref_pos = Vec::new();
    let mut j = ref_start as usize;

    for (i, pos) in ref_pos.iter().enumerate() {
        match pos {
            None => {
                if i == read_seq.len() - 1 {
                    is_methylated.push('.');
                }
                continue;
            }
            Some(pos) => {
                let pos = usize::try_from(*pos).expect("check overflow");
                is_methylated.extend(std::iter::repeat('.').take(pos - j - 1));
                j = pos;

                let ref_base = ref_seq[pos - ref_start as usize];
                let read_base = read_seq[i];

                match xg {
                    "GA" => match (ref_base, read_base) {
                        (b'g', b'A') => {
                            let context = get_context(
                                ref_seq,
                                pos - ref_start as usize - 1,
                                pos - ref_start as usize + 1,
                            );
                            is_methylated.push(context.to_ascii_lowercase());
                            ga_converted_ref_pos.push(j);
                        }
                        (b'G', b'G') => {
                            let context = get_context(
                                ref_seq,
                                pos - ref_start as usize - 1,
                                pos - ref_start as usize + 1,
                            );
                            is_methylated.push(context.to_ascii_uppercase());
                        }
                        (b'c', b'T') => {
                            is_methylated.push('.');
                            ct_converted_ref_pos.push(j);
                        }
                        _ => is_methylated.push('.'),
                    },
                    "CT" => match (ref_base, read_base) {
                        (b'c', b'T') => {
                            let context = get_context(
                                ref_seq,
                                pos - ref_start as usize - 1,
                                pos - ref_start as usize + 1,
                            );
                            is_methylated.push(context.to_ascii_lowercase());
                            ct_converted_ref_pos.push(j);
                        }
                        (b'C', b'C') => {
                            let context = get_context(
                                ref_seq,
                                pos - ref_start as usize - 1,
                                pos - ref_start as usize + 1,
                            );
                            is_methylated.push(context.to_ascii_uppercase());
                        }
                        (b'g', b'A') => {
                            is_methylated.push('.');
                            ga_converted_ref_pos.push(j);
                        }
                        _ => is_methylated.push('.'),
                    },
                    _ => is_methylated.push('.'),
                }
            }
        }
    }

    (is_methylated, ga_converted_ref_pos, ct_converted_ref_pos)
}

/// Return expanded sequence from MD tag.
///
/// The sequence includes substitutions and both insertions in the
/// reference as well as deletions to the reference sequence. Combine
/// with the cigar string to reconstitute the query or the reference
/// sequence.
///
/// Positions corresponding to `N` (skipped region from the reference)
/// or `P` (padding (silent deletion from padded reference)) in the CIGAR
/// string will not appear in the returned sequence. The MD should
/// correspondingly not contain these. Thus proper tags are:
///
///     Deletion from the reference:    cigar=5M1D5M    MD=5^C5
///     Skipped region from reference:  cigar=5M1N5M    MD=10
///     Padded region in the reference: cigar=5M1P5M    MD=10
///
/// Returns None if no MD tag is present.
fn build_alignment_sequence(read: &bam::Record) -> Result<Vec<u8>> {
    // Get MD tag
    let md_tag = if let Ok(Aux::String(md)) = read.aux(b"MD") {
        md.to_owned()
    } else {
        bail!("MD tag not present")
    };

    // Get read sequence, taking into account soft-clipping
    let start = query_start(read);
    let end = query_end(read);
    let seq = read.seq().as_bytes();
    let read_sequence = if start < seq.len() && end <= seq.len() {
        &seq[start..end]
    } else {
        bail!("Invalid sequence coordinates for read");
    };

    // Get CIGAR operations
    let cigar = read.cigar();
    if cigar.is_empty() {
        bail!("No CIGAR string available");
    }

    // Calculate max length for sequence allocation
    let alignment_length = get_alignment_length(read);
    if alignment_length == 0 {
        bail!("Could not determine alignment length");
    }

    // Prepare the result sequence
    let mut result = Vec::with_capacity(alignment_length);

    // Process CIGAR operations to build the sequence
    let mut r_idx = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for _i in 0..*len {
                    if r_idx < read_sequence.len() {
                        result.push(read_sequence[r_idx]);
                        r_idx += 1;
                    }
                }
            }
            Cigar::Del(len) => {
                for _ in 0..*len {
                    result.push(b'-');
                }
            }
            Cigar::RefSkip(_) => {
                // Skip operation, nothing to do
            }
            Cigar::Ins(len) | Cigar::Pad(len) => {
                for _i in 0..*len {
                    if r_idx < read_sequence.len() {
                        // Encode insertions into reference as lowercase
                        let base = read_sequence[r_idx];
                        let lowercase_base = if base >= b'A' && base <= b'Z' {
                            base + 32
                        } else {
                            base
                        };
                        result.push(lowercase_base);
                        r_idx += 1;
                    }
                }
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                // No action needed for clipped bases
            }
        }
    }

    // Process MD tag
    let mut md_idx = 0;
    let mut s_idx = 0;
    let mut nmatches = 0;
    let md_bytes = md_tag.as_bytes();

    // Verify MD tag length matches CIGAR length
    let md_len = get_md_reference_length(md_bytes);
    let insertions = result.iter().filter(|&&b| b.is_ascii_lowercase()).count();

    if md_len + insertions > alignment_length {
        bail!(
            "Invalid MD tag: MD length {} mismatch with CIGAR length {} and {} insertions",
            md_len,
            alignment_length,
            insertions
        );
    }

    while md_idx < md_bytes.len() {
        if md_bytes[md_idx].is_ascii_digit() {
            nmatches = nmatches * 10 + (md_bytes[md_idx] - b'0') as usize;
            md_idx += 1;
            continue;
        } else {
            // Save matches up to this point, skipping insertions
            for _ in 0..nmatches {
                while s_idx < result.len() && result[s_idx].is_ascii_lowercase() {
                    s_idx += 1;
                }
                if s_idx < result.len() {
                    s_idx += 1;
                }
            }
            while s_idx < result.len() && result[s_idx].is_ascii_lowercase() {
                s_idx += 1;
            }

            r_idx += nmatches;
            nmatches = 0;

            if md_idx < md_bytes.len() && md_bytes[md_idx] == b'^' {
                md_idx += 1;
                while md_idx < md_bytes.len() && md_bytes[md_idx].is_ascii_uppercase() {
                    if s_idx < result.len() {
                        result[s_idx] = md_bytes[md_idx];
                        s_idx += 1;
                    }
                    md_idx += 1;
                }
            } else if md_idx < md_bytes.len() {
                // Save mismatch as lowercase
                if s_idx < result.len() {
                    let c = md_bytes[md_idx];
                    let lowercase_c = if c.is_ascii_uppercase() { c + 32 } else { c };
                    result[s_idx] = lowercase_c;
                    s_idx += 1;
                    r_idx += 1;
                    md_idx += 1;
                }
            }
        }
    }

    // Save remaining matches, skipping insertions
    for _ in 0..nmatches {
        while s_idx < result.len() && result[s_idx].is_ascii_lowercase() {
            s_idx += 1;
        }
        if s_idx < result.len() {
            s_idx += 1;
        }
    }
    while s_idx < result.len() && result[s_idx].is_ascii_lowercase() {
        s_idx += 1;
    }

    Ok(result)
}

/// Helper function to get query start position
fn query_start(read: &bam::Record) -> usize {
    let mut start_offset = 0;

    for op in read.cigar().iter() {
        match op {
            Cigar::HardClip(_) => {
                // Only validate if we've already seen other operations
                if start_offset != 0 && start_offset != read.seq_len() {
                    panic!("Invalid clipping in CIGAR string");
                }
            }
            Cigar::SoftClip(len) => {
                start_offset += *len as usize;
            }
            _ => break,
        }
    }

    start_offset
}

/// Helper function to get query end position
fn query_end(read: &bam::Record) -> usize {
    let cigar = read.cigar();
    let end_offset = read.seq_len();

    // If there is no sequence, compute length from cigar string
    if end_offset == 0 {
        let mut offset = 0;
        for op in cigar.iter() {
            match op {
                Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len)
                | Cigar::SoftClip(len)
                    if offset == 0 =>
                {
                    offset += *len as usize;
                }
                _ => {}
            }
        }
        offset
    } else {
        // Walk backwards in cigar string
        let mut offset = end_offset;
        for op in cigar.iter().rev() {
            match op {
                Cigar::HardClip(_) => {
                    if offset != end_offset {
                        panic!("Invalid clipping in CIGAR string");
                    }
                }
                Cigar::SoftClip(len) => {
                    offset -= *len as usize;
                }
                _ => break,
            }
        }
        offset
    }
}

/// Helper function to get alignment length
fn get_alignment_length(read: &bam::Record) -> usize {
    let mut length = 0;

    for op in read.cigar().iter() {
        match op {
            Cigar::SoftClip(_) | Cigar::HardClip(_) => continue,
            _ => length += op.len() as usize,
        }
    }

    length
}

/// Helper function to calculate reference length from MD tag
fn get_md_reference_length(md_tag: &[u8]) -> usize {
    let mut length = 0;
    let mut md_idx = 0;
    let mut nmatches = 0;

    while md_idx < md_tag.len() {
        if md_tag[md_idx].is_ascii_digit() {
            nmatches = nmatches * 10 + (md_tag[md_idx] - b'0') as usize;
            md_idx += 1;
        } else {
            length += nmatches;
            nmatches = 0;

            if md_tag[md_idx] == b'^' {
                md_idx += 1;
                while md_idx < md_tag.len() && md_tag[md_idx].is_ascii_uppercase() {
                    md_idx += 1;
                    length += 1;
                }
            } else {
                md_idx += 1;
                length += 1;
            }
        }
    }

    length + nmatches
}

/// Return the reference sequence in the region that is covered by the
/// alignment of the read to the reference.
///
/// This method requires the MD tag to be set.
///
/// Returns None if MD tag is not present.
fn build_reference_sequence(read: &bam::Record) -> Result<Vec<u8>> {
    // Get the expanded alignment sequence using MD tag
    let ref_seq = build_alignment_sequence(read)?;

    // Get CIGAR operations
    let cigar = read.cigar();
    if cigar.is_empty() {
        bail!("No CIGAR string available");
    }

    // Pre-allocate result with capacity equal to ref_seq
    let mut result = Vec::with_capacity(ref_seq.len());

    // Process CIGAR operations to build the reference sequence
    let mut r_idx = 0; // Index for ref_seq

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                // For matches and mismatches, copy from alignment sequence
                for _ in 0..*len {
                    if r_idx < ref_seq.len() {
                        result.push(ref_seq[r_idx]);
                        r_idx += 1;
                    }
                }
            }
            Cigar::Del(len) => {
                // For deletions, also copy from alignment sequence
                for _ in 0..*len {
                    if r_idx < ref_seq.len() {
                        result.push(ref_seq[r_idx]);
                        r_idx += 1;
                    }
                }
            }
            Cigar::RefSkip(_) => {
                // Skip operation - no action needed
            }
            Cigar::Ins(len) | Cigar::Pad(len) => {
                // For insertions, skip these bases in the alignment sequence
                // as they're not part of the reference
                r_idx += *len as usize;
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                // No action needed for clips
            }
        }
    }

    Ok(result)
}

/// Retrieve the reference sequence for a given read using either
/// the provided reference sequence or by reconstructing it from
/// the MD tag.
///
/// This is a convenience function that handles both methods of obtaining
/// the reference sequence.
fn get_reference_sequence(read: &bam::Record, ref_seq: &[u8]) -> Result<Vec<u8>> {
    if !ref_seq.is_empty() {
        // Use provided reference sequence
        return Ok(ref_seq.to_vec());
    } else {
        // Use MD tag to reconstruct reference sequence
        build_reference_sequence(read)
    }
}

fn process_read_pair(pair: ReadPair, bam: &bam::Reader, trim: i32) -> Result<Option<String>> {
    let read1 = pair.read1;
    let read2 = pair.read2;

    // Get read tags
    let xr = read1
        .aux(b"XR")
        .ok()
        .and_then(|a| a.string())
        .map(|s| s.to_owned())
        .unwrap_or_default();

    let xg = read1
        .aux(b"XG")
        .ok()
        .and_then(|a| a.string())
        .map(|s| s.to_owned())
        .unwrap_or_default();

    let mut ref_seq = Vec::new();
    bam.fetch_seq(
        read1.tid(),
        read1.pos() as u64,
        read1.cigar().end_pos() as u64,
        &mut ref_seq,
    )
    .context("Failed to fetch reference sequence")?;

    let read_seq = read1.seq().as_bytes();

    // Use the built-in reference_positions_full method to get reference positions for each base in the read
    let ref_pos = read1.reference_positions_full();

    // Get the reference sequence using our new function
    let ref_sequence =
        get_reference_sequence(&read1, &ref_seq).context("Failed to get reference sequence")?;

    let (methylation, ga_pos, ct_pos) =
        call_methyl(read_seq, &ref_sequence, &ref_pos, read1.pos() as usize, &xg);

    // Format output
    let output = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        read1.qname().context("Invalid read name")?.escape_ascii(),
        read1.pos(),
        read1.cigar().to_string(),
        std::str::from_utf8(&read_seq).context("Invalid sequence")?,
        methylation,
        ga_pos
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(","),
        ct_pos
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );

    Ok(Some(output))
}

#[derive(Clone)]
struct MethylProcessor {
    bam_path: PathBuf,
    map_quality_filter: u8,
    trimming_tn5: i32,
}

#[derive(Debug, Clone, Serialize)]
struct MethylOutput {
    read_name: String,
    position: i64,
    cigar: String,
    sequence: String,
    methylation: String,
    ga_positions: String,
    ct_positions: String,
}

impl RegionProcessor for MethylProcessor {
    type P = MethylOutput;

    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        // Create a reader
        let mut reader =
            bam::IndexedReader::from_path(&self.reads).expect("Indexed Reader for region");

        // If passed add ref_fasta
        if let Some(ref_fasta) = &self.ref_fasta {
            reader.set_reference(ref_fasta).expect("Set ref");
        }

        // fetch the region of interest
        reader.fetch((tid, start, stop)).expect("Fetched a region");
        let mut results = Vec::new();
        let mut mates = HashMap::new();

        // Walk over each read, counting the starts and ends
        for record in reader
            .records()
            .map(|r| r.expect("Read record"))
            .filter(|read| self.read_filter.filter_read(&read, None))
        {
            let rec_start = u32::try_from(record.reference_start()).expect("check overflow");
            let mate_start = u32::try_from(record.mpos()).expect("check overflow");

            if rec_start < start && mate_start < start {
                continue;
            }

            if !mates.contains_key(record.qname()) {
                mates.insert(record.qname(), record);
                continue;
            }

            let mate = mates.remove(record.qname()).unwrap();

            let (read1, read2) = if record.is_first_in_template() {
                (record, mate)
            } else {
                (mate, record)
            };

            let pair = ReadPair::new(read1, read2);
            if let Ok(Some(output_str)) = process_read_pair(pair, &bam, self.trimming_tn5) {
                // Parse the tab-separated output string
                let fields: Vec<&str> = output_str.trim().split('\t').collect();
                if fields.len() >= 7 {
                    results.push(MethylOutput {
                        read_name: fields[0].to_string(),
                        position: fields[1].parse().unwrap_or(0),
                        cigar: fields[2].to_string(),
                        sequence: fields[3].to_string(),
                        methylation: fields[4].to_string(),
                        ga_positions: fields[5].to_string(),
                        ct_positions: fields[6].to_string(),
                    });
                }
            }
        }

        results
    }
}

pub fn run(args: MethylExtractArgs) -> Result<()> {
    // Set number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.cpu)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Create a processor for parallel processing
    let processor = MethylProcessor {
        bam_path: args.alignment.clone(),
        map_quality_filter: args.map_quality_filter,
        trimming_tn5: args.trimming_tn5,
    };

    // Create a chrom sizes file from the BAM header
    let mut bam = bam::Reader::from_path(&args.alignment).context("Failed to open BAM file")?;

    let header = bam.header();
    let temp_chrom_sizes = PathBuf::from("temp_chrom_sizes.txt");
    {
        let mut chrom_sizes_file = File::create(&temp_chrom_sizes)
            .context("Failed to create temporary chrom sizes file")?;

        for (tid, target) in header.target_names().iter().enumerate() {
            let target_name = std::str::from_utf8(target).unwrap_or("unknown");
            let target_len = header.target_len(tid as u32).unwrap_or(0);
            writeln!(chrom_sizes_file, "{}\t{}", target_name, target_len)
                .context("Failed to write to chrom sizes file")?;
        }
    }

    // Create ParGranges for parallel processing
    let par_granges = ParGranges::new(
        temp_chrom_sizes,
        None,
        false,
        Some(args.cpu),
        None,
        None,
        processor,
    );

    // Process in parallel
    let receiver = par_granges
        .process()
        .context("Failed to process in parallel")?;

    // Write results to output file
    let mut output = File::create(&args.output).context("Failed to create output file")?;

    // Write header
    output
        .write_all(
            b"read_name\tposition\tcigar\tsequence\tmethylation\tga_positions\tct_positions\n",
        )
        .context("Failed to write header")?;

    // Write results
    for result in receiver {
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            result.read_name,
            result.position,
            result.cigar,
            result.sequence,
            result.methylation,
            result.ga_positions,
            result.ct_positions
        )
        .context("Failed to write result to output file")?;
    }

    // Clean up temporary file
    std::fs::remove_file(temp_chrom_sizes).ok();

    Ok(())
}
