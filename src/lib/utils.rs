use anyhow::{Error, Result};
use bio::io::fastq;
use flate2::read::MultiGzDecoder;
use grep_cli::stdout;
use gzp::{deflate::Gzip, Compression, ZBuilder};
use lazy_static::lazy_static;
use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    path::Path,
};
use termcolor::ColorChoice;

/// Check if err is a broken pipe.
#[inline]
pub fn is_broken_pipe(err: &Error) -> bool {
    if let Some(io_err) = err.root_cause().downcast_ref::<io::Error>() {
        if io_err.kind() == io::ErrorKind::BrokenPipe {
            return true;
        }
    }
    false
}

pub fn get_reader<P: AsRef<Path>>(
    path: &Option<P>,
    gzipped: bool,
) -> Result<fastq::Reader<BufReader<Box<dyn Read>>>> {
    let raw_reader: Box<dyn Read> = match &path {
        Some(path) if path.as_ref().to_str().unwrap() != "-" => {
            let file = File::open(path)?;
            if gzipped {
                Box::new(MultiGzDecoder::new(file))
            } else {
                Box::new(BufReader::new(file))
            }
        }
        _ => {
            let reader = std::io::stdin();
            if gzipped {
                Box::new(MultiGzDecoder::new(reader))
            } else {
                Box::new(reader)
            }
        }
    };

    Ok(fastq::Reader::new(raw_reader))
}

pub fn get_writer<P: AsRef<Path>>(
    path: &Option<P>,
    gzipped: bool,
    threads: usize,
    compression_level: u32,
) -> Result<fastq::Writer<Box<dyn Write>>> {
    let raw_writer: Box<dyn Write> = match &path {
        Some(path) if path.as_ref().to_str().unwrap() != "-" => {
            let writer = BufWriter::new(File::create(path)?);
            if gzipped {
                Box::new(
                    ZBuilder::<Gzip, _>::new()
                        .num_threads(threads)
                        .compression_level(Compression::new(compression_level))
                        .from_writer(writer),
                )
            } else {
                Box::new(writer)
            }
        }
        _ => {
            let writer = stdout(ColorChoice::Never);
            if gzipped {
                Box::new(
                    ZBuilder::<Gzip, _>::new()
                        .num_threads(threads)
                        .compression_level(Compression::new(compression_level))
                        .from_writer(writer),
                )
            } else {
                Box::new(writer)
            }
        }
    };
    Ok(fastq::Writer::new(raw_writer))
}

lazy_static! {
    /// Return the number of cpus as an &str
    pub static ref NUM_CPU: String = num_cpus::get().to_string();
}
