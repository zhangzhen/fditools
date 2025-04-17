use anyhow::Result;
use crossbeam_channel::{bounded, Receiver};
use lazy_static::lazy_static;
use log::*;
use rayon::prelude::*;
use serde::Serialize;
use std::thread;

const CHUNKSIZE: usize = 1_000;
const NUM_CHUNKS: usize = 100;

lazy_static! {
    /// CHANNEL_SIZE_MODIFIER as a str
    pub static ref CHUNKSIZE_STR: String = CHUNKSIZE.to_string();

    /// CHUNKSIZE as a str
    pub static ref NUM_CHUNKS_STR: String = NUM_CHUNKS.to_string();
}

pub trait ChunkProcessor {
    type P: 'static + Send + Sync + Serialize;

    fn process_chunk(&self, chunk: Vec<Self::P>) -> Vec<Self::P>;
}

#[derive(Debug)]
pub struct ParChunks<R: 'static + ChunkProcessor + Send + Sync> {
    /// Number of threads this is allowed to use, uses all if None
    threads: usize,
    /// The ideal number of basepairs each worker will receive. Total bp in memory at one time = `threads` * `chunksize`
    chunksize: usize,
    chunk_rx: Receiver<Vec<R::P>>,
    num_chunks: usize,
    /// A modifier to apply to the channel size formular that is (BYTES_INA_GIGABYTE * channel_size_modifier) * threads / size_of(R::P)
    pool: rayon::ThreadPool,
    /// The implementation of [RegionProcessor] that will be used to process regions
    processor: R,
}

impl<R: 'static + ChunkProcessor + Send + Sync> ParChunks<R> {
    pub fn new(
        threads: Option<usize>,
        chunksize: Option<usize>,
        chunk_rx: Receiver<Vec<R::P>>,
        num_chunks: Option<usize>,
        processor: R,
    ) -> Self {
        let threads = if let Some(threads) = threads {
            threads
        } else {
            num_cpus::get()
        };

        // Keep two around for main thread and thread running the pool
        let threads = std::cmp::max(threads.checked_sub(2).unwrap_or(0), 1);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();

        info!("Using {} worker threads.", threads);
        Self {
            threads,
            chunk_rx,
            chunksize: chunksize.unwrap_or(CHUNKSIZE),
            num_chunks: num_chunks.unwrap_or(NUM_CHUNKS),
            pool,
            processor,
        }
    }

    pub fn process(self) -> Result<Receiver<R::P>> {
        let channel_size2 = self.chunksize * self.threads;
        let (result_tx, result_rx) = bounded(channel_size2);

        thread::spawn(move || {
            self.pool.install(|| {
                let mut chunks = Vec::with_capacity(self.num_chunks);
                let mut result = vec![];
                let mut cnt = 0;
                while let Ok(chunk) = self.chunk_rx.recv() {
                    chunks.push(chunk);

                    if chunks.len() >= self.num_chunks {
                        // Process all chunks in parallel
                        let (r, _) = rayon::join(
                            || {
                                chunks
                                    .into_par_iter()
                                    .flat_map(|chunk| self.processor.process_chunk(chunk))
                                    .collect()
                            },
                            || {
                                result.into_iter().for_each(|p| {
                                    result_tx.send(p).expect("Sent a serializable to writer")
                                })
                            },
                        );
                        cnt += 1;
                        info!("Processed {} batches", cnt);
                        result = r;
                        chunks = Vec::with_capacity(self.num_chunks);
                    }
                }

                // Process any remaining chunks
                if !chunks.is_empty() {
                    let (r, _) = rayon::join(
                        || {
                            chunks
                                .into_par_iter()
                                .flat_map(|chunk| self.processor.process_chunk(chunk))
                                .collect()
                        },
                        || {
                            result.into_iter().for_each(|p| {
                                result_tx.send(p).expect("Sent a serializable to writer")
                            })
                        },
                    );
                    cnt += 1;
                    info!("Processed {} batches", cnt);
                    result = r;
                }

                // Send final set of results
                result
                    .into_iter()
                    .for_each(|p| result_tx.send(p).expect("Sent a serializable to writer"));
            });
        });

        Ok(result_rx)
    }
}
