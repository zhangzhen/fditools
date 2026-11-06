pub mod commands;
use anyhow::Result;
use clap::{Parser, Subcommand};
use commands::*;
use env_logger::Env;
use fditools_lib::utils;
use log::*;

#[derive(Parser)]
#[command(version, author, about, rename_all = "kebab-case")]
/// Commands for generating per-base analysis
struct Args {
    #[command(subcommand)]
    subcommand: Subcmd,
}

#[derive(Subcommand)]
enum Subcmd {
    CellBarcode(cell_barcode::CellBarcode),
}

impl Subcmd {
    fn run(self) -> Result<()> {
        match self {
            Subcmd::CellBarcode(x) => x.run()?,
        }
        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    if let Err(err) = Args::parse().subcommand.run() {
        if utils::is_broken_pipe(&err) {
            std::process::exit(0);
        }
        error!("{}", err);
        std::process::exit(1);
    }
    Ok(())
}
