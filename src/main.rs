use clap::{Parser};

mod pdr;
mod pm;
mod me;
mod fdrp;
mod lpmd;
mod bamutil;
mod readutil;
mod progressbar;

fn main() {
    let args = metheor::Cli::parse();

    match &args.command {
        metheor::Commands::Pdr { input, output, min_depth, min_cpgs, min_qual, cpg_set } => {
            pdr::compute(input, output, *min_depth, *min_cpgs, *min_qual, cpg_set);
        }
        metheor::Commands::Pm { input, output, min_depth, min_qual, cpg_set } => {
            pm::compute(input, output, *min_depth, *min_qual, cpg_set);
        }
        metheor::Commands::Me { input, output, min_depth, min_qual, cpg_set } => {
            me::compute(input, output, *min_depth, *min_qual, cpg_set);
        }
        metheor::Commands::Fdrp { input, output, min_qual, max_depth, min_overlap, cpg_set } => {
            fdrp::compute(input, output, *min_qual, *max_depth, *min_overlap, cpg_set);
        }
        metheor::Commands::Lpmd { input, output, pairs, min_distance, max_distance, min_qual, cpg_set } => {
            lpmd::compute(input, output, *min_distance, *max_distance, *min_qual, cpg_set, pairs );
        }
    }
}
