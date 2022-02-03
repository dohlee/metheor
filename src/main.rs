use clap::{Parser};

mod pdr;
mod pm;
mod me;
mod lpmd;
mod bamutil;
mod readutil;
mod progressbar;

fn main() {
    let args = metheor::Cli::parse();

    match &args.command {
        metheor::Commands::Pdr { input, output, min_depth, min_cpgs, min_qual } => {
            pdr::compute(input, output, *min_depth, *min_cpgs, *min_qual);
        }
        metheor::Commands::Pm { input, output, min_depth, min_qual } => {
            pm::compute(input, output, *min_depth, *min_qual);
        }
        metheor::Commands::Me { input, output, min_depth } => {
            me::compute(input, output, *min_depth);
        }
        metheor::Commands::Lpmd { input, output, min_distance, max_distance, cpg_set, min_qual } => {
            lpmd::compute(input, output, *min_distance, *max_distance, cpg_set, *min_qual);
        }
    }
}
