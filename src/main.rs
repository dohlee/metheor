use clap::{Parser};



mod pdr;
mod pm;
mod me;
mod lpmd;
mod bamutil;
mod readutil;

fn main() {
    let args = metheor::Cli::parse();

    match &args.command {
        metheor::Commands::Pdr { input, output, min_depth, min_cpgs } => {
            pdr::compute(input, output, *min_depth, *min_cpgs);
        }
        metheor::Commands::Pm { input, output, min_depth } => {
            pm::compute(input, output, *min_depth);
        }
        metheor::Commands::Me { input, output, min_depth } => {
            me::compute(input, output, *min_depth);
        }
        metheor::Commands::Lpmd { input, output, min_distance, max_distance, cpg_set, min_qual } => {
            lpmd::compute(input, output, *min_distance, *max_distance, cpg_set, *min_qual);
        }
    }
}
