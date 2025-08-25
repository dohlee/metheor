use clap::Parser;
use rayon::ThreadPoolBuilder;

mod bamutil;
mod fdrp;
mod lpmd;
mod me;
mod mhl;
mod pdr;
mod pm;
mod progressbar;
mod qfdrp;
mod readutil;
mod tag;

fn main() {
    let args = metheor::Cli::parse();

    // Configure rayon thread pool
    let num_threads = if args.threads == 0 {
        // Auto-detect number of threads
        num_cpus::get()
    } else {
        args.threads
    };

    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to build rayon thread pool");

    match &args.command {
        metheor::Commands::Pdr {
            input,
            output,
            min_depth,
            min_cpgs,
            min_qual,
            cpg_set,
        } => {
            pdr::compute(input, output, *min_depth, *min_cpgs, *min_qual, cpg_set);
        }
        metheor::Commands::Pm {
            input,
            output,
            min_depth,
            min_qual,
            cpg_set,
        } => {
            pm::compute(input, output, *min_depth, *min_qual, cpg_set);
        }
        metheor::Commands::Me {
            input,
            output,
            min_depth,
            min_qual,
            cpg_set,
        } => {
            me::compute(input, output, *min_depth, *min_qual, cpg_set);
        }
        metheor::Commands::Fdrp {
            input,
            output,
            min_qual,
            min_depth,
            max_depth,
            min_overlap,
            cpg_set,
        } => {
            fdrp::compute_with_threshold(
                input,
                output,
                *min_qual,
                *min_depth,
                *max_depth,
                *min_overlap,
                cpg_set,
                args.parallel_threshold,
            );
        }
        metheor::Commands::Qfdrp {
            input,
            output,
            min_qual,
            min_depth,
            max_depth,
            min_overlap,
            cpg_set,
        } => {
            qfdrp::compute_with_threshold(
                input,
                output,
                *min_qual,
                *min_depth,
                *max_depth,
                *min_overlap,
                cpg_set,
                args.parallel_threshold,
            );
        }
        metheor::Commands::Mhl {
            input,
            output,
            min_depth,
            min_cpgs,
            min_qual,
            cpg_set,
        } => {
            mhl::compute(input, output, *min_depth, *min_cpgs, *min_qual, cpg_set);
        }
        metheor::Commands::Lpmd {
            input,
            output,
            pairs,
            min_distance,
            max_distance,
            min_qual,
            cpg_set,
        } => {
            lpmd::compute(
                input,
                output,
                *min_distance,
                *max_distance,
                *min_qual,
                cpg_set,
                pairs,
            );
        }
        metheor::Commands::Tag {
            input,
            output,
            genome,
        } => {
            tag::run(input, output, genome);
        }
    }
}
