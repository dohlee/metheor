use clap::Parser;

mod alignmentwriter;
mod bamutil;
mod fdrp;
mod lpmd;
mod me;
mod mhl;
mod overlap;
mod pairwise;
mod pdr;
mod pm;
mod progressbar;
mod qfdrp;
mod quartetcatalog;
mod readutil;
mod tag;

fn main() {
    let args = metheor::Cli::parse();

    match &args.command {
        metheor::Commands::Pdr {
            input,
            output,
            min_depth,
            min_cpgs,
            min_qual,
            cpg_set,
            ignore_overlaps,
        } => {
            pdr::compute(
                input,
                output,
                *min_depth,
                *min_cpgs,
                *min_qual,
                cpg_set,
                *ignore_overlaps,
            );
        }
        metheor::Commands::Pm {
            input,
            output,
            min_depth,
            min_qual,
            cpg_set,
            ignore_overlaps,
            include_counts,
            reference_cpgs,
        } => {
            pm::compute(
                input,
                output,
                *min_depth,
                *min_qual,
                cpg_set,
                *ignore_overlaps,
                *include_counts,
                reference_cpgs,
            );
        }
        metheor::Commands::Me {
            input,
            output,
            min_depth,
            min_qual,
            cpg_set,
            ignore_overlaps,
            reference_cpgs,
        } => {
            me::compute(
                input,
                output,
                *min_depth,
                *min_qual,
                cpg_set,
                *ignore_overlaps,
                reference_cpgs,
            );
        }
        metheor::Commands::Fdrp {
            input,
            output,
            min_qual,
            min_depth,
            max_depth,
            min_overlap,
            cpg_set,
            ignore_overlaps,
        } => {
            fdrp::compute(
                input,
                output,
                *min_qual,
                *min_depth,
                *max_depth,
                *min_overlap,
                cpg_set,
                *ignore_overlaps,
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
            ignore_overlaps,
        } => {
            qfdrp::compute(
                input,
                output,
                *min_qual,
                *min_depth,
                *max_depth,
                *min_overlap,
                cpg_set,
                *ignore_overlaps,
            );
        }
        metheor::Commands::Mhl {
            input,
            output,
            min_depth,
            min_cpgs,
            min_qual,
            cpg_set,
            ignore_overlaps,
        } => {
            mhl::compute(
                input,
                output,
                *min_depth,
                *min_cpgs,
                *min_qual,
                cpg_set,
                *ignore_overlaps,
            );
        }
        metheor::Commands::Lpmd {
            input,
            output,
            pairs,
            min_distance,
            max_distance,
            min_qual,
            cpg_set,
            ignore_overlaps,
        } => {
            if *ignore_overlaps {
                lpmd::compute_with_overlap_correction(
                    input,
                    output,
                    *min_distance,
                    *max_distance,
                    *min_qual,
                    cpg_set,
                    pairs,
                    true,
                );
            } else {
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
        }
        metheor::Commands::Tag {
            input,
            output,
            genome,
            threads,
        } => {
            tag::run(input, output, genome, threads.get());
        }
    }
}
