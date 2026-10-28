use clap::Parser;

mod bamutil;
mod fdrp;
mod lpmd;
mod me;
mod mhl;
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
        } => {
            pdr::compute(input, output, *min_depth, *min_cpgs, *min_qual, cpg_set);
        }
        metheor::Commands::Pm {
            input,
            output,
            min_depth,
            min_qual,
            cpg_set,
            reference_cpgs,
        } => {
            pm::compute(
                input,
                output,
                *min_depth,
                *min_qual,
                cpg_set,
                reference_cpgs,
            );
        }
        metheor::Commands::Me {
            input,
            output,
            min_depth,
            min_qual,
            cpg_set,
            reference_cpgs,
        } => {
            me::compute(
                input,
                output,
                *min_depth,
                *min_qual,
                cpg_set,
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
        } => {
            fdrp::compute(
                input,
                output,
                *min_qual,
                *min_depth,
                *max_depth,
                *min_overlap,
                cpg_set,
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
            qfdrp::compute(
                input,
                output,
                *min_qual,
                *min_depth,
                *max_depth,
                *min_overlap,
                cpg_set,
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
