use std::str;

use clap::{Parser, Subcommand};

pub mod readutil;
pub mod bamutil;
pub mod progressbar;
pub mod lpmd;

/// Summarizes the heterogeneity of DNA methylation states using BAM files.
#[derive(Parser)]
#[clap(name = "metheor")]
#[clap(about = "Summarizes the heterogeneity of DNA methylation states using BAM files.")]
#[clap(version = "0.0.8")]
#[clap(author = "Dohoon Lee. <dohlee.bioinfo@gmail.com>\nBonil Koo. <bikoo95@snu.ac.kr>")]
#[clap(arg_required_else_help = true)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Compute proportion of discordant reads (PDR).
    #[clap(arg_required_else_help = true)]
    Pdr {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Path to output table file summarizing the result of PM calculation.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum depth of CpG stretches to consider.
        #[clap(long, short='d', default_value_t=10, display_order=3)]
        min_depth: u32,

        /// Minimum number of consecutive CpGs in a CpG stretch to consider.
        #[clap(long, short='c', default_value_t=10, display_order=4)]
        min_cpgs: usize,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=5)]
        min_qual: u8,

        /// (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', required=false, display_order=6)]
        cpg_set: Option<String>,
    },
    /// Compute epipolymorphism.
    #[clap(arg_required_else_help = true)]
    Pm {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Path to output table file summarizing the result of ME calculation.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum depth of CpG quartets to consider
        #[clap(long, short='d', default_value_t=10, display_order=3)]
        min_depth: u32,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=4)]
        min_qual: u8,

        /// (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', required=false, display_order=5)]
        cpg_set: Option<String>,
    },
    /// Compute methylation entropy.
    #[clap(arg_required_else_help = true)]
    Me {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Path to output table file summarizing the result of PDR calculation.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum depth of CpG quartets to consider.
        #[clap(long, short='d', default_value_t=10, display_order=3)]
        min_depth: u32,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=4)]
        min_qual: u8,

        /// (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', required=false, display_order=5)]
        cpg_set: Option<String>,
    },
    /// Compute fraction of discordant read pairs (FDRP).
    #[clap(arg_required_else_help = true)]
    Fdrp {
        /// Path to input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Path to output table file summarizing the result of FDRP calculation.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=3)]
        min_qual: u8,

        /// Maximum number of reads to consider.
        #[clap(long, short='n', default_value_t=40, display_order=4)]
        max_depth: usize,

        /// Minimum overlap between two reads to consider in bp.
        #[clap(long, short='l', default_value_t=35, display_order=5)]
        min_overlap: i32,

        /// (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', required=false, display_order=6)]
        cpg_set: Option<String>,
    },
    /// Compute quantitative fraction of discordant read pairs (qFDRP).
    #[clap(arg_required_else_help = true)]
    Qfdrp {
        /// Path to input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Path to output table file summarizing the result of FDRP calculation.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=3)]
        min_qual: u8,

        /// Maximum number of reads to consider.
        #[clap(long, short='n', default_value_t=40, display_order=4)]
        max_depth: usize,

        /// Minimum overlap between two reads to consider in bp.
        #[clap(long, short='l', default_value_t=35, display_order=5)]
        min_overlap: i32,

        /// (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', required=false, display_order=6)]
        cpg_set: Option<String>,
    },
    /// Compute local pairwise methylation discordance (LPMD).
    #[clap(arg_required_else_help = true)]
    Lpmd {
        /// Path to input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Path to output table file summarizing the result of LPMD calculation.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// (Optional) Concordance information for all CpG pairs.
        #[clap(long, short='p', required=false, display_order=3)]
        pairs: Option<String>,

        /// Minimum distance between CpG pairs to consider.
        #[clap(long, short='m', default_value_t=2, display_order=4)]
        min_distance: i32,

        /// Maximum distance between CpG pairs to consider.
        #[clap(long, short='M', default_value_t=16, display_order=5)]
        max_distance: i32,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=6)]
        min_qual: u8,

        /// (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', required=false, display_order=7)]
        cpg_set: Option<String>,
    }
}
