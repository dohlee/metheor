use std::str;

use clap::{AppSettings, Parser, Subcommand};

pub mod readutil;
pub mod bamutil;
pub mod progressbar;
pub mod lpmd;

/// Summarizes the heterogeneity of DNA methylation states using BAM files.
#[derive(Parser)]
#[clap(name = "metheor")]
#[clap(about = "Summarizes the heterogeneity of DNA methylation states using BAM files.")]
#[clap(version = "0.0.3")]
#[clap(author = "Dohoon Lee. <dohlee.bioinfo@gmail.com>")]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Compute proportion of discordant reads (PDR).
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    Pdr {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum depth of CpG stretches to consider.
        #[clap(long, short='d', default_value_t=10, display_order=3)]
        min_depth: u32,

        /// Minimum number of consecutive CpGs in a CpG stretch to consider.
        #[clap(long, short='c', default_value_t=10, display_order=4)]
        min_cpgs: i32,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=5)]
        min_qual: u8,
    },
    /// Compute epipolymorphism.
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    Pm {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum depth of CpG quartets to consider
        #[clap(long, short='d', default_value_t=10, display_order=3)]
        min_depth: u32,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=4)]
        min_qual: u8,
    },
    /// Compute methylation entropy.
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    Me {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Minimum depth of CpG quartets to consider.
        #[clap(long, short='d', default_value_t=10, display_order=3)]
        min_depth: u32,

        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=4)]
        min_qual: u8,
    },
    /// Compute local pairwise methylation discordance (LPMD).
    Lpmd {
        /// Input BAM file.
        #[clap(long, short='i', required=true, display_order=1)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true, display_order=2)]
        output: String,

        /// Concordance information for all CpG pairs. (Optional)
        #[clap(long, short='p', display_order=3)]
        pairs: String,

        /// Minimum distance between CpG pairs to consider.
        #[clap(long, short='m', default_value_t=2, display_order=4)]
        min_distance: i32,

        /// Maximum distance between CpG pairs to consider.
        #[clap(long, short='M', default_value_t=16, display_order=5)]
        max_distance: i32,

        /// Specify a set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', default_value="", display_order=6)]
        cpg_set: String,
        
        /// Minimum quality for a read to be considered.
        #[clap(long, short='q', default_value_t=10, display_order=7)]
        min_qual: u8,
    }
}
