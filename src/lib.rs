use std::str;

use clap::{AppSettings, Parser, Subcommand};

pub mod readutil;

/// Summarizes the heterogeneity of DNA methylation states using BAM files.
#[derive(Parser)]
#[clap(name = "metheor")]
#[clap(about = "Summarizes the heterogeneity of DNA methylation states using BAM files.")]
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
        #[clap(long, short='i', required=true)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true)]
        output: String,

        /// Minimum depth of CpG stretches to consider.
        #[clap(long, short='d', default_value_t=10)]
        min_depth: i32,

        /// Minimum number of consecutive CpGs in a CpG stretch to consider.
        #[clap(long, short='c', default_value_t=10)]
        min_cpgs: i32,

    },
    /// Compute epipolymorphism.
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    Pm {
        /// Input BAM file.
        #[clap(long, short='i', required=true)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true)]
        output: String,

        /// Minimum depth of CpG quartets to consider
        #[clap(long, short='d', default_value_t=10)]
        min_depth: i32,
    },
    /// Compute methylation entropy.
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    Me {
        /// Input BAM file.
        #[clap(long, short='i', required=true)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true)]
        output: String,

        /// Minimum depth of CpG quartets to consider
        #[clap(long, short='d', default_value_t=10)]
        min_depth: i32,
    },
    Lpmd {
        /// Input BAM file.
        #[clap(long, short='i', required=true)]
        input: String,

        /// Output csv file.
        #[clap(long, short='o', required=true)]
        output: String,

        /// Minimum depth of CpG quartets to consider
        #[clap(long, short='m', default_value_t=2)]
        min_distance: i32,

        /// Minimum depth of CpG quartets to consider
        #[clap(long, short='M', default_value_t=16)]
        max_distance: i32,

        /// Specify a set of CpGs (in BED file) to be analyzed.
        #[clap(long, short='c', default_value="")]
        cpg_set: String,
        
        /// Minimum quality for a read to be considered.
        #[clap(long, short='c', default_value_t=10)]
        min_qual: u8,
    }
}
