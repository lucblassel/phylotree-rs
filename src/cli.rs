use clap::{Parser, Subcommand};
use std::path::PathBuf;

use crate::distr::Distr;

/// A simple command line tool to generate a random phylogenetic tree
#[derive(Parser, Debug)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}


#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Generate random tree(s)
    Generate {
        /// Number of tips in the generated tree
        #[arg(short, long, default_value_t = 20)]
        tips: usize,

        /// Generate uniform branch lengths
        #[arg(short, long)]
        branch_lengths: bool,

        /// Output file (directory if generating multiple trees)
        #[arg(short, long)]
        output: PathBuf,

        /// Number of trees to generate
        #[arg(short='n', long)]
        trees: Option<usize>,

        /// Distribution of branch lengths
        #[arg(value_enum, short, long, default_value_t=Distr::Uniform)]
        distribution: Distr
    },

    /// Get statistics about a tree
    Stats {
        /// Input newick file of the tree
        trees: Vec<PathBuf>
    },

    /// Compare two trees
    Compare {
        /// Reference tree
        reftree: PathBuf,
        /// Tree to compare to reference
        tocompare: PathBuf,
    },
    /// Compute RF between two tree
    RF {
        /// Reference tree
        reftree: PathBuf,
        /// Tree to compare to reference
        tocompare: PathBuf,
    }

}