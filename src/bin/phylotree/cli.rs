use clap::{Parser, Subcommand};
use std::path::PathBuf;

use phylotree::distr::Distr;

/// A simple command line tool to manipulate phylogenetic trees
#[derive(Parser, Debug)]
pub struct Args {
    #[command(subcommand)]
    /// The command to execute
    pub command: Commands,
}

/// The available commands in the `phylotree` tool
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
        output: Option<PathBuf>,

        /// Number of trees to generate
        #[arg(short = 'n', long)]
        trees: Option<usize>,

        /// Distribution of branch lengths
        #[arg(value_enum, short, long, default_value_t=Distr::Uniform)]
        distribution: Distr,
    },

    /// Get statistics about a tree
    Stats {
        /// Input newick file of the tree
        trees: Vec<PathBuf>,
    },

    /// Compare two trees
    Compare {
        /// Reference tree
        reftree: PathBuf,
        /// Tree to compare to reference
        tocompare: PathBuf,
    },
    /// Output the phylogenetic distance matrix of the tree
    Matrix {
        /// The phylogenetic tree
        tree: PathBuf,
        /// Output a square matrix instead of a triangular one
        #[arg(short, long)]
        square: bool,
        /// FIle to save the matrix to
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Collapse branches that are under a certain branch length threshold
    Collapse {
        /// The phylogenetic tree
        tree: PathBuf,
        /// The threshold under which to collapse the branch
        threshold: phylotree::tree::Edge,
        /// Print the number of collapsed branches at the end
        #[arg(short, long)]
        verbose: bool,
        /// Don't collapse external branches
        #[arg(short, long)]
        exclude_tips: bool,
        /// File to save the tree to
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

