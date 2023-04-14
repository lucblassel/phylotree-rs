use clap::Parser;
use std::path::PathBuf;

/// A simple command line tool to generate a random Phylogenentic tree
#[derive(Parser, Debug)]
pub struct Args {
    /// Number of leaves in the tree
    #[arg(short, long, default_value_t = 20)]
    pub size: usize,

    /// Generate uniform branch lengths
    #[arg(short, long)]
    pub branch_lengths: bool,

    /// File to save the tree to
    #[arg(short, long)]
    pub output: PathBuf,
}