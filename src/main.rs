use clap::Parser;
use phylotree::*;
use std::path::PathBuf;

/// A simple command line tool to generate a random Phylogenentic tree
#[derive(Parser, Debug)]
struct Args {
    /// Number of leaves in the tree
    #[arg(short, long, default_value_t = 20)]
    size: usize,

    /// Generate uniform branch lengths
    #[arg(short, long)]
    branch_lengths: bool,

    /// File to save the tree to
    #[arg(short, long)]
    output: PathBuf,
}

fn main() {
    let args = Args::parse();
    println!("{:?}", args);

    let random = generate_tree(args.size, args.branch_lengths);
    random.to_file(&args.output).unwrap()
}
