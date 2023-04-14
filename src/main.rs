use phylotree::*;
use clap::Parser;

mod cli;

fn main() {
    let args = cli::Args::parse();

    let random = generate_tree(args.size, args.branch_lengths);
    random.to_file(&args.output).unwrap()
}
