use clap::Parser;
use phylotree::*;

mod cli;

fn main() {
    match cli::Args::parse().command {
        cli::Commands::Generate {
            tips,
            branch_lengths,
            output,
        } => {
            let random = generate_tree(tips, branch_lengths);
            random.to_file(&output).unwrap()
        }
        cli::Commands::Stats { tree } => {
            let mut tree = Tree::from_file(&tree).unwrap();
            println!("Height: {:?}", tree.height());
            println!("Number of nodes: {:?}", tree.size());
            println!("Number of tips: {:?}", tree.get_leaves().len());
            println!("Is rooted ? {}", tree.is_rooted());
            println!("Is binary ? {}", tree.is_binary());
            println!("Sackin: {:?}", tree.sackin());
        }
    }
}
