#![warn(missing_docs)]
//! The `phylotree` binary is a command line tool, using the `[phylotree]` crate.
//! It is made to execute common operations on phylogenetic trees directly in the terminal.

use clap::Parser;
use phylotree::{generate_tree, tree::Tree};
use std::path::Path;

/// contains the struct representing the command line arguments
/// parsed by [`clap`] and used to execute this binary
pub mod cli;

fn print_stats_header() {
    println!("height\tnodes\ttips\trooted\tbinary\tncherries\tcolless\tsackin")
}

fn print_stats(path: &Path) {
    let tree = Tree::from_file(path).unwrap();
    println!(
        "{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
        tree.height(),
        tree.size(),
        tree.get_leaves().len(),
        tree.is_rooted(),
        tree.is_binary(),
        tree.cherries(),
        tree.colless(),
        tree.sackin()
    )
}
fn main() {
    match cli::Args::parse().command {
        cli::Commands::Generate {
            tips,
            branch_lengths,
            output,
            trees,
            distribution,
        } => {
            if let Some(ntrees) = trees {
                // Create output directory if it's missing
                assert!(
                    output.is_some(),
                    "If you are generating multiple trees you must specify an output directory"
                );
                let output = output.unwrap();
                std::fs::create_dir_all(&output).unwrap();

                for i in 1..=ntrees {
                    let output = output.join(format!("{i}_{tips}_tips.nwk"));
                    let random = generate_tree(tips, branch_lengths, distribution).unwrap();
                    random.to_file(&output).unwrap()
                }
            } else {
                let random = generate_tree(tips, branch_lengths, distribution).unwrap();
                if let Some(output) = output {
                    random.to_file(&output).unwrap()
                } else {
                    println!("{}", random.to_newick().unwrap())
                }
            }
        }
        cli::Commands::Stats { trees } => {
            print_stats_header();
            for tree in trees {
                print_stats(&tree)
            }
        }
        cli::Commands::Compare { reftree, tocompare } => {
            let reftree = Tree::from_file(&reftree).unwrap();
            let compare = Tree::from_file(&tocompare).unwrap();

            let ref_parts = reftree.get_partitions().unwrap();
            let other_parts = compare.get_partitions().unwrap();

            let common = ref_parts.intersection(&other_parts).count();
            let rf = reftree.robinson_foulds(&compare).unwrap();

            println!("tree\treference\tcommon\tcompared\trf");
            println!(
                "0\t{}\t{}\t{}\t{}",
                ref_parts.len() - common,
                common,
                other_parts.len() - common,
                rf,
            )
        }
        cli::Commands::Matrix {
            tree,
            square,
            output,
        } => {
            let tree = Tree::from_file(&tree).unwrap();
            let dm = tree.distance_matrix().unwrap();
            if let Some(output) = output {
                dm.to_file(&output, square).unwrap();
            } else {
                println!("{}", dm.to_phylip(square).unwrap())
            }
        }
    }
}
