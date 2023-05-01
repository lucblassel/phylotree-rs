use clap::Parser;
use phylotree::*;
use std::path::Path;

mod cli;

fn print_header() {
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
                std::fs::create_dir_all(&output).unwrap();

                for i in 1..=ntrees {
                    let output = output.join(format!("{i}_{tips}_tips.nwk"));
                    let random = generate_tree(tips, branch_lengths, distribution);
                    random.to_file(&output).unwrap()
                }
            } else {
                let random = generate_tree(tips, branch_lengths, distribution);
                random.to_file(&output).unwrap()
            }
        }
        cli::Commands::Stats { trees } => {
            print_header();
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

            dbg!("Ref: {}", reftree.leaf_index);
            dbg!("Cmp: {}", compare.leaf_index);
            eprintln!("Ref: {ref_parts:#?}");
            eprintln!("Cmp: {other_parts:#?}");

            println!("tree\treference\tcommon\tcompared\trf");
            println!(
                "0\t{}\t{}\t{}\t{}",
                ref_parts.len() - common,
                common,
                other_parts.len() - common,
                rf,
            )
        }
        cli::Commands::RF { reftree, tocompare } => {
            let reftree = Tree::from_file(&reftree).unwrap();
            let compare = Tree::from_file(&tocompare).unwrap();

            let rf = reftree.robinson_foulds_new(&compare).unwrap();

            println!("{rf}")
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
