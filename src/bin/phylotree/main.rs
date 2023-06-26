#![warn(missing_docs)]
//! The `phylotree` binary is a command line tool, using the `[phylotree]` crate.
//! It is made to execute common operations on phylogenetic trees directly in the terminal.

use clap::Parser;
use itertools::Itertools;
use phylotree::{
    generate_tree,
    tree::{
        draw::{self, Layout, Node},
        Tree,
    },
};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    fs::File,
    io,
    io::{BufWriter, Write},
    path::Path,
};
use tinytemplate::TinyTemplate;

/// contains the struct representing the command line arguments
/// parsed by [`clap`] and used to execute this binary
pub mod cli;

fn print_stats_header() {
    println!("height\tdiameter\tnodes\ttips\trooted\tbinary\tncherries\tcolless\tsackin")
}

fn print_stats(path: &Path) {
    let tree = Tree::from_file(path).unwrap();

    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        tree.height()
            .map_or_else(|_| "-".to_string(), |v| format!("{}", v)),
        tree.diameter().unwrap_or(0.0),
        tree.size(),
        tree.n_leaves(),
        tree.is_rooted().unwrap(),
        tree.is_binary().unwrap(),
        tree.cherries().unwrap(),
        tree.colless()
            .map(|v| format!("{}", v))
            .unwrap_or("-".to_string()),
        tree.sackin()
            .map(|v| format!("{}", v))
            .unwrap_or("-".to_string()),
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

            let stats = reftree.compare_topologies(&compare).unwrap();

            println!("tree\treference\tcommon\tcompared\trf\trf_w\tbranch_score");
            println!(
                "0\t{}\t{}\t{}\t{}\t{}\t{}",
                ref_parts.len() - common,
                common,
                other_parts.len() - common,
                stats.rf,
                stats.weighted_rf,
                stats.branch_score,
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
        cli::Commands::Collapse {
            tree,
            threshold,
            verbose,
            exclude_tips,
            output,
        } => {
            let mut tree = Tree::from_file(&tree).unwrap();
            let mut n = 0;
            for node_idx in tree.preorder(&tree.get_root().unwrap()).unwrap().iter() {
                let node = tree.get_mut(node_idx).unwrap();
                if exclude_tips && node.is_tip() {
                    continue;
                }
                let parent_idx = node.parent;
                let mut collapsed = false;
                if let Some(len) = node.parent_edge {
                    if len < threshold {
                        node.set_parent(parent_idx.unwrap(), Some(0.0));
                        collapsed = true;
                        n += 1;
                        if verbose {
                            eprintln!("Collapsed {node:?}")
                        }
                    }
                }
                if collapsed {
                    let parent = tree.get_mut(&parent_idx.unwrap()).unwrap();
                    parent.set_child_edge(node_idx, Some(0.0))
                }
            }

            if let Some(output) = output {
                tree.to_file(&output).unwrap()
            } else {
                println!("{}", tree.to_newick().unwrap());
            }

            if verbose {
                eprintln!("{n}")
            }
        }
        cli::Commands::Remove { tree, tips, output } => {
            let mut tree = Tree::from_file(&tree).unwrap();
            for tip_name in tips.iter() {
                let node = tree.get_by_name(tip_name).unwrap();
                if !node.is_tip() {
                    panic!("{} is not a tip node!", tip_name)
                }
                let id = node.id;
                tree.prune(&id).unwrap();
            }

            tree.compress().unwrap();

            if let Some(output) = output {
                tree.to_file(&output).unwrap()
            } else {
                println!("{}", tree.to_newick().unwrap())
            }
        }
        cli::Commands::Distance { tree, tips, output } => {
            let tree = Tree::from_file(&tree).unwrap();
            let mut writer = BufWriter::new(match output {
                Some(path) => Box::new(File::create(&path).unwrap()) as Box<dyn Write>,
                None => Box::new(io::stdout()) as Box<dyn Write>,
            });

            writer.write("Seq1\tSeq2\tDistance\n".as_bytes()).unwrap();

            for pair in tips.iter().combinations(2) {
                let (n1, n2) = (pair[0], pair[1]);
                let i1 = tree.get_by_name(n1).unwrap().id;
                let i2 = tree.get_by_name(n2).unwrap().id;
                if let (Some(d), _) = tree.get_distance(&i1, &i2).unwrap() {
                    writer
                        .write(format!("{n1}\t{n2}\t{d}\n").as_bytes())
                        .unwrap();
                } else {
                    panic!("The tree is missing some branch lengths between {n1} and {n2}")
                }
            }
        }
        cli::Commands::Deduplicate {
            tree,
            alignment,
            collapse,
            output,
            verbose,
        } => {
            let mut tree = Tree::from_file(&tree).unwrap();

            // Get duplicated sequence names
            let mut alignment = needletail::parse_fastx_file(alignment).unwrap();
            let mut duplicate_sequences = BTreeMap::new();
            while let Some(r) = alignment.next() {
                let record = r.unwrap();
                let seq = String::from_utf8(record.seq().to_vec()).unwrap();
                let id = String::from_utf8(record.id().to_vec()).unwrap();

                let entry = duplicate_sequences.entry(seq.clone()).or_insert(vec![]);
                (*entry).push(id.trim().to_string());
            }

            let mut n = 0;

            let duplicates: Vec<_> = duplicate_sequences
                .into_iter()
                .filter(|(_, v)| v.len() > 1)
                .map(|(_, mut v)| {
                    v.sort();
                    v
                })
                .collect();

            if duplicates.is_empty() {}

            let mut removed = vec![];

            // Remove/Collapse tips corresponding to duplicate_sequences
            for duplicated in duplicates {
                if duplicated.len() == 1 {
                    continue;
                }
                for (i, taxa) in duplicated.iter().enumerate() {
                    let idx = tree.get_by_name(taxa).unwrap().id;
                    let parent_idx = tree.get_by_name(taxa).unwrap().parent.unwrap();

                    if collapse {
                        removed.push(taxa.clone());
                        tree.get_mut(&idx)
                            .unwrap()
                            .set_parent(parent_idx, Some(0.0));
                        tree.get_mut(&parent_idx)
                            .unwrap()
                            .set_child_edge(&idx, Some(0.0));
                    } else if i == 0 {
                        // If we are deleting duplicates we
                        // still want to keep one instance
                        continue;
                    } else {
                        removed.push(taxa.clone());
                        tree.prune(&idx).unwrap();
                    }
                    n += 1;
                }
            }

            if !collapse && n > 0 {
                tree.compress().unwrap();
            }

            // Check if we removed a full cherry, in which case we
            // Also need to remove the parent branch
            let mut pruned = false;
            while !pruned {
                pruned = true;
                for leaf_idx in tree.get_leaves() {
                    let leaf = tree.get(&leaf_idx).unwrap();
                    let unnamed = leaf.name.is_none();
                    if unnamed {
                        pruned = false;
                        tree.prune(&leaf_idx).unwrap();
                    }
                }
                tree.compress().unwrap();
            }

            if let Some(path) = output {
                tree.to_file(&path).unwrap();
            } else {
                println!("{}", tree.to_newick().unwrap())
            }

            if verbose > 0 {
                eprint!("{}", n);
                if verbose > 1 {
                    eprint!(": {}", removed.join(" "));
                }
                eprintln!()
            }
        }
        cli::Commands::Draw {
            tree,
            width,
            height,
            transparent,
            padding,
            output,
        } => {
            let tree = Tree::from_file(&tree).unwrap();
            let mut layout = draw::radial_layout(&tree).unwrap();

            let mut min = f64::INFINITY;
            let mut max = f64::NEG_INFINITY;
            for Node { x, y, label: _ } in layout.nodes.iter() {
                let i = x.min(*y);
                let a = x.max(*y);

                min = min.min(i);
                max = max.max(a);
            }

            let scale = width.min(height) / (max - min);
            layout.rescale(scale);
            let padding_scale = (100. - padding) / 100.;

            let ctx = Context {
                xmin: -(width / 2.0) as i32,
                ymin: -(height / 2.0) as i32,
                width: width as i32,
                height: height as i32,
                scale: padding_scale,
                transparent,
                layout,
            };

            let mut writer = BufWriter::new(match output {
                Some(path) => Box::new(File::create(&path).unwrap()) as Box<dyn Write>,
                None => Box::new(io::stdout()) as Box<dyn Write>,
            });

            let mut tt = TinyTemplate::new();
            tt.add_template("svg", SVG).unwrap();
            writer
                .write(tt.render("svg", &ctx).unwrap().as_bytes())
                .unwrap();
        }
    }
}

#[derive(Serialize)]
struct Context {
    xmin: i32,
    ymin: i32,
    width: i32,
    height: i32,
    scale: f64,
    transparent: bool,
    layout: Layout,
}

static SVG : &'static str =  "\
<?xml version=\"1.0\" standalone=\"no\"?>
<svg viewBox=\"{xmin} {ymin} {width} {height}\" width='100%' height='100%' xmlns='http://www.w3.org/2000/svg'>
    {{ if not transparent }}
    <rect fill=\"white\" x=\"{xmin}\" y=\"{ymin}\" width=\"100%\" height=\"100%\"/>
    {{ endif }}
    <g transform=\"scale({scale})\">
    {{ for branch in layout.branches }}
        <path stroke=\"black\" stroke-width=\"1\" fill=\"none\" d=\"M {branch.xstart} {branch.ystart} L {branch.xend} {branch.yend}\"/>
    {{ endfor }}
    {{ for node in layout.nodes }}
        {{ if node.label }}
        <text x=\"{node.x}\" y=\"{node.y}\" class=\"small\" font-size=\"2px\">{node.label}</text>
        {{ endif }}
    {{ endfor }}
    </g>
</svg>
";
