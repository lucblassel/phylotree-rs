use std::collections::VecDeque;

use distr::{Distr, Sampler};
use ptree::{print_tree, TreeBuilder};
use rand::prelude::*;

use tree::{Node, Tree};

pub mod distance;
pub mod distr;
// pub mod errors;
pub mod tree;

type Error = Box<dyn std::error::Error>;
type Result<T> = std::result::Result<T, Error>;


/// Genereates a random binary tree of a given size. Branch lengths are uniformly distributed
pub fn generate_tree(n_leaves: usize, brlens: bool, sampler_type: Distr) -> Result<Tree> {
    let mut tree = Tree::new();
    let mut rng = thread_rng();

    let sampler = Sampler::new(sampler_type);

    let mut next_deq = VecDeque::new();
    next_deq.push_back(0);

    for _ in 0..(n_leaves - 1) {
        let parent_id = if rng.gen_bool(0.5) {
            next_deq.pop_front()
        } else {
            next_deq.pop_back()
        }
        .unwrap();
        let l1: Option<f64> = if brlens {
            Some(sampler.sample(&mut rng))
        } else {
            None
        };
        let l2: Option<f64> = if brlens {
            Some(sampler.sample(&mut rng))
        } else {
            None
        };
        next_deq.push_back(tree.add_child(Node::new(), parent_id, l1)?);
        next_deq.push_back(tree.add_child(Node::new(), parent_id, l2)?);
    }

    for (i, id) in next_deq.iter().enumerate() {
        tree.get_mut(id)?.set_name(format!("Tip_{i}"));
    }

    Ok(tree)
}


/// Generates a caterpillar tree by adding children to the last node addesd to the tree
/// until we reach the desired numebr of leaves. Branch lengths are uniformly distributed
pub fn generate_caterpillar(n_leaves: usize, brlens: bool) -> Result<Tree> {
    let mut tree = Tree::new();
    let mut rng = thread_rng();

    let mut parent = 0;
    for i in 1..n_leaves {
        let parent_bkp = parent;
        let l1: Option<f64> = if brlens { Some(rng.gen()) } else { None };
        let l2: Option<f64> = if brlens { Some(rng.gen()) } else { None };
        if i == n_leaves - 1 {
            // Adding tip
            tree.add_child(Node::new_named(&format!("Tip_{i}")), parent, l1)?;
            tree.add_child(Node::new_named(&format!("Tip_{}", i + 1)), parent, l2)?;
        } else {
            // Adding parent node
            parent = tree.add_child(Node::new(), parent, l1)?;
            tree.add_child(Node::new_named(&format!("Tip_{i}")), parent_bkp, l2)?;
        }
    }

    Ok(tree)
}
