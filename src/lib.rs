#![warn(missing_docs)]

//! The `phylotree` crate aims to be useful when dealing with phylogenetic trees.
//! It can be used to build such trees or read then from newick files. this crate
//! can also be used to compare trees.  
//!
//! Since phylogenetic trees and phylolgenetic distance matrices are so closely related
//! this crate can also be used to extract such matrices from phylolgenetic trees as
//! well as read and write phylip distance matrix files.  
//!
//! # A note on implementation
//! Recursive data structures can be a pain in rust, which is why this crate exists:  
//!   
//! **so you don't have to implement it...**  
//!   
//! To avoid this problem here the tree is stored as a vector
//! of nodes, each node has an identifier and accessing and mutating
//! nodes in a tree is done using these identifiers. As such we can have a
//! non-recursive data structure representing the tree but also easily
//! implement recursive algorithms *(like simple tree traversals)* on this tree.
//!
//! # Using `phylotree`
//! Most of the functionality is implemented in [`crate::tree`]. The
//! [`crate::distance`] module is used to dealt with phylolgenetic distance matrices.
//! [`crate::distr`] is a helper module to provide different branch
//! length distributions when generating random phylogenetic trees.
//!
//! ## Building trees
//! The simplest way to build a tree is to create an empty tree, add a root node and
//! then add children to the various added nodes:
//!
//! ```
//! use phylotree::tree::{Tree, Node};
//!
//! let mut tree = Tree::new();
//!
//! // Add the root node
//! let root = tree.add(Node::new());
//!
//! // Add a child to the root
//! let child1 = tree.add_child(Node::new_named("Child_1"), root, None).unwrap();
//! // Add a child to the root with a branch length
//! let child2 = tree.add_child(Node::new_named("Child_2"), root, Some(0.5)).unwrap();
//!
//! // Add more children
//! let child3 = tree.add_child(Node::new_named("Child_3"), child1, None).unwrap();
//!
//! // Get depth of child
//! assert_eq!(tree.get(&child3).unwrap().get_depth(), 2)
//! ```
//!
//! ## Reading and writing trees
//! This library can build trees strings (or files) encoded in the
//! [newick](https://en.wikipedia.org/wiki/Newick_format) format:
//! ```
//! use phylotree::tree::Tree;
//!
//! let newick_str = "((A:0.1,B:0.2)F:0.6,(C:0.3,D:0.4)E:0.5)G;";
//! let tree = Tree::from_newick(newick_str).unwrap();
//!
//! assert_eq!(tree.to_newick().unwrap(), newick_str)
//! ```
//!
//! ## Traversing trees
//! Several traversals are implemented to visit nodes in a particular order. pre-order,
//! post-order and level-order traversals are implemented on all trees. In-order traversls
//! are implemented only for binary trees. A traversals returns a [`Vec`] of [`tree::NodeId`]
//! in the order they are to be visited in.
//! ```
//! use phylotree::tree::{Tree, Node};
//!
//! //          |
//! //     -----G-----
//! //    |          |
//! // ---C---    ---F---
//! // |     |    |     |
//! // A     B    D     E
//!
//! let newick_str = "((A,B)C,(D,E)F)G;";
//! let mut tree = Tree::from_newick(newick_str).unwrap();
//! let root = tree.get_root().unwrap();
//!
//! let preorder: Vec<_> = tree.preorder(&root).unwrap()
//!     .iter()
//!     .map(|node_id| tree.get(node_id).unwrap().name.clone().unwrap())
//!     .collect();
//!
//! assert_eq!(preorder, vec!["G", "C", "A", "B", "F", "D", "E"]);
//!
//! // Add a child node to F so the tree is no longer binary
//! let f_idx = tree.get_by_name("F").unwrap().id;
//! tree.add_child(Node::new_named("third_child"), f_idx, None).unwrap();
//!
//! assert!(tree.inorder(&root).is_err())
//! ```
//!
//!
//! ## Comparing trees
//! A number of metrics taking into account topology and branch lenghts are implemented
//! in order to compare trees with each other:
//! ```
//! use phylotree::tree::Tree;
//!
//! // The second tree is just a random rotation of the first,
//! // they represent the same phylogeney
//! let newick_orig = "((A:0.1,B:0.2)F:0.6,(C:0.3,D:0.4)E:0.5)G;";
//! let newick_rota = "((D:0.3,C:0.4)E:0.5,(B:0.2,A:0.1)F:0.6)G;";
//!
//! let tree_orig = Tree::from_newick(newick_orig).unwrap();
//! let tree_rota = Tree::from_newick(newick_rota).unwrap();
//!
//! let rf = tree_orig.robinson_foulds(&tree_rota).unwrap();
//!
//! assert_eq!(rf, 0)
//! ```
//!
//! ## Computing distances between nodes in a tree
//! We can get the distance (either as number of edges or sum of edge lengths) betweem
//! nodes of the tree as well as compute the whole phyhlogenetic distance matrix
//! of a tree.
//! ```
//! use phylotree::tree::Tree;
//!
//! // The following tree is encoded by the newick string:
//! //          |
//! //     +----+----+
//! //     |         |
//! //    0.3        |
//! //     |         |
//! //     |        0.6
//! //   --+--       |
//! //   |   |       |
//! //  0.2 0.2      |
//! //   |   |    ---+---
//! //   T3  T1   |     |
//! //            |     |
//! //           0.4   0.5
//! //            |     |
//! //            |     |
//! //            T2    |
//! //                  T0
//!
//! let newick = "((T3:0.2,T1:0.2):0.3,(T2:0.4,T0:0.5):0.6);";
//! let tree = Tree::from_newick(newick).unwrap();
//!
//! let t0 = tree.get_by_name("T0").unwrap();
//! let t3 = tree.get_by_name("T3").unwrap();
//!
//! let (edge_sum, num_edges) = tree.get_distance(&t0.id, &t3.id).unwrap();
//!
//! assert_eq!(num_edges, 4);
//! assert_eq!(edge_sum, Some(0.5 + 0.6 + 0.3 + 0.2));
//!
//! // Compute the whole distance matrix
//! let matrix = tree.distance_matrix_recursive().unwrap();
//! let phylip="\
//! 4
//! T0    0  1.6  0.9  1.6
//! T1    1.6  0  1.5  0.4
//! T2    0.9  1.5  0  1.5
//! T3    1.6  0.4  1.5  0
//! ";
//!
//! assert_eq!(matrix.to_phylip(true).unwrap(), phylip)
//! ```
//!
use std::collections::VecDeque;

use clap::ValueEnum;
use distr::{Distr, Sampler};
use rand::prelude::*;

use tree::{Node, Tree, TreeError};

#[cfg(feature = "python")]
pub mod python;

pub mod distance;
pub mod distr;
pub mod tree;

// type Error = Box<dyn std::error::Error>;
// type Result<T> = std::result::Result<T, Error>;

/// Shape of random trees to generate
#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum TreeShape {
    /// Yule model tree shape
    Yule,
    /// Caterpillar tree shape
    Caterpillar,
    /// Ete3 Tree.populate replicate
    Ete3,
}

/// Genereates a random binary tree of a given size.
pub fn generate_tree(
    n_leaves: usize,
    brlens: bool,
    sampler_type: Distr,
) -> Result<Tree, TreeError> {
    let mut tree = Tree::new();
    // Add root
    tree.add(Node::default());

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

/// Generate a random binary tree under the Yule model.
pub fn generate_yule(
    n_leaves: usize,
    brlens: bool,
    sampler_type: Distr,
) -> Result<Tree, TreeError> {
    // Initialize tree
    let mut tree = Tree::new();
    let root = tree.add(Node::default());

    let mut rng = thread_rng();
    let sampler = Sampler::new(sampler_type);

    let mut parent_candidates = vec![root];

    while tree.n_leaves() != n_leaves {
        // Choose parent
        let parent = parent_candidates
            .choose(&mut rng)
            .expect("No parent candidate")
            .clone();

        // Generate child node
        let edge1: Option<f64> = brlens.then_some(sampler.sample(&mut rng));
        let edge2: Option<f64> = brlens.then_some(sampler.sample(&mut rng));
        let child1 = tree.add_child(Node::default(), parent, edge1)?;
        let child2 = tree.add_child(Node::default(), parent, edge2)?;
        parent_candidates.push(child1);
        parent_candidates.push(child2);

        let pos = parent_candidates.iter().position(|n| *n == parent).unwrap();
        parent_candidates.swap_remove(pos);
    }

    // Assign names to tips
    for (i, tip_idx) in tree.get_leaves().iter().cloned().enumerate() {
        tree.get_mut(&tip_idx)?.set_name(format!("Tip_{i}"));
    }

    Ok(tree)
}

/// Generates a caterpillar tree by adding children to the last node addesd to the tree
/// until we reach the desired numebr of leaves.
pub fn generate_caterpillar(
    n_leaves: usize,
    brlens: bool,
    sampler_type: Distr,
) -> Result<Tree, TreeError> {
    let mut tree = Tree::new();
    tree.add(Node::default());

    let mut rng = thread_rng();
    let sampler = Sampler::new(sampler_type);

    let mut parent = 0;
    for i in 1..n_leaves {
        let parent_bkp = parent;

        let l1: Option<f64> = brlens.then_some(sampler.sample(&mut rng));
        let l2: Option<f64> = brlens.then_some(sampler.sample(&mut rng));

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
