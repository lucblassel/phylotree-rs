#![warn(missing_docs)]

//! Phylogenetic trees, topology comparisons, and distance matrices.
//!
//! - [`tree`]: parse Newick, traverse or edit trees, and compare topologies.
//! - [`simulation`]: seeded topology and branch-length priors.
//! - [`distance`]: compact distance matrices, PHYLIP I/O, and UPGMA.
//! - [`distr`]: fixed branch-length distributions for the legacy generators.
//!
//! Trees store nodes in a vector rather than a recursive data structure. Access
//! and modify nodes through tree-local [`tree::NodeId`] values. Use
//! [`tree::TreeList`] to read, write, and summarize collections of trees.
//!
//! # Building trees
//!
//! Start with an empty tree, add a root, then attach children. Both names and
//! branch lengths are optional; [`tree::Tree::add_child`] maintains the links
//! between parent and child nodes.
//!
//! ```
//! use phylotree::tree::{Node, Tree};
//!
//! let mut tree = Tree::new();
//! let root = tree.add(Node::new_named("root"));
//! let a = tree.add_child(Node::new_named("A"), root, Some(0.5)).unwrap();
//! tree.add_child(Node::new_named("B"), root, Some(0.2)).unwrap();
//!
//! assert_eq!(tree.get(a).unwrap().get_depth(), 1);
//! assert_eq!(tree.n_leaves(), 2);
//! ```
//!
//! # Reading, writing, and traversing trees
//!
//! Parse a Newick string with [`tree::Tree::from_newick`] and serialize it with
//! [`tree::Tree::to_newick`]. File-based equivalents are also available.
//! Traversals return node IDs in visiting order: preorder, postorder, and
//! level order work on general trees; inorder requires a binary tree.
//!
//! ```
//! use phylotree::tree::Tree;
//!
//! //          |
//! //     -----G-----
//! //     |         |
//! //  ---C---   ---F---
//! //  |     |   |     |
//! //  A     B   D     E
//! let newick = "((A,B)C,(D,E)F)G;";
//! let tree = Tree::from_newick(newick).unwrap();
//! let root = tree.root().unwrap();
//! let names: Vec<_> = tree.preorder(root).unwrap().into_iter()
//!     .map(|id| tree.get(id).unwrap().name.as_deref().unwrap())
//!     .collect();
//!
//! assert_eq!(names, vec!["G", "C", "A", "B", "F", "D", "E"]);
//! assert_eq!(tree.to_newick().unwrap(), newick);
//! ```
//!
//! # Comparing trees
//!
//! Robinson–Foulds (RF) distance compares groupings of tips. Choose rooted
//! comparisons for descendant clades or unrooted comparisons for splits.
//! Rotating children does not change the topology:
//!
//! ```
//! use phylotree::tree::Tree;
//!
//! let first = Tree::from_newick("((A,B),(C,D));").unwrap();
//! let rotated = Tree::from_newick("((D,C),(B,A));").unwrap();
//! assert_eq!(first.unrooted_robinson_foulds(&rotated).unwrap(), 0);
//! ```
//!
//! [`tree::Tree::compare_topologies`] also computes normalized RF, weighted RF,
//! and the Kuhner–Felsenstein branch score. Retain a [`tree::BipartitionProfile`]
//! when comparing one tree against many others.
//!
//! # Distances and distance matrices
//!
//! [`tree::Tree::get_distance`] returns both the branch-length sum (when all
//! lengths on the path are present) and the number of edges. A distance matrix
//! stores all pairwise tip distances and supports PHYLIP input and output.
//!
//! ```
//! use phylotree::tree::Tree;
//!
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
//! let tree = Tree::from_newick("((T3:0.2,T1:0.2):0.3,(T2:0.4,T0:0.5):0.6);").unwrap();
//! let t0 = tree.get_by_name("T0").unwrap().id;
//! let t3 = tree.get_by_name("T3").unwrap().id;
//! let (length, edges) = tree.get_distance(t0, t3).unwrap();
//!
//! assert_eq!(edges, 4);
//! assert!((length.unwrap() - 1.6).abs() < 1e-12);
//!
//! let matrix = tree.distance_matrix().unwrap();
//! assert!((*matrix.get("T0", "T3").unwrap() - 1.6).abs() < 1e-12);
//! let phylip = matrix.to_phylip(true).unwrap();
//! assert!(phylip.starts_with("4\n"));
//! ```
//!
//! To obtain just the longest tip-to-tip distance, use [`tree::Tree::diameter`];
//! it runs in linear time without constructing a distance matrix.
//!
//! # Simulating trees
//!
//! Supply your own RNG for reproducible simulations. Choose a topology prior
//! independently of its branch-length prior:
//!
//! ```
//! use phylotree::simulation::{simulate, BranchLengths, Topology};
//! use rand::{rngs::StdRng, SeedableRng};
//!
//! let mut rng = StdRng::seed_from_u64(42);
//! let mut tree = simulate(20, Topology::UniformUnrooted,
//!     BranchLengths::Gamma { shape: 2.0, scale: 0.05 }, &mut rng).unwrap();
//! tree.rescale_diameter(1.0).unwrap();
//! assert!((tree.diameter().unwrap() - 1.0).abs() < 1e-12);
//! ```
//!
//! See [`simulation`] for the uniform, ETE3, Yule, and birth-death priors.
//!
use std::collections::VecDeque;

use clap::ValueEnum;
use distr::{Distr, Sampler};
use rand::prelude::*;

use tree::{Node, NodeId, Tree, TreeError};

#[cfg(feature = "python")]
pub mod python;

pub mod distance;
pub mod distr;
pub mod simulation;
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

/// Generate an ETE3-style binary topology using the thread-local RNG.
///
/// Requires at least three tips. `brlens` controls whether lengths are sampled
/// from `sampler_type`; see [`simulation::simulate`] for seedable priors.
///
/// ```
/// use phylotree::{generate_tree, distr::Distr};
/// let tree = generate_tree(10, true, Distr::Exponential).unwrap();
/// assert_eq!(tree.n_leaves(), 10);
/// ```
pub fn generate_tree(
    n_leaves: usize,
    brlens: bool,
    sampler_type: Distr,
) -> Result<Tree, TreeError> {
    if n_leaves < 3 {
        return Err(TreeError::InvalidTreeSize(n_leaves));
    }

    let mut tree = Tree::new();
    // Add root
    tree.add(Node::default());

    let mut rng = thread_rng();

    let sampler = Sampler::new(sampler_type);

    let mut next_deq: VecDeque<NodeId> = VecDeque::new();
    next_deq.push_back(0.into());

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

    for (i, &id) in next_deq.iter().enumerate() {
        tree.get_mut(id)?.set_name(format!("Tip_{i}"));
    }

    Ok(tree)
}

/// Generate a Yule binary topology using the thread-local RNG.
///
/// Requires at least three tips; `brlens` enables the selected branch prior.
/// See [`simulation::simulate`] for caller-owned RNG state.
///
/// ```
/// use phylotree::{generate_yule, distr::Distr};
/// let tree = generate_yule(10, true, Distr::Exponential).unwrap();
/// assert_eq!(tree.n_leaves(), 10);
/// ```
pub fn generate_yule(
    n_leaves: usize,
    brlens: bool,
    sampler_type: Distr,
) -> Result<Tree, TreeError> {
    if n_leaves < 3 {
        return Err(TreeError::InvalidTreeSize(n_leaves));
    }

    // Initialize tree
    let mut tree = Tree::new();
    let root = tree.add(Node::default());

    let mut rng = thread_rng();
    let sampler = Sampler::new(sampler_type);

    let mut parent_candidates = vec![root];

    while tree.n_leaves() != n_leaves {
        // Choose parent
        let parent = *parent_candidates
            .choose(&mut rng)
            .expect("No parent candidate");
        // .clone();

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
    for (i, tip_idx) in tree.leaves().iter().cloned().enumerate() {
        tree.get_mut(tip_idx)?.set_name(format!("Tip_{i}"));
    }

    Ok(tree)
}

/// Generate a caterpillar topology with at least three tips.
///
/// `brlens` enables lengths sampled from `sampler_type` using the thread-local RNG.
///
/// ```
/// use phylotree::{generate_caterpillar, distr::Distr};
/// let tree = generate_caterpillar(10, true, Distr::Exponential).unwrap();
/// assert_eq!(tree.n_leaves(), 10);
/// ```
pub fn generate_caterpillar(
    n_leaves: usize,
    brlens: bool,
    sampler_type: Distr,
) -> Result<Tree, TreeError> {
    if n_leaves < 3 {
        return Err(TreeError::InvalidTreeSize(n_leaves));
    }

    let mut tree = Tree::new();
    tree.add(Node::default());

    let mut rng = thread_rng();
    let sampler = Sampler::new(sampler_type);

    let mut parent: NodeId = 0.into();
    for i in 0..n_leaves - 1 {
        let parent_bkp = parent;

        let l1: Option<f64> = brlens.then_some(sampler.sample(&mut rng));
        let l2: Option<f64> = brlens.then_some(sampler.sample(&mut rng));

        if i == n_leaves - 2 {
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

#[cfg(test)]
mod tests {
    use crate::{
        distr::Distr,
        generate_caterpillar, generate_tree, generate_yule,
        tree::{Tree, TreeError},
    };
    use std::collections::HashSet;

    #[test]
    fn test_invalid_sizes() {
        for size in [0, 1, 2] {
            matches!(
                generate_tree(size, false, Distr::Uniform),
                Err(TreeError::InvalidTreeSize(_)),
            );
            matches!(
                generate_yule(size, false, Distr::Uniform),
                Err(TreeError::InvalidTreeSize(_)),
            );
            matches!(
                generate_caterpillar(size, false, Distr::Uniform),
                Err(TreeError::InvalidTreeSize(_)),
            );
        }
    }

    #[test]
    fn test_requested_size() {
        for size in [10, 100, 3, 3000] {
            let tree = generate_tree(size, false, Distr::Uniform).unwrap();
            assert_eq!(size, tree.n_leaves());

            let tree = generate_yule(size, false, Distr::Uniform).unwrap();
            assert_eq!(size, tree.n_leaves());

            let tree = generate_caterpillar(size, false, Distr::Uniform).unwrap();
            assert_eq!(size, tree.n_leaves());
        }
    }

    fn n_unique_names(tree: &Tree) -> usize {
        let names: HashSet<String> =
            HashSet::from_iter(tree.leaf_names().into_iter().map(|v| v.unwrap()));
        names.len()
    }

    #[test]
    fn test_unique_tip_names() {
        for size in [10, 100, 3, 3000] {
            let tree = generate_tree(size, false, Distr::Uniform).unwrap();
            assert_eq!(size, n_unique_names(&tree));

            let tree = generate_yule(size, false, Distr::Uniform).unwrap();
            assert_eq!(size, n_unique_names(&tree));

            let tree = generate_caterpillar(size, false, Distr::Uniform).unwrap();
            assert_eq!(size, n_unique_names(&tree));
        }
    }
}
