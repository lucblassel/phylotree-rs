//! Build and manioulate phylogenetic trees.  
//!  
//! This module defines the two essential structs to represent phylogenetic trees:  
//!  - The [`Node`] struct that represents a node of a phylogenetic tree.
//!  - The [`Tree`] struct that holds a collection of [`Node`] objects.
//!

/// A module to draw phylogenetic trees
pub mod draw;
mod node;
mod tree_impl;

pub use self::node::{Node, NodeError};
pub use self::tree_impl::{Comparison, NewickParseError, Tree, TreeError};

/// A type that represents Identifiers of [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type NodeId = usize;

/// A type that represents branch lengths between [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type EdgeLength = f64;

/// A type that represents the depth (i.e. distance from the root) o        
/// given edge within a phylogenetic [`Tree`] object.
pub type EdgeDepth = usize;

/// Newick output format
#[derive(Debug, Copy, Clone)]
pub enum NewickFormat {
    /// Output all supported and available fields
    AllFields,
    /// Only output topology
    Topology,
    /// Output all fields except for comments
    NoComments,
    /// Output node names
    OnlyNames,
    /// Output branch lengths
    OnlyLengths,
    /// Output leaf branch lenghs + all node names
    LeafLengthsAllNames,
    /// Output leaf branch lengths + leaf node names
    LeafLengthsLeafNames,
    /// Output internal branch lenghts + leaf node names
    InternalLengthsLeafNames,
    /// Output all branch lenghts + leaf names
    AllLengthsLeafNames,
}
