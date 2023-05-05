//! Build and manioulate phylogenetic trees.  
//!  
//! This module defines the two essential structs to represent phylogenetic trees:  
//!  - The [`Node`] struct that represents a node of a phylogenetic tree.
//!  - The [`Tree`] struct that holds a collection of [`Node`] objects.
//! 

mod node;
mod tree;

pub use self::node::{Node, NodeError};
pub use self::tree::{NewickParseError, Tree, TreeError, Comparison};

/// A type that represents Identifiers of [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type NodeId = usize;

/// A type that represents branch lengths between [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type Edge = f64;
