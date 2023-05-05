//! This module defines the two essential structs to represent phylogenetic trees:  
//!  - The [`Node`] struct that represents a node of a phylogenetic tree.
//!  - The [`Tree`] struct that holds a collection of [`Node`] objects. 

mod node;
pub use self::node::{Node, NodeError};

mod tree;
pub use self::tree::{Tree, TreeError, ParseError};

pub type NodeId = usize;
pub type Edge = f64;