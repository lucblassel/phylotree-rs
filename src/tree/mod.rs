//! Build and manioulate phylogenetic trees.
//!
//! This module defines the two essential structs to represent phylogenetic trees:
//!  - The [`Node`] struct that represents a node of a phylogenetic tree.
//!  - The [`Tree`] struct that holds a collection of [`Node`] objects.
//!

/// A module to draw phylogenetic trees
pub mod draw;
mod newick_serializer;
mod node;
mod parser;
mod tokenizer;
mod tree_impl;

pub use self::newick_serializer::NewickFormat;
pub use self::node::{Node, NodeError};
pub use self::parser::NewickParseError;
pub use self::tokenizer::{NewickToken, Span, TokenizerError};
pub use self::tree_impl::{Comparison, Tree, TreeError};

/// A type that represents Identifiers of [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type NodeId = usize;

/// A type that represents branch lengths between [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type EdgeLength = f64;

/// A type that represents the depth (i.e. distance from the root) o
/// given edge within a phylogenetic [`Tree`] object.
pub type EdgeDepth = usize;
