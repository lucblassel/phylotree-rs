//! Build and manioulate phylogenetic trees.
//!
//! This module defines the essential structs used to represent phylogenetic trees:
//!  - The [`Node`] struct that represents a node of a phylogenetic tree.
//!  - The [`Tree`] struct that holds a collection of [`Node`] objects.
//!  - The [`TreeList`] struct that holds an ordered collection of [`Tree`] objects.
//!

mod bipartition;
/// A module to draw phylogenetic trees
pub mod draw;
mod newick;
mod node;
mod taxon_index;
mod tree_impl;
mod tree_list;

pub use self::bipartition::{BipartitionProfile, RobinsonFouldsMode};
pub use self::newick::{
    InternalNodeLabelMode, NewickBranchLengthFormat, NewickCommentFormat, NewickFormat,
    NewickNameFormat, NewickParseError, NewickParseOptions, NewickSerializeOptions,
    NewickSupportFormat, NewickToken, Span, SupportCommentMode, TokenizerError,
};
pub use self::node::{Node, NodeError};
pub use self::taxon_index::TaxonIndex;
pub use self::tree_impl::{Comparison, Tree, TreeError};
pub use self::tree_list::{Bipartition, TopologyFrequency, TreeList, TreeListError};

/// A type that represents Identifiers of [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type NodeId = usize;

/// A type that represents branch lengths between [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type EdgeLength = f64;

/// A type that represents the depth (i.e. distance from the root) o
/// given edge within a phylogenetic [`Tree`] object.
pub type EdgeDepth = usize;
