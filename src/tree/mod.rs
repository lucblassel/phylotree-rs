//! Build and manipulate phylogenetic trees.
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
mod traversal;
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
pub use self::traversal::{
    InOrderTreeTraversal, LevelOrderTreeTraversal, NodeTraversal, NodeTraversalMut,
    PostOrderTreeTraversal, PreOrderTreeTraversal,
};
pub use self::tree_impl::{Comparison, Tree, TreeError};
pub use self::tree_list::{Bipartition, TopologyFrequency, TreeList, TreeListError};

/// An identifier for a [`Node`] within a phylogenetic [`Tree`].
///
/// Node IDs are local to the tree that created them. Use [`NodeId::index`] when
/// an underlying storage index is required, or the [`From`] conversions to move
/// between `NodeId` and [`usize`].
///
/// ```
/// use phylotree::tree::NodeId;
///
/// let id = NodeId::from(3);
/// assert_eq!(id.index(), 3);
/// assert_eq!(usize::from(id), 3);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(usize);

impl NodeId {
    /// Returns the underlying storage index.
    ///
    /// ```
    /// use phylotree::tree::NodeId;
    /// assert_eq!(NodeId::from(3).index(), 3);
    /// ```
    #[inline]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl From<usize> for NodeId {
    #[inline]
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<NodeId> for usize {
    #[inline]
    fn from(value: NodeId) -> Self {
        value.index()
    }
}

impl std::fmt::Display for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.index().fmt(f)
    }
}
/// A type that represents branch lengths between [`Node`] objects
/// within phylogenetic [`Tree`] object.
pub type EdgeLength = f64;

/// A type that represents the depth (i.e. distance from the root) o
/// given edge within a phylogenetic [`Tree`] object.
pub type EdgeDepth = usize;
