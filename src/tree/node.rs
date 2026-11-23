use std::{
    cell::RefCell,
    collections::HashMap,
    fmt::{Debug, Display},
};

use thiserror::Error;
use vec_map::VecMap;

use super::{EdgeLength, NewickFormat, NodeId};

/// Errors that can occur when manipulating [`Node`] structs.
#[derive(Error, Debug)]
pub enum NodeError {
    /// We are trying to access the an unexisting child of the node
    #[error("Node {parent} does not have child {child}.")]
    HasNoChild {
        /// Id of the parent the parent node
        parent: NodeId,
        /// Id of the inexistant child node
        child: NodeId,
    },
    /// We are trying to access the parent of a parentless node
    #[error("Node {0} does not have a parent")]
    HasNoParent(NodeId),
    /// We are trying to read an edge length that is missing
    #[error("Missing edge length between nodes {parent} and {child}")]
    MissingEdgeLength {
        /// Id of the parent the parent node
        parent: NodeId,
        /// Id of the inexistant child node
        child: NodeId,
    },
    /// We are trying to add the same child twice
    #[error("Node already has child {0}.")]
    DuplicateChild(NodeId),
    /// Trying to add a non-finite support value
    #[error("Trying to set non-finite support value {support} to node {node}.")]
    NonFiniteSupport {
        /// Id of the current node
        node: NodeId,
        /// Problematic support value
        support: f64,
    },
}

#[derive(Clone)]
/// A node of the Tree
pub struct Node {
    /// Index of the node
    pub(crate) id: NodeId,
    /// Name of the node
    pub(crate) name: Option<String>,
    /// Index of the parent node
    pub(crate) parent: Option<NodeId>,
    /// Indices of child nodes
    pub(crate) children: Vec<NodeId>,
    /// Length of the branch between this node and its parent.
    pub(crate) parent_edge: Option<EdgeLength>,
    /// Support value associated with the branch leading to this node.
    ///
    /// The value is stored without normalization because both fractions and
    /// percentages are common in phylogenetic formats.
    pub(crate) support: Option<f64>,
    /// Optional comment attached to node
    pub(crate) comment: Option<String>,
    /// lenght of branches between node and children
    pub(crate) child_edges: Option<HashMap<NodeId, EdgeLength>>,
    /// Distance to descendants of this node
    pub(crate) subtree_distances: RefCell<Option<VecMap<EdgeLength>>>,
    /// Number of edges to root
    pub(crate) depth: usize,
    // Whether the node is deleted or not
    pub(crate) deleted: bool,
}

impl Node {
    /// Creates a new Node
    ///
    /// ```
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// assert!(node.is_tip());
    /// assert!(node.is_root());
    /// ```
    pub fn new() -> Self {
        Self {
            id: NodeId(0),
            name: None,
            parent: None,
            children: vec![],
            parent_edge: None,
            support: None,
            child_edges: None,
            subtree_distances: RefCell::new(None),
            comment: None,
            depth: 0,
            deleted: false,
        }
    }

    /// Creates a new named Node
    ///
    /// ```
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// let node = Node::new_named("A");
    /// assert_eq!(node.name(), Some("A"));
    /// ```
    pub fn new_named(name: &str) -> Self {
        Self {
            id: NodeId(0),
            name: Some(String::from(name)),
            parent: None,
            children: vec![],
            parent_edge: None,
            support: None,
            child_edges: None,
            subtree_distances: RefCell::new(None),
            comment: None,
            depth: 0,
            deleted: false,
        }
    }

    /// Get Id of the node
    pub fn id(&self) -> NodeId {
        self.id
    }

    /// Get name of node if it exists
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Get the id of the nodes parent if it exists
    pub fn parent(&self) -> Option<NodeId> {
        self.parent
    }

    /// Get length of parent edge if it exists
    pub fn parent_edge(&self) -> Option<EdgeLength> {
        self.parent_edge
    }

    /// Get the ids of the node's children
    pub fn children(&self) -> &[NodeId] {
        &self.children[..]
    }

    /// Get the Node's comment
    pub fn comment(&self) -> Option<&str> {
        self.comment.as_deref()
    }

    /// Get the Node's parent branch support value
    pub fn support(&self) -> Option<f64> {
        self.support
    }

    /// Gets the depth of the node
    pub fn depth(&self) -> usize {
        self.depth
    }

    /// Sets the Node's name
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = Some(name.into());
    }

    /// Clears the Node's name
    pub fn clear_name(&mut self) {
        self.name = None;
    }

    /// Set a comment value on the current Node
    pub fn set_comment(&mut self, comment: impl Into<String>) {
        self.comment = Some(comment.into());
    }

    /// Clears the Node's comment
    pub fn clear_comment(&mut self) {
        self.comment = None;
    }

    /// Sets the support value for the parent branch of this node
    pub fn set_support(&mut self, support: f64) -> Result<(), NodeError> {
        if !support.is_finite() {
            return Err(NodeError::NonFiniteSupport {
                node: self.id,
                support,
            });
        }

        self.support = Some(support);
        Ok(())
    }

    /// Clears the support value for this Node's parent branch
    /// This function will not fail if there is no parent branch.
    pub fn clear_support(&mut self) {
        self.support = None;
    }

    /// Sets the internal Node id
    pub(crate) fn set_id(&mut self, id: NodeId) {
        self.id = id;
    }

    /// Set this node's parent and incoming length.
    ///
    /// Update the parent's children separately; prefer [`super::Tree::add_child`]
    /// when editing a tree.
    pub(crate) fn set_parent(&mut self, parent: NodeId, parent_edge: Option<EdgeLength>) {
        self.parent = Some(parent);
        self.parent_edge = parent_edge;
    }

    /// Sets the depth of the node
    pub(crate) fn set_depth(&mut self, depth: usize) {
        self.depth = depth;
    }

    /// Empties the node and sets it as deleted
    pub(crate) fn delete(&mut self) {
        let id = self.id;
        *self = Self::new();
        self.id = id;
        self.deleted = true;
    }

    /// Adds a child to the node
    /// ```
    /// use phylotree::tree::{Node, Tree};
    /// let mut tree = Tree::new();
    /// let root = tree.add_root(Node::new()).unwrap();
    /// let child = tree.add_child(Node::new(), root, Some(0.1)).unwrap();
    /// assert_eq!(tree.get(root).unwrap().children(), &[child]);
    /// assert_eq!(tree.get(root).unwrap().get_child_edge(child), Some(0.1));
    /// assert_eq!(tree.get(child).unwrap().parent_edge(), Some(0.1));
    /// ```
    pub(crate) fn add_child(
        &mut self,
        child: NodeId,
        child_edge: Option<EdgeLength>,
    ) -> Result<(), NodeError> {
        if self.children.contains(&child) {
            return Err(NodeError::DuplicateChild(child));
        }
        self.children.push(child);
        self.set_child_edge(child, child_edge);

        Ok(())
    }

    /// Get the Edge between node and a child
    pub fn get_child_edge(&self, child: NodeId) -> Option<EdgeLength> {
        if let Some(edges) = &self.child_edges {
            edges.get(&child).copied()
        } else {
            None
        }
    }

    /// Sets the Edge between a node and its child
    ///
    /// ```
    /// use phylotree::tree::{Node, Tree};
    /// let mut tree = Tree::new();
    /// let root = tree.add_root(Node::new()).unwrap();
    /// let child = tree.add_child(Node::new(), root, Some(0.5)).unwrap();
    /// tree.set_edge(root, child, Some(1.0)).unwrap();
    /// assert_eq!(tree.get(root).unwrap().get_child_edge(child), Some(1.0));
    /// assert_eq!(tree.get(child).unwrap().parent_edge(), Some(1.0));
    /// // Clearing the edge removes it from both representations.
    /// tree.set_edge(root, child, None).unwrap();
    /// assert_eq!(tree.get(root).unwrap().get_child_edge(child), None);
    /// assert_eq!(tree.get(child).unwrap().parent_edge(), None);
    /// ```
    pub(crate) fn set_child_edge(&mut self, child: NodeId, edge: Option<EdgeLength>) {
        match edge {
            Some(e) => {
                self.child_edges
                    .get_or_insert_with(Default::default)
                    .insert(child, e);
            }
            None => {
                self.child_edges.as_mut().map(|edges| edges.remove(&child));
            }
        }
    }

    /// Removes the child from the node
    ///
    /// ```
    /// use phylotree::tree::{Node, Tree};
    /// let mut tree = Tree::new();
    /// let root = tree.add_root(Node::new()).unwrap();
    /// let child = tree.add_child(Node::new(), root, Some(0.5)).unwrap();
    /// tree.remove_child(root, child).unwrap();
    /// assert!(tree.get(root).unwrap().is_tip());
    /// // The removed child's storage slot is deleted.
    /// assert!(tree.get(child).is_err());
    /// ```
    pub(crate) fn remove_child(&mut self, child: NodeId) -> Result<(), NodeError> {
        let vec_index = match self.children.iter().position(|node_id| *node_id == child) {
            Some(idx) => idx,
            None => {
                return Err(NodeError::HasNoChild {
                    parent: self.id,
                    child,
                })
            }
        };

        self.children.remove(vec_index);
        self.child_edges.as_mut().map(|edges| edges.remove(&child));

        Ok(())
    }

    /// Scale this node's stored parent/child edge lengths and clear its distance cache.
    ///
    /// Apply to every node to keep both ends of every branch consistent. Mutable
    /// tree iterators also clear all tree caches before yielding any nodes.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("(A:1,B:2)R;")?;
    /// tree.iter_preorder_nodes_mut()?.try_for_each(|node| {
    ///     node?.rescale_edges(2.0);
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// assert_eq!(tree.to_newick()?, "(A:2,B:4)R;");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn rescale_edges(&mut self, factor: f64) {
        self.subtree_distances.get_mut().take();
        self.parent_edge = self.parent_edge.map(|edge| edge * factor);
        if let Some(edges) = &mut self.child_edges {
            for v in edges.values_mut() {
                *v *= factor;
            }
        }
    }

    /// Check if the node is a tip node
    pub fn is_tip(&self) -> bool {
        self.children.is_empty()
    }

    /// Check if the node is a root node
    pub fn is_root(&self) -> bool {
        self.parent.is_none()
    }

    /// Returns this node's fields in Newick format, excluding its children.
    ///
    /// ```
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// use phylotree::tree::NewickFormat;
    /// node.set_name("A");
    /// assert_eq!(node.to_newick(NewickFormat::AllFields).unwrap(), "A");
    /// ```
    pub fn to_newick(&self, format: NewickFormat) -> Result<String, super::TreeError> {
        super::newick::serialize_node(self, format)
    }
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        match (self.parent, other.parent) {
            (None, None) | (Some(_), Some(_)) => {}
            _ => return false,
        }

        if (self.parent.is_none() && other.parent.is_some())
            || (self.parent.is_some() && other.parent.is_none())
        {
            return false;
        }

        let parent_edges_equal = match (self.parent_edge, other.parent_edge) {
            (None, None) => true,
            (Some(l1), Some(l2)) => (l1 - l2).abs() < f64::EPSILON,
            _ => false,
        };

        let supports_equal = match (self.support, other.support) {
            (None, None) => true,
            (Some(s1), Some(s2)) => (s1 - s2).abs() < f64::EPSILON,
            _ => false,
        };

        self.name == other.name
            && self.children.len() == other.children.len()
            && parent_edges_equal
            && supports_equal
    }
}

impl Default for Node {
    fn default() -> Self {
        Self::new()
    }
}

impl Eq for Node {}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.parent_edge {
            Some(l) => write!(f, "({l:.3}) {:?}", self.name),
            None => write!(f, "{:?}", self.name),
        }
    }
}

impl Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({:?}) {:?} Id[{}] Parent[{:?}] Depth[{:?}] Support[{:?}] Comments[{:?}] Children({:?})",
            self.parent_edge,
            self.name,
            self.id,
            self.parent,
            self.depth,
            self.support,
            self.comment,
            self.children,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node_with_parent() -> Node {
        let mut node = Node::new();
        node.set_id(NodeId::from(7));
        node.set_parent(NodeId::from(0), Some(0.5));
        node
    }

    #[test]
    fn set_and_clear_name() {
        let mut node = Node::new();
        assert_eq!(node.name(), None);
        node.set_name("A");
        assert_eq!(node.name(), Some("A"));
        node.set_name(String::from("B"));
        assert_eq!(node.name(), Some("B"));
        node.clear_name();
        assert_eq!(node.name(), None);
    }

    #[test]
    fn set_and_clear_comment() {
        let mut node = Node::new();
        assert_eq!(node.comment(), None);
        node.set_comment("note");
        assert_eq!(node.comment(), Some("note"));
        node.set_comment(String::from("other"));
        assert_eq!(node.comment(), Some("other"));
        node.clear_comment();
        assert_eq!(node.comment(), None);
    }

    #[test]
    fn set_and_clear_support() {
        let mut node = node_with_parent();
        assert_eq!(node.support(), None);
        node.set_support(95.0).unwrap();
        assert_eq!(node.support(), Some(95.0));
        node.clear_support();
        assert_eq!(node.support(), None);
        node.set_support(0.5).unwrap();
        node.clear_support();
        assert_eq!(node.support(), None);
    }

    #[test]
    fn set_support_rejects_non_finite_values_and_leaves_node_unchanged() {
        let mut node = node_with_parent();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                node.set_support(bad),
                Err(NodeError::NonFiniteSupport { .. })
            ));
            assert_eq!(node.support(), None);
        }
    }

    #[test]
    fn set_support_allows_root_and_clear_support_is_infallible() {
        // Support is node metadata (e.g. root placement), not tied to a parent branch.
        let mut node = Node::new();
        node.set_support(1.0).unwrap();
        assert_eq!(node.support(), Some(1.0));
        node.clear_support();
        assert_eq!(node.support(), None);
    }
}
