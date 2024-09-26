use std::{
    cell::RefCell,
    collections::HashMap,
    fmt::{Debug, Display},
};

use thiserror::Error;

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
}

use crate::tree::tree_impl::IdentityHasher;
type BuildIdentityHasher = core::hash::BuildHasherDefault<IdentityHasher>;

#[derive(Clone)]
/// A node of the Tree
pub struct Node {
    /// Index of the node
    pub id: NodeId,
    /// Name of the node
    pub name: Option<String>,
    /// Index of the parent node
    pub parent: Option<NodeId>,
    /// Indices of child nodes
    pub children: Vec<NodeId>,
    /// length of branch between parent and node
    pub parent_edge: Option<EdgeLength>,
    /// Optional comment attached to node
    pub comment: Option<String>,
    /// lenght of branches between node and children
    pub(crate) child_edges: Option<HashMap<NodeId, EdgeLength>>,
    /// Distance to descendants of this node
    pub(crate) subtree_distances: RefCell<Option<HashMap<NodeId, EdgeLength, BuildIdentityHasher>>>,
    /// Number of edges to root
    pub(crate) depth: usize,
    // Whether the node is deleted or not
    pub(crate) deleted: bool,
}

impl Node {
    /// Creates a new Node
    pub fn new() -> Self {
        Self {
            id: 0,
            name: None,
            parent: None,
            children: vec![],
            parent_edge: None,
            child_edges: None,
            subtree_distances: RefCell::new(None),
            comment: None,
            depth: 0,
            deleted: false,
        }
    }

    /// Creates a new named Node
    pub fn new_named(name: &str) -> Self {
        Self {
            id: 0,
            name: Some(String::from(name)),
            parent: None,
            children: vec![],
            parent_edge: None,
            child_edges: None,
            subtree_distances: RefCell::new(None),
            comment: None,
            depth: 0,
            deleted: false,
        }
    }

    /// Sets the internal Node name
    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }

    /// Sets the internal Node id
    pub fn set_id(&mut self, id: NodeId) {
        self.id = id;
    }

    /// Set the parent node
    /// See `add_child` for example usage
    pub fn set_parent(&mut self, parent: NodeId, parent_edge: Option<EdgeLength>) {
        self.parent = Some(parent);
        self.parent_edge = parent_edge;
    }

    /// Sets the depth of the node
    pub fn set_depth(&mut self, depth: usize) {
        self.depth = depth;
    }

    /// Gets the depth of the node
    pub fn get_depth(&self) -> usize {
        self.depth
    }

    /// Empties the node and sets it as deleted
    pub(crate) fn delete(&mut self) {
        *self = Self::new();
        self.deleted = true;
    }

    /// Adds a child to the node
    /// ```
    /// use phylotree::tree::Node;
    ///
    /// let mut parent = Node::new();
    /// parent.id = 0;
    /// let mut child = Node::new();
    /// child.id = 1;
    ///
    /// let l = 0.1;
    ///
    /// child.set_parent(parent.id, Some(l));
    /// parent.add_child(child.id, Some(l));
    ///
    /// assert_eq!(child.parent_edge, Some(l));
    /// assert_eq!(parent.get_child_edge(&child.id), Some(l));
    /// ```
    pub fn add_child(&mut self, child: NodeId, child_edge: Option<EdgeLength>) {
        self.children.push(child);
        self.set_child_edge(&child, child_edge);
    }

    /// Get the Edge between node and a child
    pub fn get_child_edge(&self, child: &NodeId) -> Option<EdgeLength> {
        if let Some(edges) = &self.child_edges {
            edges.get(child).copied()
        } else {
            None
        }
    }

    /// Sets the Edge between a node and its child
    pub fn set_child_edge(&mut self, child: &NodeId, edge: Option<EdgeLength>) {
        if let Some(edge) = edge {
            if self.child_edges.is_none() {
                self.child_edges = Some(HashMap::new());
            }
            self.child_edges
                .as_mut()
                .map(|edges| edges.insert(*child, edge));
        }
    }

    /// Removes the child from the node
    pub fn remove_child(&mut self, child: &NodeId) -> Result<(), NodeError> {
        let vec_index = match self.children.iter().position(|node_id| node_id == child) {
            Some(idx) => idx,
            None => {
                return Err(NodeError::HasNoChild {
                    parent: self.id,
                    child: *child,
                })
            }
        };

        self.children.remove(vec_index);
        self.child_edges.as_mut().map(|edges| edges.remove(child));

        Ok(())
    }

    /// Rescales parent and child edges by a factor.
    pub(crate) fn rescale_edges(&mut self, factor: f64) {
        self.parent_edge = self.parent_edge.map(|edge| edge * factor);
        if let Some(edges) = &mut self.child_edges {
            for (_, v) in edges.iter_mut() {
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

    fn format_name(&self) -> String {
        self.name.clone().unwrap_or_default()
    }

    fn format_length(&self) -> String {
        self.parent_edge
            .map(|v| format!(":{v}"))
            .unwrap_or_default()
    }

    fn format_comment(&self) -> String {
        self.comment
            .clone()
            .map(|v| format!("[{v}]"))
            .unwrap_or_default()
    }

    /// Returns String with node in newick format
    pub fn to_newick(&self, format: NewickFormat) -> String {
        let mut repr = String::new();

        match format {
            NewickFormat::AllFields
            | NewickFormat::NoComments
            | NewickFormat::OnlyNames
            | NewickFormat::LeafLengthsAllNames => repr += &self.format_name(),
            NewickFormat::LeafLengthsLeafNames
            | NewickFormat::InternalLengthsLeafNames
            | NewickFormat::AllLengthsLeafNames => {
                if self.is_tip() {
                    repr += &self.format_name()
                }
            }
            _ => (),
        }

        match format {
            NewickFormat::AllFields
            | NewickFormat::NoComments
            | NewickFormat::OnlyLengths
            | NewickFormat::AllLengthsLeafNames => repr += &self.format_length(),
            NewickFormat::InternalLengthsLeafNames => {
                if !self.is_tip() {
                    repr += &self.format_length()
                }
            }
            NewickFormat::LeafLengthsLeafNames | NewickFormat::LeafLengthsAllNames => {
                if self.is_tip() {
                    repr += &self.format_length()
                }
            }
            _ => (),
        }

        if let NewickFormat::AllFields = format {
            repr += &self.format_comment()
        }

        repr
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

        self.name == other.name && self.children.len() == other.children.len() && parent_edges_equal
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
            "({:?}) {:?} Id[{}] Parent[{:?}] Depth[{:?}] Comments[{:?}] Children({:?})",
            self.parent_edge,
            self.name,
            self.id,
            self.parent,
            self.depth,
            self.comment,
            self.children,
        )
    }
}
