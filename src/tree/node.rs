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
}

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
    /// Length of the branch between this node and its parent.
    pub parent_edge: Option<EdgeLength>,
    /// Support value associated with the branch leading to this node.
    ///
    /// The value is stored without normalization because both fractions and
    /// percentages are common in phylogenetic formats.
    pub support: Option<f64>,
    /// Optional comment attached to node
    pub comment: Option<String>,
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
    /// assert_eq!(node.name.as_deref(), Some("A"));
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

    /// Sets the internal Node name
    ///
    /// ```
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// node.set_name("A".into());
    /// assert_eq!(node.name.as_deref(), Some("A"));
    /// ```
    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }

    /// Sets the internal Node id
    pub fn set_id(&mut self, id: NodeId) {
        self.id = id;
    }

    /// Set this node's parent and incoming length.
    ///
    /// Update the parent's children separately; prefer [`super::Tree::add_child`]
    /// when editing a tree.
    ///
    /// ```
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// node.set_parent(NodeId::from(0), Some(0.5));
    /// assert_eq!(node.parent_edge, Some(0.5));
    /// ```
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
    /// use phylotree::tree::{Node, NodeId};
    ///
    /// let mut parent = Node::new();
    /// parent.id = NodeId::from(0);
    /// let mut child = Node::new();
    /// child.id = NodeId::from(1);
    ///
    /// let l = 0.1;
    ///
    /// child.set_parent(parent.id, Some(l));
    /// parent.add_child(child.id, Some(l));
    ///
    /// assert_eq!(child.parent_edge, Some(l));
    /// assert_eq!(parent.get_child_edge(child.id), Some(l));
    /// ```
    pub fn add_child(&mut self, child: NodeId, child_edge: Option<EdgeLength>) {
        self.children.push(child);
        self.set_child_edge(child, child_edge);
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
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// let child = NodeId::from(1);
    /// node.add_child(child, Some(0.5));
    /// node.set_child_edge(child, Some(1.0));
    /// assert_eq!(node.get_child_edge(child), Some(1.0));
    /// ```
    pub fn set_child_edge(&mut self, child: NodeId, edge: Option<EdgeLength>) {
        if let Some(edge) = edge {
            self.child_edges
                .get_or_insert_with(Default::default)
                .insert(child, edge);
        }
    }

    /// Removes the child from the node
    ///
    /// ```
    /// use phylotree::tree::{Node, NodeId};
    /// let mut node = Node::new();
    /// let child = NodeId::from(1);
    /// node.add_child(child, Some(0.5));
    /// node.remove_child(child).unwrap();
    /// assert!(node.is_tip());
    /// ```
    pub fn remove_child(&mut self, child: NodeId) -> Result<(), NodeError> {
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
    /// node.set_name("A".into());
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
