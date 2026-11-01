use accurate::sum::NaiveSum;
use accurate::traits::*;

use itertools::Itertools;
#[cfg(not(target_arch = "wasm32"))]
use ptree::{print_tree, TreeBuilder};
use rand::seq::SliceRandom;
use std::collections::VecDeque;
use std::collections::{HashMap, HashSet};
use std::iter::zip;
#[cfg(not(target_arch = "wasm32"))]
use std::{fs, path::Path};
use vec_map::VecMap;

use thiserror::Error;

use super::bipartition::{
    BipartitionProfile, Partition, PartitionMap, PartitionSet, RobinsonFouldsMode,
};
use super::newick::{NewickParser, NewickSerializer, NewickTokenizer};
use super::node::{Node, NodeError};
use super::taxon_index::TaxonIndex;
use super::tree_list::Bipartition;
use super::{EdgeDepth, EdgeLength, NewickFormat, NewickParseError, NewickParseOptions, NodeId};

use crate::distance::{tril_to_rowvec_index, DistanceMatrix, MatrixError};

/// Errors that can occur when reading, writing and manipulating [`Tree`] structs.
#[derive(Error, Debug)]
pub enum TreeError {
    /// The tree is not binary and we are trying to do something
    /// only possible on binary trees
    #[error("This tree is not Binary.")]
    IsNotBinary,
    /// The tree is not rooted and we are trying to do something
    /// only possible on rooted trees
    #[error("This tree is not rooted.")]
    IsNotRooted,
    /// The tree is empty and we are trying to do something that require at least one node
    #[error("This tree is empty.")]
    IsEmpty,
    /// No root node was found in the tree and we are trying to do something
    /// that requires a root node
    #[error("No root node found")]
    RootNotFound,
    /// Some of the leaves in the tree have no name
    #[error("All your leaf nodes must ne named.")]
    UnnamedLeaves,
    /// Some of the leaves in the tree share the same name
    #[error("Your leaf names must be unique.")]
    DuplicateLeafNames,
    /// The leaf index is not initialized *(the leaf index is used when comparing tree topologies)*
    #[error("The leaf index of the tree is not initialized.")]
    LeafIndexNotInitialized,
    /// Some branches of the tree have no length
    #[error("The tree must have all branch lengths.")]
    MissingBranchLengths,
    /// The trees we want to compare have different tips
    #[error("The trees have different tips indices.")]
    DifferentTipIndices,
    /// Rooted comparison was requested for at least one unrooted tree.
    #[error("Rooted comparison requires two rooted trees.")]
    RootedComparisonRequiresRootedTrees,
    /// The requested node with index [`NodeId`] does not exist in the tree
    #[error("There is no node with index: {0}")]
    NodeNotFound(NodeId),
    /// The node with index [`NodeId`] could not be compressed
    #[error("Could not compress node {0}, it does not have exactly one parent and one child")]
    CouldNotCompressNode(NodeId),
    /// The two nodes could not be merged into a single parent
    #[error("Cound not merge nodes {0} and {1} since they are not siblings")]
    MergingNonSiblingNodes(NodeId, NodeId),
    /// There was a [`std::io::Error`] when writin the tree to a file
    #[error("Error writing tree to file")]
    IoError(#[from] std::io::Error),
    /// There was a [`NodeError`] when operating on a node
    #[error("Could operate on Node")]
    NodeError(#[from] NodeError),
    /// There was a [`MatrixError`] when extracting distance matrix
    #[error("Could not convert to matrix")]
    MatrixError(#[from] MatrixError),
    /// An internal node name and support value were both selected for the Newick label position.
    #[error("Node {0} cannot serialize both its name and support as an internal-node label")]
    ConflictingNewickNodeLabel(NodeId),
    /// A preserved support annotation conflicts with the node's support value.
    #[error("Node {0} has conflicting Newick support values")]
    ConflictingNewickSupport(NodeId),
    /// A preserved support annotation selected for serialization is not numeric.
    #[error("Node {0} has an invalid Newick support annotation")]
    InvalidNewickSupportComment(NodeId),
    /// A custom support key contains characters that cannot be serialized safely.
    #[error("Invalid Newick support comment key '{0}'")]
    InvalidNewickSupportKey(String),
    /// Tried to generate a tree with too little tips
    #[error("Tree must have at least 3 tips, got {0}")]
    InvalidTreeSize(usize),
    /// General error
    #[error("Encountered an error: {0}")]
    GeneralError(&'static str),
}

/// Struct to hold tree comparison metrics
#[derive(Debug, Clone)]
pub struct Comparison {
    /// Robinson Foulds metric
    pub rf: f64,
    /// Normalized Robinson Foulds
    pub norm_rf: f64,
    /// Weighted Robinson Foulds
    pub weighted_rf: f64,
    /// Kuhner Felsenstein Branch Score
    pub branch_score: f64,
}

/// Used to hold compared tree edges
type EdgeCompare = (
    Vec<(EdgeDepth, EdgeLength)>,
    Vec<(EdgeDepth, EdgeLength)>,
    Vec<((EdgeDepth, EdgeLength), (EdgeDepth, EdgeLength))>,
);

/// A Phylogenetic tree
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    taxon_index: Option<TaxonIndex>,
}

/// Base methods to add and get [`Node`] objects to and from the [`Tree`].
///
/// ----
/// ----
impl Tree {
    /// Create a new empty Tree object
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::new();
    /// assert_eq!(tree.size(), 0);
    /// ```
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            taxon_index: None,
        }
    }

    pub(crate) fn from_nodes(nodes: Vec<Node>) -> Self {
        let mut tree = Self {
            nodes,
            taxon_index: None,
        };
        tree.taxon_index = tree.compute_taxon_index().ok();
        tree
    }

    // ############################
    // # ADDING AND GETTING NODES #
    // ############################

    /// Add a new node to the tree.
    ///
    /// ```
    /// use phylotree::tree::{Tree, Node};
    /// let mut tree = Tree::new();
    /// let root = tree.add(Node::new_named("root"));
    /// assert_eq!(tree.root().unwrap(), root);
    /// ```
    pub fn add(&mut self, node: Node) -> NodeId {
        self.reset_bipartition_cache();
        let idx = NodeId(self.nodes.len());
        let mut node = node;
        node.id = idx;
        self.nodes.push(node);

        idx
    }

    /// Add a child to one of the tree's nodes.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::{Tree,Node};
    ///
    /// // Create the tree and add a root node
    /// let mut tree = Tree::new();
    /// let root_id = tree.add(Node::new());
    ///
    /// // Add children to the root
    /// let left = tree.add_child(Node::new(), root_id, None).unwrap();
    /// let right = tree.add_child(Node::new(), root_id, Some(0.1)).unwrap();
    ///
    /// assert_eq!(tree.get(root_id).unwrap().children.len(), 2);
    ///
    /// // The depths of child nodes are derived from the parent node
    /// assert_eq!(tree.get(left).unwrap().get_depth(), 1);
    /// assert_eq!(tree.get(right).unwrap().get_depth(), 1);
    ///
    /// // If an edge length is specified then it is set in both child and parent
    /// assert_eq!(tree.get(right).unwrap().parent_edge, Some(0.1));
    /// assert_eq!(tree.get(root_id).unwrap().get_child_edge(right), Some(0.1));
    /// ```
    pub fn add_child(
        &mut self,
        node: Node,
        parent: NodeId,
        edge: Option<EdgeLength>,
    ) -> Result<NodeId, TreeError> {
        if parent.index() >= self.nodes.len() {
            return Err(TreeError::NodeNotFound(parent));
        }

        let mut node = node;

        node.set_parent(parent, edge);
        node.set_depth(self.get(parent)?.depth + 1);

        let id = self.add(node);

        self.get_mut(id)?.set_id(id);
        self.get_mut(parent)?.add_child(id, edge);

        Ok(id)
    }

    /// Get a reference to a specific Node of the tree
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// let id = tree.root().unwrap();
    /// assert!(tree.get(id).unwrap().is_root());
    /// ```
    pub fn get(&self, id: NodeId) -> Result<&Node, TreeError> {
        if id.index() >= self.nodes.len() {
            return Err(TreeError::NodeNotFound(id));
        }
        let node = &self.nodes[id.index()];
        if node.deleted {
            return Err(TreeError::NodeNotFound(id));
        }

        Ok(node)
    }

    /// Get a mutable reference to a specific Node of the tree
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// let id = tree.root().unwrap();
    /// tree.get_mut(id).unwrap().set_name("root".into());
    /// assert_eq!(tree.get(id).unwrap().name.as_deref(), Some("root"));
    /// ```
    pub fn get_mut(&mut self, id: NodeId) -> Result<&mut Node, TreeError> {
        if id.index() >= self.nodes.len() || self.nodes[id.index()].deleted {
            return Err(TreeError::NodeNotFound(id));
        }
        self.reset_bipartition_cache();
        Ok(&mut self.nodes[id.index()])
    }

    /// Get a reference to a node in the tree by name.
    /// Note that this does not check for name unicity, if several nodes
    /// match a name this funciton will return the first match in the tree.
    /// If you want to find all nodes matching a name in a given tree,
    /// use [`Tree::search_nodes`].
    /// ```
    /// use phylotree::tree::{Tree, Node};
    ///
    /// let mut tree = Tree::new();
    /// let root_idx = tree.add(Node::new_named("root"));
    /// let child_idx = tree.add_child(Node::new_named("child"), root_idx, None).unwrap();
    ///
    /// assert_eq!(tree.get_by_name("child"), Some(tree.get(child_idx).unwrap()));
    /// ```
    pub fn get_by_name(&self, name: &str) -> Option<&Node> {
        self.nodes
            .iter()
            .find(|node| node.name.as_deref().is_some_and(|n| n == name))
    }

    /// Get a mutable reference to a node in the tree by name
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// tree.get_by_name_mut("A").unwrap().set_name("X".into());
    /// assert!(tree.get_by_name("X").is_some());
    /// ```
    pub fn get_by_name_mut(&mut self, name: &str) -> Option<&mut Node> {
        let index = self
            .nodes
            .iter()
            .position(|node| node.name.as_deref().is_some_and(|n| n == name))?;
        self.reset_bipartition_cache();
        Some(&mut self.nodes[index])
    }

    /// Search nodes in the tree with a closure.
    /// ```
    /// use phylotree::tree::{Tree, Node};
    ///
    /// let mut tree = Tree::new();
    /// let root_idx = tree.add(Node::new_named("root"));
    /// let mut indices = vec![];
    ///
    /// for name in ["A", "B", "A"] {
    ///     let idx = tree.add_child(Node::new_named(name), root_idx, None).unwrap();
    ///     if name == "A" { indices.push(idx) }
    /// }
    ///
    /// let found = tree.search_nodes(|node| node.name == Some("A".into()));
    /// assert_eq!(found, indices);
    /// ```
    pub fn search_nodes(&self, cond: impl Fn(&Node) -> bool) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|node| cond(node))
            .map(|node| node.id)
            .collect()
    }

    /// Gets the root node. In the case of unrooted trees this node is a "virtual root"
    /// that has exactly 3 children.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.get(tree.root().unwrap()).unwrap().is_root());
    /// ```
    pub fn root(&self) -> Result<NodeId, TreeError> {
        self.nodes
            .iter()
            .filter(|&node| node.parent.is_none())
            .map(|node| node.id)
            .next()
            .ok_or(TreeError::RootNotFound)
    }

    /// Returns a [`Vec`] containing the Node IDs of leaf nodes of the tree
    /// ```
    /// use phylotree::tree::{Tree, Node};
    ///
    /// let mut tree = Tree::new();
    /// let root_idx = tree.add(Node::new());
    /// let left = tree.add_child(Node::new(), root_idx, None).unwrap();
    /// let right = tree.add_child(Node::new(), root_idx, None).unwrap();
    ///
    /// assert_eq!(tree.leaves(), vec![left, right]);
    /// ```
    pub fn leaves(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|&node| !node.deleted && node.is_tip())
            .map(|node| node.id)
            .collect()
    }

    /// Returns a [`Vec`] containing the Names of the leaf nodes of the tree
    /// ```
    /// use phylotree::tree::{Tree, Node};
    ///
    /// let mut tree = Tree::new();
    /// let root_idx = tree.add(Node::new());
    /// let _ = tree.add_child(Node::new_named("left"), root_idx, None).unwrap();
    /// let _ = tree.add_child(Node::new_named("right"), root_idx, None).unwrap();
    ///
    /// let names: Vec<_> = tree.leaf_names()
    ///     .into_iter()
    ///     .flatten()
    ///     .collect();
    /// assert_eq!(names, vec!["left", "right"]);
    /// ```
    pub fn leaf_names(&self) -> Vec<Option<String>> {
        self.leaves()
            .iter()
            .map(|leaf_id| self.get(*leaf_id).unwrap().name.clone())
            .collect()
    }

    /// Gets the node ids of all the nodes in the subtree rooted at the specified node
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let subtree: Vec<_> = tree.subtree(sub_root.id)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(subtree, vec!["E", "C", "D"])
    /// ```
    pub fn subtree(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![root];

        for child in self.get(root)?.children.iter() {
            indices.extend(self.subtree(*child)?);
        }

        Ok(indices)
    }

    /// Gets the node ids of all the nodes in the subtree rooted at the specified node, except the root
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let subtree: Vec<_> = tree.descendants(sub_root.id)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(subtree, vec!["C", "D"])
    /// ```
    pub fn descendants(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];

        for child in self.get(root)?.children.iter() {
            indices.extend(self.subtree(*child)?);
        }

        Ok(indices)
    }

    /// Gets the node ids of all the leaves in the subtree rooted at the specified node
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let sub_leaves: Vec<_> = tree.subtree_leaves(sub_root.id)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(sub_leaves, vec!["C", "D"])
    /// ```
    pub fn subtree_leaves(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        Ok(self
            .subtree(root)?
            .into_iter()
            .filter(|id| self.get(*id).unwrap().is_tip())
            .collect())
    }
}

/// Methods to traverse the [`Tree`]
///
/// ----
/// ----
impl Tree {
    // ###################
    // # TREE TRAVERSALS #
    // ###################

    /// Returns a vector containing node ids in the same order as the
    /// [preorder](https://en.wikipedia.org/wiki/Tree_traversal#Pre-order,_NLR) tree traversal
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let preorder: Vec<_> = tree.preorder(tree.root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(preorder, vec!["F", "B", "A", "D", "C", "E", "G", "I", "H"])
    /// ```
    pub fn preorder(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![root];
        for child in self.get(root)?.children.iter() {
            indices.extend(self.preorder(*child)?)
        }

        Ok(indices)
    }

    /// Returns a vector containing node ids in the same order as the
    /// [postorder](https://en.wikipedia.org/wiki/Tree_traversal#Post-order,_LRN ) tree traversal
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let postorder: Vec<_> = tree.postorder(tree.root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(postorder, vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"])
    /// ```
    pub fn postorder(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];
        for child in self.get(root)?.children.iter() {
            indices.extend(self.postorder(*child)?)
        }
        indices.push(root);

        Ok(indices)
    }

    /// Returns a vector containing node ids in the same order as the
    /// [inorder](https://en.wikipedia.org/wiki/Tree_traversal#In-order,_LNR) tree traversal.
    /// This assumes that the tree is binary.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let inorder: Vec<_> = tree.inorder(tree.root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(inorder, vec!["A", "B", "C", "D", "E", "F", "H", "I", "G"])
    /// ```
    /// *N.B.: the resulting traversal in the example is different than the one in the
    /// [wikipedia](https://en.wikipedia.org/wiki/Tree_traversal) page. That is because
    /// the [`Node`] struct cannot have a right child only, so in our tree I is the left
    /// child of G.*
    pub fn inorder(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];
        let children = &self.get(root)?.children;

        // Tree is not binary
        if children.len() > 2 {
            return Err(TreeError::IsNotBinary);
        }

        // There is a left child
        if !children.is_empty() {
            indices.extend(self.inorder(children[0])?)
        }

        indices.push(root);

        // There is a right child
        if children.len() > 1 {
            indices.extend(self.inorder(children[1])?)
        }

        Ok(indices)
    }

    /// Returns a vector containing node ids in the same order as the
    /// [levelorder](https://en.wikipedia.org/wiki/Tree_traversal#Breadth-first_search) tree traversal
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let levelorder: Vec<_> = tree.levelorder(tree.root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(levelorder, vec!["F", "B", "G", "A", "D", "I", "C", "E", "H"])
    /// ```
    pub fn levelorder(&self, root: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];
        let mut queue = VecDeque::new();
        queue.push_back(root);
        while !queue.is_empty() {
            let root = queue.pop_front().unwrap();
            indices.push(root);
            for child in self.get(root)?.children.iter() {
                queue.push_back(*child)
            }
        }

        Ok(indices)
    }
}

/// Methods that compute characteristics and measures to describe the [`Tree`]
///
/// ----
/// ----
impl Tree {
    // #######################################
    // # GETTING CHARACTERISTICS OF THE TREE #
    // #######################################

    /// Check if the tree is Binary
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.is_binary().unwrap());
    /// ```
    pub fn is_binary(&self) -> Result<bool, TreeError> {
        for node in self.nodes.iter() {
            // Root of the tree
            if node.parent.is_none() {
                if self.is_rooted()? && node.children.len() > 2 {
                    return Ok(false);
                // The virtual root of unrooted trees can has up to 3 children
                } else if !self.is_rooted()? && node.children.len() > 3 {
                    return Ok(false);
                }
            } else if node.parent.is_some() && node.children.len() > 2 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Checks if the tree is rooted (i.e. the root node exists and has exactly 2 children)
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.is_rooted().unwrap());
    /// ```
    pub fn is_rooted(&self) -> Result<bool, TreeError> {
        let root_id = self.root()?;

        Ok(!self.nodes.is_empty() && self.get(root_id)?.children.len() == 2)
    }

    /// Checks if all the tips have unique names (This check assumes that all tips have a name)
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.has_unique_tip_names().unwrap());
    /// ```
    pub fn has_unique_tip_names(&self) -> Result<bool, TreeError> {
        let mut names = HashSet::new();
        for name in self.leaf_names() {
            if let Some(name) = name {
                names.insert(name);
            } else {
                return Err(TreeError::UnnamedLeaves);
            }
        }

        Ok(names.len() == self.n_leaves())
    }

    /// Returns the number of nodes in the tree
    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the number of leaves in the tree
    pub fn n_leaves(&self) -> usize {
        self.nodes.iter().filter(|&node| node.is_tip()).count()
    }

    /// Returns the height of the tree
    /// (i.e. the number of edges or branch length sum from the root to the deepest tip)
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A:0.1,B:0.2)G:0.1,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// assert_eq!(tree.height().unwrap(), 0.9);
    ///
    /// let tree_no_brlen = Tree::from_newick("((A,B)G,(C,D)E)F;").unwrap();
    /// assert_eq!(tree_no_brlen.height().unwrap(), 2.);
    /// ```
    pub fn height(&self) -> Result<EdgeLength, TreeError> {
        if !self.is_rooted()? {
            return Err(TreeError::IsNotRooted);
        }

        let root = self.root()?;

        self.leaves()
            .iter()
            .map(|&leaf| {
                let (edge_sum, num_edges) = self.get_distance(root, leaf).unwrap();
                match edge_sum {
                    Some(height) => height,
                    None => num_edges as f64,
                }
            })
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .ok_or(TreeError::IsEmpty)
    }

    /// Returns the diameter of the tree
    /// (i.e. longest tip to tip distance), in O(nodes) time and space.
    /// Paths with any missing branch length use their edge count instead.
    /// The root stem is ignored. Fewer than two tips returns `IsEmpty`.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// assert_eq!(tree.diameter().unwrap(), 1.1);
    ///
    /// let tree_no_brlen = Tree::from_newick("(A,B,(C,D)E)F;").unwrap();
    /// assert_eq!(tree_no_brlen.diameter().unwrap(), 3.);
    /// ```
    pub fn diameter(&self) -> Result<EdgeLength, TreeError> {
        if self.leaves().len() < 2 {
            return Err(TreeError::IsEmpty);
        }
        // Iterative postorder avoids recursive traversal on deep trees.
        let mut pending = vec![self.root()?];
        let mut order = Vec::new();
        while let Some(id) = pending.pop() {
            order.push(id);
            pending.extend(self.get(id)?.children.iter().rev().copied());
        }
        let mut weighted = vec![None::<f64>; self.nodes.len()];
        let mut hops = vec![0_usize; self.nodes.len()];
        let mut missing = vec![None::<usize>; self.nodes.len()];
        let mut diameter = None::<f64>;
        for id in order.into_iter().rev() {
            let node = self.get(id)?;
            if node.is_tip() {
                weighted[id.index()] = Some(0.0);
                continue;
            }
            let mut best_weighted = None::<f64>;
            let mut best_hops = None::<usize>;
            let mut best_missing = None::<usize>;
            for &child in &node.children {
                let edge = self.get(child)?.parent_edge;
                if edge.is_some_and(|length| !length.is_finite()) {
                    return Err(TreeError::GeneralError("Nonfinite branch length"));
                }
                let distance = edge.zip(weighted[child.index()]).map(|(a, b)| a + b);
                let depth = hops[child.index()] + 1;
                // Preserve the pairwise API: any missing edge makes the entire
                // tip-to-tip path topological, rather than mixing units.
                let missing_depth = if edge.is_none() {
                    Some(depth)
                } else {
                    missing[child.index()].map(|d| d + 1)
                };
                let candidates = [
                    best_weighted.zip(distance).map(|(a, b)| a + b),
                    best_missing.map(|d| (d + depth) as f64),
                    best_hops.zip(missing_depth).map(|(a, b)| (a + b) as f64),
                ];
                for candidate in candidates.into_iter().flatten() {
                    if !candidate.is_finite() {
                        return Err(TreeError::GeneralError("Nonfinite diameter"));
                    }
                    diameter = Some(diameter.map_or(candidate, |d| d.max(candidate)));
                }
                // Only actual leaf endpoints count, including negative-length
                // paths; an internal node must not introduce a zero endpoint.
                if let Some(d) = distance {
                    best_weighted = Some(best_weighted.map_or(d, |old| old.max(d)));
                }
                best_hops = Some(best_hops.map_or(depth, |old| old.max(depth)));
                if let Some(d) = missing_depth {
                    best_missing = Some(best_missing.map_or(d, |old| old.max(d)));
                }
            }
            weighted[id.index()] = best_weighted;
            hops[id.index()] = best_hops.unwrap_or(0);
            missing[id.index()] = best_missing;
        }
        diameter.ok_or(TreeError::IsEmpty)
    }

    /// Returns the lenght of the trees
    /// (i.e. the sum of branch lenghts)
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// assert_eq!(tree.length().unwrap(), 1.5);
    /// ```
    pub fn length(&self) -> Result<EdgeLength, TreeError> {
        let s = self
            .nodes
            .iter()
            .filter(|n| !n.is_root())
            .map(|n| n.parent_edge)
            .collect::<Option<Vec<_>>>();
        match s {
            Some(v) => Ok(v.iter().sum()),
            None => Err(TreeError::MissingBranchLengths),
        }
    }

    /// Checks if the tree is rooted and binary
    fn check_rooted_binary(&self) -> Result<(), TreeError> {
        if !self.is_rooted()? {
            Err(TreeError::IsNotRooted)
        } else if !self.is_binary()? {
            Err(TreeError::IsNotBinary)
        } else {
            Ok(())
        }
    }

    /// Computes the number of cherries in a tree
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.n_cherries().unwrap(), 2);
    /// ```
    pub fn n_cherries(&self) -> Result<usize, TreeError> {
        if !self.is_binary()? {
            return Err(TreeError::IsNotBinary);
        }
        if !self.nodes.is_empty() {
            let mut n = 0;
            for node in self.nodes.iter() {
                if node.children.len() == 2
                    && self.get(node.children[0])?.is_tip()
                    && self.get(node.children[1])?.is_tip()
                {
                    n += 1;
                }
            }
            Ok(n)
        } else {
            Err(TreeError::IsEmpty)
        }
    }

    /// Sum absolute differences in descendant tip counts at each bifurcation.
    /// Requires a rooted binary tree; zero indicates a balanced topology.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.colless().unwrap(), 0);
    /// ```
    pub fn colless(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;

        let mut colless = 0;

        for node in self.nodes.iter().filter(|node| !node.is_tip()) {
            let left = self.subtree_leaves(node.children[0])?.len();
            let right = if node.children.len() > 1 {
                self.subtree_leaves(node.children[1])?.len()
            } else {
                0
            };

            colless += left.abs_diff(right);
        }

        Ok(colless)
    }

    /// Return Colless normalized against the Yule expectation.
    ///
    /// Uses `(Ic - n * ln(n) - n * (gamma - 1 - ln(2))) / n`, where `gamma`
    /// is Euler's constant; see [`Self::colless`].
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.colless_yule().unwrap().is_finite());
    /// ```
    pub fn colless_yule(&self) -> Result<f64, TreeError> {
        self.colless().map(|i_c| {
            let n = self.n_leaves() as f64;
            let e_i_c = n * n.ln() + (0.57721566 - 1. - f64::ln(2.0)) * n;

            (i_c as f64 - e_i_c) / n
        })
    }

    /// Return Colless divided by `n^(3/2)` for a PDA null model.
    /// Requires a rooted binary tree; see [`Self::colless`].
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.colless_pda().unwrap(), 0.0);
    /// ```
    pub fn colless_pda(&self) -> Result<f64, TreeError> {
        self.colless()
            .map(|i_c| i_c as f64 / f64::powf(self.n_leaves() as f64, 3.0 / 2.0))
    }

    /// Sum the depths of all tips (the Sackin imbalance index).
    /// Requires a rooted binary tree; smaller values indicate greater balance.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.sackin().unwrap(), 8);
    /// ```
    pub fn sackin(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;

        Ok(self
            .leaves()
            .iter()
            .map(|&tip_idx| self.get(tip_idx).unwrap().depth)
            .sum())
    }

    /// Return Sackin normalized against the Yule expectation.
    /// Uses `(Is - 2n * sum(1/j, j=2..n)) / n`; see [`Self::sackin`].
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.sackin_yule().unwrap().is_finite());
    /// ```
    pub fn sackin_yule(&self) -> Result<f64, TreeError> {
        self.sackin().map(|i_n| {
            let n = self.n_leaves();
            let sum: f64 = (2..=n).map(|i| 1.0 / (i as f64)).sum();

            (i_n as f64 - 2.0 * (n as f64) * sum) / n as f64
        })
    }

    /// Return Sackin divided by `n^(3/2)` for a PDA null model.
    /// Requires a rooted binary tree; see [`Self::sackin`].
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.sackin_pda().unwrap(), 1.0);
    /// ```
    pub fn sackin_pda(&self) -> Result<f64, TreeError> {
        self.sackin()
            .map(|i_n| i_n as f64 / f64::powf(self.n_leaves() as f64, 3.0 / 2.0))
    }
}

/// Methods that compute edge bipartitions and compare [`Tree`] objects with each other.
///
/// ----
/// ----
impl Tree {
    // #########################
    // # GET EDGES IN THE TREE #
    // #########################

    fn compute_taxon_index(&self) -> Result<TaxonIndex, TreeError> {
        if self.nodes.is_empty() {
            return Err(TreeError::IsEmpty);
        }
        if !self.has_unique_tip_names()? {
            return Err(TreeError::DuplicateLeafNames);
        }

        let taxa = self.leaf_names().into_iter().flatten().sorted().collect();
        Ok(TaxonIndex::from_sorted_taxa(taxa))
    }

    /// Returns the stored taxon index, if it is currently available.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// tree.rebuild_taxon_index().unwrap();
    /// assert_eq!(tree.cached_taxon_index().unwrap().len(), 4);
    /// ```
    pub fn cached_taxon_index(&self) -> Option<&TaxonIndex> {
        self.taxon_index.as_ref()
    }

    /// Returns the mapping between taxon labels and partition bit positions.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((C,A),B);").unwrap();
    /// let index = tree.taxon_index().unwrap();
    ///
    /// assert_eq!(index.taxa(), &["A", "B", "C"]);
    /// assert_eq!(index.index_of("C"), Some(2));
    /// ```
    pub fn taxon_index(&self) -> Result<TaxonIndex, TreeError> {
        self.taxon_index
            .clone()
            .map(Ok)
            .unwrap_or_else(|| self.compute_taxon_index())
    }

    /// Recomputes and stores the mapping between taxon labels and bit positions.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.rebuild_taxon_index().unwrap().len(), 4);
    /// ```
    pub fn rebuild_taxon_index(&mut self) -> Result<&TaxonIndex, TreeError> {
        self.taxon_index = None;
        self.taxon_index = Some(self.compute_taxon_index()?);
        Ok(self.taxon_index.as_ref().unwrap())
    }

    pub(crate) fn bind_taxon_index(&mut self, taxon_index: TaxonIndex) -> Result<(), TreeError> {
        if !self.taxon_index()?.is_compatible(&taxon_index) {
            return Err(TreeError::DifferentTipIndices);
        }
        self.taxon_index = Some(taxon_index);
        Ok(())
    }

    /// Converts the stored side of a partition to concatenated taxon labels.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// let partition = tree.get_partitions().unwrap().into_iter().next().unwrap();
    /// assert_eq!(tree.partition_to_leaves(&partition).unwrap().chars().count(), 2);
    /// ```
    pub fn partition_to_leaves(&self, partition: &Partition) -> Result<String, TreeError> {
        let taxon_index = self.taxon_index()?;
        Ok(partition
            .ones()
            .map(|index| taxon_index.taxon(index).unwrap())
            .collect())
    }

    /// Computes an immutable bipartition profile for repeated comparisons.
    ///
    /// Creating a profile does not populate or otherwise modify this tree's
    /// optional taxon-index cache. The storage strategy controls whether the
    /// canonical unrooted split projection is retained by the profile.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// let profile = tree.bipartition_profile().unwrap();
    /// assert_eq!(profile.taxon_index().len(), 4);
    /// ```
    pub fn bipartition_profile(&self) -> Result<BipartitionProfile, TreeError> {
        BipartitionProfile::from_tree(self)
    }

    /// Get all partitions of a tree.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.get_partitions().unwrap().len(), 1);
    /// ```
    pub fn get_partitions(&self) -> Result<PartitionSet, TreeError> {
        self.bipartition_profile()?
            .partition_set(RobinsonFouldsMode::Unrooted)
    }

    /// Returns all nontrivial bipartitions with their taxon index.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,B),(C,D));").unwrap();
    /// let partitions = tree.bipartitions().unwrap();
    /// let partition = partitions.iter().next().unwrap();
    ///
    /// assert_eq!(partition.taxa().collect::<Vec<_>>(), vec!["A", "B"]);
    /// assert_eq!(partition.complement_taxa().collect::<Vec<_>>(), vec!["C", "D"]);
    /// ```
    pub fn bipartitions(&self) -> Result<HashSet<Bipartition>, TreeError> {
        let profile = self.bipartition_profile()?;
        let taxon_index = profile.taxon_index().clone();
        Ok(profile
            .partition_set(RobinsonFouldsMode::Unrooted)?
            .into_iter()
            .map(|bits| Bipartition::new(taxon_index.clone(), bits))
            .collect())
    }

    /// Get all partitions of a tree along with corresponding branch lengths and branch depths.
    pub(crate) fn get_partitions_with_lengths(&self) -> Result<PartitionMap, TreeError> {
        self.bipartition_profile()?.partitions_with_lengths()
    }

    /// Invalidates the cached taxon index.
    ///
    /// Bipartitions are no longer cached on trees; use [`Tree::bipartition_profile`]
    /// to retain a computed representation for repeated comparisons.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// tree.rebuild_taxon_index().unwrap();
    /// tree.reset_bipartition_cache();
    /// assert!(tree.cached_taxon_index().is_none());
    /// ```
    pub fn reset_bipartition_cache(&mut self) {
        self.taxon_index = None;
    }

    // #################
    // # COMPARE TREES #
    // #################

    /// Computes the Robinson–Foulds distance using an explicit rooted or
    /// unrooted interpretation.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// use phylotree::tree::RobinsonFouldsMode;
    /// assert_eq!(tree.robinson_foulds(&tree, RobinsonFouldsMode::Unrooted).unwrap(), 0);
    /// ```
    pub fn robinson_foulds(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<usize, TreeError> {
        self.bipartition_profile()?
            .robinson_foulds(&other.bipartition_profile()?, mode)
    }

    /// Computes rooted Robinson–Foulds distance from descendant clades.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.rooted_robinson_foulds(&tree).unwrap(), 0);
    /// ```
    pub fn rooted_robinson_foulds(&self, other: &Self) -> Result<usize, TreeError> {
        self.robinson_foulds(other, RobinsonFouldsMode::Rooted)
    }

    /// Computes unrooted Robinson–Foulds distance from canonical splits.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.unrooted_robinson_foulds(&tree).unwrap(), 0);
    /// ```
    pub fn unrooted_robinson_foulds(&self, other: &Self) -> Result<usize, TreeError> {
        self.robinson_foulds(other, RobinsonFouldsMode::Unrooted)
    }

    /// Return RF divided by the total number of groupings in both trees.
    /// The result lies in `[0, 1]`; two empty grouping sets have distance zero.
    /// See [`Self::robinson_foulds`] for the rooted/unrooted interpretation.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// use phylotree::tree::RobinsonFouldsMode;
    /// assert_eq!(tree.robinson_foulds_norm(&tree, RobinsonFouldsMode::Unrooted).unwrap(), 0.0);
    /// ```
    pub fn robinson_foulds_norm(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<f64, TreeError> {
        self.bipartition_profile()?
            .robinson_foulds_norm(&other.bipartition_profile()?, mode)
    }

    /// Sum absolute branch-length differences over unrooted splits.
    /// Unmatched splits contribute their full length; requires matching taxa
    /// and lengths on the relevant internal branches.
    /// [Robinson & Foulds (1979)](https://doi.org/10.1007/BFb0102690).
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.weighted_robinson_foulds(&tree).unwrap(), 0.0);
    /// ```
    pub fn weighted_robinson_foulds(&self, other: &Self) -> Result<f64, TreeError> {
        self.bipartition_profile()?
            .weighted_robinson_foulds(&other.bipartition_profile()?, RobinsonFouldsMode::Unrooted)
    }

    /// Return the Euclidean branch score over unrooted splits.
    /// Unmatched splits are compared with length zero; requires matching taxa
    /// and lengths on the relevant internal branches.
    /// [Kuhner & Felsenstein (1994)](https://doi.org/10.1093/oxfordjournals.molbev.a040126).
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert_eq!(tree.kuhner_felsenstein(&tree).unwrap(), 0.0);
    /// ```
    pub fn kuhner_felsenstein(&self, other: &Self) -> Result<f64, TreeError> {
        self.bipartition_profile()?
            .kuhner_felsenstein(&other.bipartition_profile()?, RobinsonFouldsMode::Unrooted)
    }

    /// Compute several the RF metric, the weighted and normalized RF metrics and
    /// the KF branch score in one pass. This is more efficient than calling the
    /// different functions separately.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree1 = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let tree2 = Tree::from_newick("(A:0.1,D:0.2,(C:0.3,B:0.4)E:0.5)F;").unwrap();
    ///
    /// let rf = tree1.unrooted_robinson_foulds(&tree2).unwrap() as f64;
    /// let norm_rf = tree1.robinson_foulds_norm(
    ///     &tree2,
    ///     phylotree::tree::RobinsonFouldsMode::Unrooted,
    /// ).unwrap();
    /// let weighted_rf = tree1.weighted_robinson_foulds(&tree2).unwrap();
    /// let branch_score = tree1.kuhner_felsenstein(&tree2).unwrap();
    ///
    /// let comparison = tree1.compare_topologies(&tree2).unwrap();
    ///
    /// assert_eq!(rf, comparison.rf);
    /// assert_eq!(norm_rf, comparison.norm_rf);
    /// assert_eq!(weighted_rf, comparison.weighted_rf);
    /// assert_eq!(branch_score, comparison.branch_score);
    /// ```
    pub fn compare_topologies(&self, other: &Self) -> Result<Comparison, TreeError> {
        self.bipartition_profile()?
            .compare(&other.bipartition_profile()?, RobinsonFouldsMode::Unrooted)
    }

    /// Compare sets of branches between 2 trees. This will return 3
    /// sets of branch lengths:
    ///  - branches exclusive to the `self` tree
    ///  - branches exclusive to the `other` tree
    ///  - 2-uples for common branches: (self_len, other_len)
    ///
    ///  You can specify wether to include terminal branches in this calculation or not
    ///  with the `include_terminal` parameter.
    ///
    ///  # Example
    ///  ```
    ///  use phylotree::tree::Tree;
    ///
    ///  // Ref Tree:
    ///  //               0.3
    ///  //         0.2 +----- A
    ///  //   0.1 +-----| 0.4
    ///  // +-----|     +----- B
    ///  // |     | 0.5
    ///  // |     +----- C
    ///  // |       0.7
    ///  // | 0.6 +----- D
    ///  // +-----| 0.8
    ///  //       +----- E
    ///
    ///  // Compared Tree:
    ///  //               0.3
    ///  //         0.2 +----- A
    ///  //   0.1 +-----| 0.4
    ///  // +-----|     +----- C
    ///  // |     | 0.5
    ///  // |     +----- B
    ///  // |       7.0
    ///  // | 0.6 +----- D
    ///  // +-----| 0.8
    ///  //       +----- E
    ///  //
    ///  // We just switched the B and C labels, so we should have 1
    ///  // exlusive branch per tree. The compared tree also has one branch
    ///  // that is 7.0 instead of 0.7.
    ///
    ///  let reftree = Tree::from_newick("(((A:0.3,B:0.4):0.2,C:0.5):0.1,(D:0.7,E:0.8):0.6);").unwrap();
    ///  let cmptree = Tree::from_newick("(((A:0.3,C:0.4):0.2,B:0.5):0.1,(D:7.0,E:0.8):0.6);").unwrap();
    ///
    ///  // Get branch comparison including terminal branches
    ///  let (refset, cmpset, common) = reftree.compare_branch_lengths(&cmptree, true).unwrap();
    ///
    ///  assert_eq!(refset, vec![(2,0.2)]); // AB,CDE exclusive to reftree
    ///  assert_eq!(cmpset, vec![(2,0.2)]); // AC,BDE exclusive to cmptree
    ///
    ///  let expected_common = [
    ///     (0.3, 0.3), // A,...
    ///     (0.4, 0.4), // B,...
    ///     (0.5, 0.5), // C,...
    ///     (0.7, 0.7), // ABC,DE
    ///     (0.7, 7.0), // D,...
    ///     (0.8, 0.8), // E,...
    ///  ];
    ///
    ///  assert_eq!(expected_common.len(), common.len());
    ///  for val in common.into_iter().map(|((_d1,l1),(_d2,l2))| (l1,l2)) {
    ///     let mut found = false;
    ///     for exp in expected_common.iter() {
    ///         if val.0 - exp.0 < f64::EPSILON && val.1 - exp.1 < f64::EPSILON {
    ///             found = true;
    ///         }
    ///     }
    ///     assert!(found);
    ///  }
    ///  ```
    pub fn compare_branch_lengths(
        &self,
        other: &Self,
        include_tips: bool,
    ) -> Result<EdgeCompare, TreeError> {
        let partitions_s = self.get_partitions_with_lengths()?;
        let partitions_o = other.get_partitions_with_lengths()?;

        let mut self_branches = vec![];
        let mut other_branches = vec![];
        let mut common_branches = vec![];

        for (edge, len_s) in partitions_s.iter() {
            if let Some(len_o) = partitions_o.get(edge) {
                common_branches.push((*len_s, *len_o));
            } else {
                self_branches.push(*len_s);
            }
        }

        for (edge, len_o) in partitions_o.iter() {
            if !partitions_s.contains_key(edge) {
                other_branches.push(*len_o);
            }
        }

        if include_tips {
            let selftips = self.get_terminal_branches()?;
            let othertips = other.get_terminal_branches()?;

            for (tipname, (depth_s, len_s)) in selftips.iter() {
                let len_s = len_s.ok_or(TreeError::MissingBranchLengths)?;
                if let Some((depth_o, len_o)) = othertips.get(tipname) {
                    let len_o = len_o.ok_or(TreeError::MissingBranchLengths)?;
                    common_branches.push(((*depth_s, len_s), (*depth_o, len_o)));
                } else {
                    self_branches.push((*depth_s, len_s));
                }
            }

            for (tipname, (depth_o, len_o)) in othertips.iter() {
                if !selftips.contains_key(tipname) {
                    let len_o = len_o.ok_or(TreeError::MissingBranchLengths)?;
                    other_branches.push((*depth_o, len_o));
                }
            }
        }

        Ok((self_branches, other_branches, common_branches))
    }

    // Get terminal branch lengths of a tree keyed by tip name
    fn get_terminal_branches(
        &self,
    ) -> Result<HashMap<String, (usize, Option<EdgeLength>)>, TreeError> {
        if !self.has_unique_tip_names()? {
            return Err(TreeError::DuplicateLeafNames);
        }

        Ok(HashMap::from_iter(self.leaves().iter().map(|&idx| {
            let node = self.get(idx).unwrap();
            (node.name.clone().unwrap(), (node.depth, node.parent_edge))
        })))
    }
}

/// Methods to find paths in a [`Tree`] as well as measure distances between [`Node`] objects.
///
/// ----
/// ----
impl Tree {
    // ##########################
    // # FIND PATHS IN THE TREE #
    // ##########################

    /// Returns the path from the node to the root
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let path: Vec<_> = tree.get_path_from_root(5.into())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|&id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(path, vec!["F", "B", "D", "E"])
    /// ```
    pub fn get_path_from_root(&self, node: NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut path = vec![];
        let mut current_node = node;
        loop {
            path.push(current_node);
            match self.get(current_node)?.parent {
                Some(parent) => current_node = parent,
                None => break,
            }
        }

        Ok(path.into_iter().rev().collect())
    }

    /// Gets the most recent common ancestor between two tree nodes
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let ancestor = tree.get_common_ancestor(
    ///     tree.get_by_name("A").unwrap().id,
    ///     tree.get_by_name("D").unwrap().id,
    /// ).unwrap();
    ///
    /// assert_eq!(tree.get(ancestor).unwrap().name, Some("B".to_owned()))
    /// ```
    pub fn get_common_ancestor(&self, source: NodeId, target: NodeId) -> Result<NodeId, TreeError> {
        if source == target {
            return Ok(source);
        }
        let root_to_source = self.get_path_from_root(source)?;
        let root_to_target = self.get_path_from_root(target)?;

        let cursor = zip(root_to_source.iter(), root_to_target.iter())
            .enumerate()
            .filter(|(_, (s, t))| s != t)
            .map(|(idx, _)| idx)
            .next()
            .unwrap_or_else(|| {
                // One node is a child of the other
                root_to_source.len().min(root_to_target.len())
            });

        Ok(root_to_source[cursor - 1])
    }

    /// Gets the distance between 2 nodes, returns the sum of branch lengths (if all
    /// branches in the path have lengths) and the number of edges in the path.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let (sum_edge_lengths, num_edges) = tree.get_distance(
    ///     tree.get_by_name("A").unwrap().id,
    ///     tree.get_by_name("I").unwrap().id,
    /// ).unwrap();
    ///
    /// assert_eq!(num_edges, 4);
    /// assert!(sum_edge_lengths.is_none());
    /// ```
    pub fn get_distance(
        &self,
        source: NodeId,
        target: NodeId,
    ) -> Result<(Option<f64>, usize), TreeError> {
        let mut dist = 0.0;
        let mut branches = 0;
        let mut all_dists = true;

        if source == target {
            return Ok((Some(0.0), 0));
        }

        let root_to_source = self.get_path_from_root(source)?;
        let root_to_target = self.get_path_from_root(target)?;

        let cursor = zip(root_to_source.iter(), root_to_target.iter())
            .enumerate()
            .filter(|(_, (s, t))| s != t)
            .map(|(idx, _)| idx)
            .next()
            .unwrap_or_else(|| {
                // One node is a child of the other
                root_to_source.len().min(root_to_target.len())
            });

        for list in [root_to_source, root_to_target] {
            for &node in list.iter().skip(cursor) {
                if let Some(d) = self.get(node)?.parent_edge {
                    dist += d;
                } else {
                    all_dists = false;
                }
                branches += 1;
            }
        }

        if all_dists {
            Ok((Some(dist), branches))
        } else {
            Ok((None, branches))
        }
    }

    // Implementation of recursive distance matrix computation
    fn distance_matrix_recursive_impl(
        &self,
        current: NodeId,
        prev: Option<NodeId>,
        lengths: &mut [f64],
        currlength: f64,
    ) -> Result<(), TreeError> {
        if prev.is_some() && self.get(current)?.is_tip() {
            lengths[current.index()] = currlength;
            return Ok(());
        }

        let children = self.get(current)?.children.clone();
        let mut neighbors: Vec<_> = children
            .iter()
            .map(|&idx| (idx, self.get(idx).unwrap().parent_edge))
            .collect();

        if let Some(parent) = self.get(current)?.parent {
            neighbors.push((parent, self.get(current).unwrap().parent_edge))
        }

        for (neighbor, brlen) in neighbors {
            if Some(neighbor) != prev {
                if let Some(brlen) = brlen {
                    self.distance_matrix_recursive_impl(
                        neighbor,
                        Some(current),
                        lengths,
                        currlength + brlen,
                    )?
                } else {
                    return Err(TreeError::MissingBranchLengths);
                }
            }
        }

        Ok(())
    }

    /// Computes the distance matrix of the tree, implemented in a recursive manner.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((T3:0.2,T1:0.2):0.3,(T2:0.4,T0:0.5):0.6);").unwrap();
    /// let matrix = tree.distance_matrix_recursive().unwrap();
    ///
    /// let phylip="\
    /// 4
    /// T0    0  1.6  0.9  1.6
    /// T1    1.6  0  1.5  0.4
    /// T2    0.9  1.5  0  1.5
    /// T3    1.6  0.4  1.5  0
    /// ";
    ///
    /// assert_eq!(phylip, matrix.to_phylip(true).unwrap())
    /// ```
    pub fn distance_matrix_recursive(&self) -> Result<DistanceMatrix<EdgeLength>, TreeError> {
        let size = self.nodes.len();
        let taxon_index = self.taxon_index()?;
        let mut matrix = DistanceMatrix::zeros(taxon_index.taxa().to_vec())?;
        let mut cache: Vec<Vec<_>> = vec![vec![f64::INFINITY; size]; size];

        for &tip in self.leaves().iter() {
            self.distance_matrix_recursive_impl(tip, None, &mut cache[tip.index()], 0.0)?
        }

        for pair in self.leaves().iter().combinations(2) {
            let (&i1, &i2) = (pair[0], pair[1]);
            let d = cache[i1.index()][i2.index()];
            let name1 = self.get(i1)?.name.clone().unwrap();
            let name2 = self.get(i2)?.name.clone().unwrap();

            matrix.set(&name1, &name2, d)?;
        }

        Ok(matrix)
    }

    /// Computes the distance matrix of the tree.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((T3:0.2,T1:0.2):0.3,(T2:0.4,T0:0.5):0.6);").unwrap();
    /// let matrix = tree.distance_matrix().unwrap();
    ///
    /// let phylip="\
    /// 4
    /// T0    0  1.6  0.9  1.6
    /// T1    1.6  0  1.5  0.4
    /// T2    0.9  1.5  0  1.5
    /// T3    1.6  0.4  1.5  0
    /// ";
    ///
    /// assert_eq!(phylip, matrix.to_phylip(true).unwrap())
    /// ```
    pub fn distance_matrix(&self) -> Result<DistanceMatrix<f64>, TreeError> {
        let mut leaf_order = self.leaves();
        leaf_order.sort_by(|&a, &b| self.get(a).unwrap().name.cmp(&self.get(b).unwrap().name));

        let n = self.n_leaves();
        let n_pairs = DistanceMatrix::<f64>::pair_count(n)?;
        let mut pairwise_vec = vec![NaiveSum::zero(); n_pairs];

        let leaf_idx_to_leaf_order = self
            .nodes
            .iter()
            .enumerate()
            .map(|(idx, _)| leaf_order.iter().position(|&v| v == NodeId(idx)))
            .collect_vec();

        // Converts the node index of a leaf to its index in the leaf_order array
        let get_leaf_index = |leaf: NodeId| -> Result<usize, TreeError> {
            leaf_idx_to_leaf_order[leaf.index()].ok_or(TreeError::NodeNotFound(leaf))
        };

        for &current_node in self.levelorder(self.root()?)?.iter().rev() {
            let mut node_cache = VecMap::new();

            let parent = self.get(current_node)?;
            if parent.is_tip() {
                node_cache.insert(current_node.index(), 0.);
            }

            // Compute distances from current node to descendant leaves
            for &child in parent.children.iter() {
                let child = self.get(child)?;

                // Use topological distance if no edge length
                let child_len = child.parent_edge.unwrap_or(1.0);

                for (leaf, distance) in child
                    .subtree_distances
                    .borrow()
                    .as_ref()
                    .ok_or(TreeError::MissingBranchLengths)?
                    .iter()
                {
                    let len = child_len + distance;
                    node_cache.insert(leaf, len);
                }
            }

            // Compute distances between leaves
            for subtree_roots in parent.children.iter().combinations(2) {
                let subtree1 = self.get(*subtree_roots[0])?;
                let subtree2 = self.get(*subtree_roots[1])?;

                for leaf1 in subtree1
                    .subtree_distances
                    .borrow()
                    .as_ref()
                    .ok_or(TreeError::MissingBranchLengths)?
                    .keys()
                {
                    for leaf2 in subtree2
                        .subtree_distances
                        .borrow()
                        .as_ref()
                        .ok_or(TreeError::MissingBranchLengths)?
                        .keys()
                    {
                        let distance1 = node_cache.get(leaf1).unwrap();
                        let distance2 = node_cache.get(leaf2).unwrap();

                        let mut i = get_leaf_index(NodeId(leaf1))?;
                        let mut j = get_leaf_index(NodeId(leaf2))?;
                        if j < i {
                            std::mem::swap(&mut i, &mut j);
                        }
                        // Compute the index of the pair in the vector representing
                        // the upper triangular matrix
                        //let vec_idx = ((2 * n - 3 - i) * i) / 2 + j - 1;
                        let vec_idx = tril_to_rowvec_index(n, i, j);

                        pairwise_vec[vec_idx] += distance1 + distance2;
                    }
                }
            }

            // Save distance between current node and descendant leaves
            (*(self.get(current_node)?).subtree_distances.borrow_mut()) = Some(node_cache);
        }

        let matrix = DistanceMatrix::from_precomputed(
            leaf_order
                .iter()
                .map(|&i| self.get(i).unwrap().clone().name.unwrap())
                .collect_vec(),
            pairwise_vec.iter().map(|v| v.sum()).collect_vec(),
        );

        Ok(matrix?)
    }
}

/// Methods to manipulate and alter the [`Tree`] object.
///
/// ----
/// ----
impl Tree {
    // ##################
    // # ALTER THE TREE #
    // ##################

    /// Prune the subtree starting at a given root node.
    /// # Example
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let mut tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let root_idx = tree.get_by_name("G").unwrap().id;
    ///
    /// tree.prune(root_idx);
    ///
    /// assert_eq!(tree.to_newick().unwrap(), "((A,(C,E)D)B)F;")
    /// ```
    pub fn prune(&mut self, root: NodeId) -> Result<(), TreeError> {
        for child in self.get(root)?.children.clone() {
            self.prune(child)?
        }

        if let Some(parent) = self.get(root)?.parent {
            self.get_mut(parent)?.remove_child(root)?;
        }

        self.get_mut(root)?.delete();

        Ok(())
    }

    // Removes a single node
    fn compress_node(&mut self, id: NodeId) -> Result<(), TreeError> {
        let node = self.get(id)?;

        if node.parent.is_none() || node.children.len() != 1 {
            return Err(TreeError::CouldNotCompressNode(id));
        }

        let parent = node.parent.unwrap();
        let child = node.children[0];
        let to_remove = node.id;

        let parent_edge = node.parent_edge;
        let child_edge = node.get_child_edge(child);

        let new_edge = match (parent_edge, child_edge) {
            (Some(p), Some(c)) => Some(p + c),
            (None, None) => None,
            _ => return Err(TreeError::MissingBranchLengths),
        };

        self.get_mut(child)?.set_parent(parent, new_edge);
        self.get_mut(parent)?.add_child(child, new_edge);
        self.get_mut(parent)?.remove_child(to_remove)?;

        self.get_mut(to_remove)?.delete();

        Ok(())
    }

    /// Compress the tree (i.e. remove nodes with exactly 1 parent and 1 child and fuse branches together)
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let mut tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// // Compress F->G->I->H to F->H
    /// tree.compress().unwrap();
    ///
    /// assert_eq!(tree.to_newick().unwrap(), "((A,(C,E)D)B,H)F;")
    /// ```
    pub fn compress(&mut self) -> Result<(), TreeError> {
        let to_compress: Vec<_> = self
            .nodes
            .iter()
            .filter(|node| !node.deleted && node.parent.is_some() && node.children.len() == 1)
            .map(|node| node.id)
            .collect();

        for id in to_compress {
            self.compress_node(id)?;
        }

        Ok(())
    }

    /// Rescale the branch lenghts of the tree
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let mut tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// // Double all branch lengths
    /// tree.rescale(2.0);
    ///
    /// assert_eq!(
    ///     tree.to_newick().unwrap(),
    ///     "(A:0.2,B:0.4,(C:0.6,D:0.8)E:1)F;"
    /// )
    /// ```
    pub fn rescale(&mut self, factor: f64) {
        for node in self.nodes.iter_mut() {
            node.rescale_edges(factor);
            node.subtree_distances.borrow_mut().take();
        }
    }

    /// Rescale all branch lengths to a positive target tip-to-tip diameter.
    ///
    /// Requires finite, nonnegative lengths on every non-root edge and a
    /// positive current diameter. Invalid inputs leave the tree unchanged.
    /// The root stem does not contribute to the diameter.
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("(A:1,B:2);").unwrap();
    /// tree.rescale_diameter(6.0).unwrap();
    /// assert_eq!(tree.diameter().unwrap(), 6.0);
    /// ```
    pub fn rescale_diameter(&mut self, target: EdgeLength) -> Result<(), TreeError> {
        if !target.is_finite() || target <= 0.0 {
            return Err(TreeError::GeneralError("Positive finite diameter required"));
        }
        for node in self
            .nodes
            .iter()
            .filter(|n| !n.deleted && n.parent.is_some())
        {
            let length = node.parent_edge.ok_or(TreeError::MissingBranchLengths)?;
            if !length.is_finite() || length < 0.0 {
                return Err(TreeError::GeneralError(
                    "Finite nonnegative branch lengths required",
                ));
            }
        }
        let current = self.diameter()?;
        if current <= 0.0 {
            return Err(TreeError::GeneralError("Positive diameter required"));
        }
        let factor = target / current;
        if !factor.is_finite() || factor <= 0.0 {
            return Err(TreeError::GeneralError(
                "Diameter scale is not representable",
            ));
        }
        if self.nodes.iter().filter(|n| !n.deleted).any(|n| {
            n.parent_edge
                .is_some_and(|length| !(length * factor).is_finite())
        }) {
            return Err(TreeError::GeneralError(
                "Rescaled branch length is not representable",
            ));
        }
        self.rescale(factor);
        Ok(())
    }

    /// Randomly resolve multifurcations to binarize the tree
    ///
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let mut tree = Tree::from_newick("((A:0.1,B:0.2):0.3, (C:0.1,D:0.2,E:0.4)F:0.5)G;").unwrap();
    /// assert!(!tree.is_binary().unwrap());
    ///
    /// tree.resolve();
    /// assert!(tree.is_binary().unwrap());
    /// ```
    pub fn resolve(&mut self) -> Result<(), TreeError> {
        let rng = &mut rand::thread_rng();
        let mut to_binarize = vec![];
        for node in self.nodes.iter() {
            if node.children.len() > 2 {
                to_binarize.push(node.id);
            }
        }

        for &node_id in to_binarize.iter() {
            loop {
                let mut children = self.get(node_id)?.children.clone();
                children.shuffle(rng);

                let parent = self.add_child(Node::new(), node_id, Some(0.0))?;

                for _ in 0..2 {
                    let child = children.pop().unwrap();
                    let edge = self.get(child)?.parent_edge;
                    self.get_mut(parent)?.add_child(child, edge);
                    self.get_mut(child)?.set_parent(parent, edge);
                    self.get_mut(node_id)?.remove_child(child)?;
                }

                children.push(parent);

                if children.len() <= 2 {
                    break;
                }
            }
        }
        Ok(())
    }

    /// Sort children of a node by number of descendants
    ///
    /// ```
    ///use phylotree::tree::Tree;
    ///
    ///let mut tree = Tree::from_newick("(A,(((D,(E,F)),C),B));").unwrap();
    ///tree.ladderize();
    ///
    ///assert_eq!("(A,(B,(C,(D,(E,F)))));", tree.to_newick().unwrap());
    ///
    /// ```
    pub fn ladderize(&mut self) -> Result<(), TreeError> {
        let mut descendant_counter = vec![0; self.nodes.len()];
        let root = self.root()?;
        // Go from tips to root
        for node_id in self.levelorder(root)?.into_iter().rev() {
            let node = self.get_mut(node_id)?;
            for child in node.children.iter() {
                descendant_counter[node_id.index()] += descendant_counter[child.index()] + 1;
            }
            node.children
                .sort_by_key(|node_id| descendant_counter[node_id.index()]);
        }

        Ok(())
    }

    // recusrive implementation of depth recomputation
    fn reset_depth_impl(&mut self, root: NodeId, depth: usize) -> Result<(), TreeError> {
        let root = self.get_mut(root)?;
        root.set_depth(depth);

        for &child in root.children.clone().iter() {
            self.reset_depth_impl(child, depth + 1)?
        }

        Ok(())
    }

    /// Recompute node depths and set them correctly.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let mut tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// tree.reset_depths().unwrap();
    /// assert_eq!(tree.get(tree.root().unwrap()).unwrap().get_depth(), 0);
    /// ```
    pub fn reset_depths(&mut self) -> Result<(), TreeError> {
        let root = self.root()?;
        self.reset_depth_impl(root, 0)
    }

    /// Merge 2 sibling nodes into a new parent node.
    /// Useful for agglomerative tree building / polytomy resolution
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// // Initialize star tree
    /// let mut tree = Tree::from_newick("(A,B,C);").unwrap();
    /// let a = tree.get_by_name("A").unwrap().id;
    /// let b = tree.get_by_name("B").unwrap().id;
    ///
    /// // Merge A and B into node D
    /// tree.merge_children(a, b, None, None, None, Some("D".into()));
    ///
    /// let expected = Tree::from_newick("((A,B)D, C);").unwrap();
    /// assert_eq!(tree.unrooted_robinson_foulds(&expected).unwrap(), 0);
    /// ```
    pub fn merge_children(
        &mut self,
        child1: NodeId,
        child2: NodeId,
        edge1: Option<EdgeLength>,
        edge2: Option<EdgeLength>,
        parent_edge: Option<EdgeLength>,
        parent_name: Option<String>,
    ) -> Result<NodeId, TreeError> {
        // Check that nodes are siblings
        let parent = self.get(child1)?.parent;
        if parent != self.get(child2)?.parent {
            return Err(TreeError::MergingNonSiblingNodes(child1, child2));
        }

        // Add new parent node as child of current parent
        let parent = match parent {
            Some(parent_id) => {
                // Remove merged nodes as children of current parent
                let parent_node = self.get_mut(parent_id)?;
                parent_node.remove_child(child1)?;
                parent_node.remove_child(child2)?;
                // Add new parent
                self.add_child(Node::new(), parent_id, parent_edge)?
            }
            None => self.add(Node::new()),
        };

        // Set parent/child relationships between merged nodes and new parent node
        let p = self.get_mut(parent)?;
        p.add_child(child1, edge1);
        p.add_child(child2, edge2);
        p.name = parent_name;

        // Set new parent in child nodes
        self.get_mut(child1)?.set_parent(parent, edge1);
        self.get_mut(child2)?.set_parent(parent, edge2);

        Ok(parent)
    }
}

/// Methods to read and write [`Tree`] objects to and from files or [`String`] objects.
///
/// ----
/// ----
impl Tree {
    // ########################
    // # READ AND WRITE TREES #
    // ########################

    /// Writes the tree as a newick formatted string
    /// # Example
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F:0.6;";
    /// let tree = Tree::from_newick(newick).unwrap();
    ///
    /// assert_eq!(tree.to_newick().unwrap(), newick);
    /// ```
    pub fn to_newick(&self) -> Result<String, TreeError> {
        NewickSerializer::new(self, NewickFormat::AllFields.options()).serialize()
    }

    /// Writes the tree as a Newick-formatted string using a [`NewickFormat`]
    /// compatibility preset or custom serialization options.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::{Tree, NewickFormat};
    ///
    /// let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F:0.6;";
    /// let tree = Tree::from_newick(newick).unwrap();
    ///
    /// assert_eq!(tree.to_formatted_newick(NewickFormat::Topology).unwrap(), "(,,(,));");
    /// assert_eq!(
    ///     tree.to_formatted_newick(NewickFormat::OnlyNames).unwrap(),
    ///     "(A,B,(C,D)E)F;"
    /// );
    /// assert_eq!(
    ///     tree.to_formatted_newick(NewickFormat::InternalLengthsLeafNames).unwrap(),
    ///     "(A,B,(C,D):0.5):0.6;"
    /// );
    /// ```
    ///
    /// Custom options can independently select attributes and support output:
    /// ```
    /// use phylotree::tree::{
    ///     InternalNodeLabelMode, NewickFormat, NewickParseOptions,
    ///     NewickSerializeOptions, NewickSupportFormat, Tree,
    /// };
    ///
    /// let tree = Tree::from_newick_with_options(
    ///     "(A,B)95;",
    ///     NewickParseOptions {
    ///         internal_node_labels: InternalNodeLabelMode::NumericSupport,
    ///         ..Default::default()
    ///     },
    /// ).unwrap();
    /// let format = NewickFormat::Custom(NewickSerializeOptions {
    ///     support: NewickSupportFormat::InternalNodeLabel,
    ///     ..Default::default()
    /// });
    ///
    /// assert_eq!(tree.to_formatted_newick(format).unwrap(), "(A,B)95;");
    /// ```
    pub fn to_formatted_newick(&self, format: NewickFormat) -> Result<String, TreeError> {
        NewickSerializer::new(self, format.options()).serialize()
    }

    /// Read a newick formatted string and build a [`Tree`] struct from it.
    /// # Example
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;";
    /// let tree = Tree::from_newick(newick).unwrap();
    ///
    /// assert_eq!(tree.size(), 6);
    /// assert_eq!(tree.n_leaves(), 4);
    /// assert_eq!(tree.is_rooted().unwrap(), false);
    /// ```
    pub fn from_newick(newick: &str) -> Result<Self, NewickParseError> {
        let mut tokenizer = NewickTokenizer::new(newick.as_bytes());
        let nodes = NewickParser::new(&mut tokenizer).parse()?;
        Ok(Self::from_nodes(nodes))
    }

    /// Reads a Newick-formatted string using explicit parsing options.
    ///
    /// Unlike [`Tree::from_newick`], this method can interpret convention-based
    /// internal labels and comments as branch support values.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::{InternalNodeLabelMode, NewickParseOptions, Tree};
    ///
    /// let options = NewickParseOptions {
    ///     internal_node_labels: InternalNodeLabelMode::NumericSupport,
    ///     ..Default::default()
    /// };
    /// let tree = Tree::from_newick_with_options("(A,B)95;", options).unwrap();
    /// let root = tree.root().unwrap();
    ///
    /// assert_eq!(tree.get(root).unwrap().support, Some(95.0));
    /// ```
    pub fn from_newick_with_options(
        newick: &str,
        options: NewickParseOptions,
    ) -> Result<Self, NewickParseError> {
        let mut tokenizer = NewickTokenizer::new(newick.as_bytes());
        let nodes = NewickParser::with_options(&mut tokenizer, options).parse()?;
        Ok(Self::from_nodes(nodes))
    }

    /// Writes the tree to a newick file
    ///
    /// ```no_run
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// tree.to_file(std::path::Path::new("tree.nwk")).unwrap();
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
    pub fn to_file(&self, path: &Path) -> Result<(), TreeError> {
        match fs::write(path, self.to_newick()?) {
            Ok(_) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    /// Creates a tree from a newick file
    ///
    /// ```no_run
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_file(std::path::Path::new("tree.nwk")).unwrap();
    /// println!("{} tips", tree.n_leaves());
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &Path) -> Result<Self, NewickParseError> {
        Self::from_file_with_options(path, NewickParseOptions::default())
    }

    /// Creates a tree from a Newick file using explicit parsing options.
    ///
    /// ```no_run
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_file_with_options(std::path::Path::new("tree.nwk"), Default::default()).unwrap();
    /// println!("{} tips", tree.n_leaves());
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file_with_options(
        path: &Path,
        options: NewickParseOptions,
    ) -> Result<Self, NewickParseError> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_newick_with_options(&newick_string, options)
    }

    /// Outputs a Nexus formatted string of the tree
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    /// assert!(tree.to_nexus().unwrap().contains("#NEXUS"));
    /// ```
    pub fn to_nexus(&self) -> Result<String, TreeError> {
        let nwk = self.to_newick()?;
        let n = self.n_leaves();
        let labels = self
            .nodes
            .iter()
            .filter_map(|node| {
                if node.is_tip() {
                    node.name.clone()
                } else {
                    None
                }
            })
            .join(" ");

        Ok(format!(
            "#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX={n};
    TAXLABELS {labels};
END;
BEGIN TREES;
    TREE tree1 = {nwk}
END;
"
        ))
    }

    /// Recursive function that adds node representation to a printable tree builder
    #[cfg(not(target_arch = "wasm32"))]
    fn print_nodes(
        &self,
        root_idx: NodeId,
        output_tree: &mut TreeBuilder,
        debug: bool,
    ) -> Result<(), TreeError> {
        let root = self.get(root_idx)?;
        let label = if debug {
            format!("{root:?}")
        } else {
            format!("{root}")
        };

        if root.children.is_empty() {
            output_tree.add_empty_child(label);
        } else {
            output_tree.begin_child(label);
            for &child_idx in root.children.iter() {
                self.print_nodes(child_idx, output_tree, debug)?;
            }
            output_tree.end_child();
        }

        Ok(())
    }

    /// Print a debug view of the tree to the console
    #[cfg(not(target_arch = "wasm32"))]
    pub fn print_debug(&self) -> Result<(), TreeError> {
        let root = self.root()?;
        let mut builder = TreeBuilder::new(format!("{:?}", root));
        for &child_idx in self.get(root)?.children.iter() {
            self.print_nodes(child_idx, &mut builder, true)?;
        }
        let tree = builder.build();
        print_tree(&tree)?;
        Ok(())
    }

    /// Print the tree to the console
    #[cfg(not(target_arch = "wasm32"))]
    pub fn print(&self) -> Result<(), TreeError> {
        let root = self.root()?;
        let mut builder = TreeBuilder::new(format!("{:?}", root));
        for &child_idx in self.get(root)?.children.iter() {
            self.print_nodes(child_idx, &mut builder, false)?;
        }
        let tree = builder.build();
        print_tree(&tree)?;
        Ok(())
    }
}

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
// #[allow(clippy::excessive_precision)]
mod tests {

    use core::f64;

    use super::*;

    /// Generates example tree from the tree traversal wikipedia page
    /// https://en.wikipedia.org/wiki/Tree_traversal#Depth-first_search
    /// The difference is that I is the left child of G since this tree structure
    /// cannot represent a right child only.
    fn build_simple_tree() -> Result<Tree, TreeError> {
        let mut tree = Tree::new();
        tree.add(Node::new_named("F")); // 0
        tree.add_child(Node::new_named("B"), NodeId(0), None)?; // 1
        tree.add_child(Node::new_named("G"), NodeId(0), None)?; // 2
        tree.add_child(Node::new_named("A"), NodeId(1), None)?; // 3
        tree.add_child(Node::new_named("D"), NodeId(1), None)?; // 4
        tree.add_child(Node::new_named("I"), NodeId(2), None)?; // 5
        tree.add_child(Node::new_named("C"), NodeId(4), None)?; // 6
        tree.add_child(Node::new_named("E"), NodeId(4), None)?; // 7
        tree.add_child(Node::new_named("H"), NodeId(5), None)?; // 8

        Ok(tree)
    }

    /// Generates example tree from the newick format wikipedia page
    /// https://en.wikipedia.org/wiki/Newick_format#Examples
    fn build_tree_with_lengths() -> Result<Tree, TreeError> {
        let mut tree = Tree::new();
        tree.add(Node::new_named("F")); // 0
        tree.add_child(Node::new_named("A"), NodeId(0), Some(0.1))?; // 1
        tree.add_child(Node::new_named("B"), NodeId(0), Some(0.2))?; // 2
        tree.add_child(Node::new_named("E"), NodeId(0), Some(0.5))?; // 3
        tree.add_child(Node::new_named("C"), NodeId(3), Some(0.3))?; // 4
        tree.add_child(Node::new_named("D"), NodeId(3), Some(0.4))?; // 5

        Ok(tree)
    }

    fn build_tree_without_lengths() -> Result<Tree, TreeError> {
        let mut tree = Tree::new();
        tree.add(Node::new_named("F")); // 0
        tree.add_child(Node::new_named("A"), NodeId(0), None)?; // 1
        tree.add_child(Node::new_named("B"), NodeId(0), None)?; // 2
        tree.add_child(Node::new_named("E"), NodeId(0), None)?; // 3
        tree.add_child(Node::new_named("C"), NodeId(3), None)?; // 4
        tree.add_child(Node::new_named("D"), NodeId(3), None)?; // 5

        Ok(tree)
    }
    fn get_values(indices: &[NodeId], tree: &Tree) -> Vec<Option<String>> {
        indices
            .iter()
            .map(|&idx| tree.get(idx).unwrap().name.clone())
            .collect()
    }

    fn node_ids_from_usize(indices: &[usize]) -> Vec<NodeId> {
        indices.iter().map(|&v| NodeId(v)).collect()
    }

    #[test]
    fn test_tips() {
        let mut tree = Tree::new();
        tree.add(Node::new_named("root"));
        assert_eq!(tree.leaves(), vec![NodeId(0)]);

        tree.add_child(Node::new_named("A"), NodeId(0), Some(0.1))
            .unwrap(); // 1
        tree.add_child(Node::new_named("B"), NodeId(0), Some(0.2))
            .unwrap(); // 2
        tree.add_child(Node::new_named("E"), NodeId(0), Some(0.5))
            .unwrap(); // 3

        assert_eq!(tree.leaves(), node_ids_from_usize(&[1, 2, 3]));

        tree.add_child(Node::new_named("C"), NodeId(3), Some(0.3))
            .unwrap(); // 4
        tree.add_child(Node::new_named("D"), NodeId(3), Some(0.4))
            .unwrap(); // 5

        assert_eq!(tree.leaves(), node_ids_from_usize(&[1, 2, 4, 5]));
    }

    #[test]
    fn test_binary() {
        let mut tree = Tree::new();
        tree.add(Node::new_named("root"));

        tree.add_child(Node::new_named("0L"), NodeId(0), None)
            .unwrap(); //1
        tree.add_child(Node::new_named("0R"), NodeId(0), None)
            .unwrap(); //2

        assert!(tree.is_binary().unwrap());

        tree.add_child(Node::new_named("1L"), NodeId(1), None)
            .unwrap(); //3
        tree.add_child(Node::new_named("1R"), NodeId(1), None)
            .unwrap(); //4

        assert!(tree.is_binary().unwrap());

        tree.add_child(Node::new_named("3L"), NodeId(3), None)
            .unwrap(); //5
        tree.add_child(Node::new_named("3R"), NodeId(3), None)
            .unwrap(); //6
        assert!(tree.is_binary().unwrap());

        tree.add_child(Node::new_named("3?"), NodeId(3), None)
            .unwrap(); //7
        assert!(!tree.is_binary().unwrap());
    }

    #[test]
    fn binary_from_newick() {
        let test_cases = vec![
            ("((A,B,C)D,E)F;", false),   // Rooted non binary
            ("(A,B,(C,D)E)F;", true),    // Unrooted binary
            ("((D,E)B,(F,G)C)A;", true), // rooted binary
        ];

        for (newick, is_binary) in test_cases {
            assert_eq!(
                Tree::from_newick(newick).unwrap().is_binary().unwrap(),
                is_binary
            )
        }
    }

    #[test]
    fn prune_tree() {
        let mut tree = build_simple_tree().unwrap();
        tree.prune(NodeId(4)).unwrap(); // prune D subtree

        assert_eq!(tree.to_newick().unwrap(), "((A)B,((H)I)G)F;");
    }

    #[test]
    fn path_from_root() {
        let tree = build_simple_tree().unwrap();
        let values: Vec<_> = get_values(&(tree.get_path_from_root(NodeId(7)).unwrap()), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["F", "B", "D", "E"])
    }

    #[test]
    fn last_common_ancestor() {
        let test_cases = vec![
            ((NodeId(3), NodeId(7)), NodeId(1)), // (A,E) -> B
            ((NodeId(6), NodeId(8)), NodeId(0)), // (C,H) -> F
            ((NodeId(3), NodeId(3)), NodeId(3)), // (A,A) -> A
            ((NodeId(8), NodeId(5)), NodeId(5)), // (H,I) -> I
            ((NodeId(4), NodeId(7)), NodeId(4)), // (D,E) -> D
        ];
        let tree = build_simple_tree().unwrap();
        for ((source, target), ancestor) in test_cases {
            println!(
                "Testing: ({:?}, {:?}) -> {:?}",
                tree.get(source).unwrap().name,
                tree.get(target).unwrap().name,
                tree.get(ancestor).unwrap().name
            );
            assert_eq!(ancestor, tree.get_common_ancestor(source, target).unwrap());
        }
    }

    #[test]
    fn get_distances_lengths() {
        let test_cases = vec![
            ((NodeId(1), NodeId(3)), (0.6, 2)), // (A,E)
            ((NodeId(1), NodeId(4)), (0.9, 3)), // (A,C)
            ((NodeId(4), NodeId(5)), (0.7, 2)), // (C,D)
            ((NodeId(5), NodeId(2)), (1.1, 3)), // (D,B)
            ((NodeId(2), NodeId(5)), (1.1, 3)), // (B,D)
            ((NodeId(0), NodeId(2)), (0.2, 1)), // (F,B)
            ((NodeId(1), NodeId(1)), (0.0, 0)), // (A,A)
        ];
        let tree = build_tree_with_lengths().unwrap();
        let tree2 = build_tree_without_lengths().unwrap();

        for ((idx_s, idx_t), (dist, branches)) in test_cases {
            let (d_pred, b_pred) = tree.get_distance(idx_s, idx_t).unwrap();
            let (d_pred2, b_pred2) = tree2.get_distance(idx_s, idx_t).unwrap();
            assert_eq!(branches, b_pred);
            assert_eq!(branches, b_pred2);

            assert!(d_pred.is_some());
            assert!((d_pred.unwrap() - dist).abs() < f64::EPSILON);
            if idx_s != idx_t {
                assert!(d_pred2.is_none(), "{d_pred2:?}");
            }
        }
    }

    #[test]
    fn get_correct_leaves() {
        let tree = build_simple_tree().unwrap();
        let values: Vec<_> = get_values(&(tree.leaves()), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["A", "C", "E", "H"])
    }

    #[test]
    fn to_newick() {
        let tree = build_tree_with_lengths().unwrap();
        dbg!(&tree);
        assert_eq!(
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            tree.to_newick().unwrap()
        );
    }

    #[test]
    fn to_formatted_newick() {
        let newick = "(A:0.1[Comment_1],B:0.2,(C:0.3,D:0.4)E:0.5[Comment_2])F;";
        let tree = Tree::from_newick(newick).unwrap();

        let cases = [
            (NewickFormat::AllFields, newick),
            (
                NewickFormat::NoComments,
                "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            ),
            (NewickFormat::Topology, "(,,(,));"),
            (NewickFormat::OnlyNames, "(A,B,(C,D)E)F;"),
            (NewickFormat::OnlyLengths, "(:0.1,:0.2,(:0.3,:0.4):0.5);"),
            (
                NewickFormat::LeafLengthsAllNames,
                "(A:0.1,B:0.2,(C:0.3,D:0.4)E)F;",
            ),
            (
                NewickFormat::LeafLengthsLeafNames,
                "(A:0.1,B:0.2,(C:0.3,D:0.4));",
            ),
            (NewickFormat::InternalLengthsLeafNames, "(A,B,(C,D):0.5);"),
            (
                NewickFormat::AllLengthsLeafNames,
                "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
            ),
        ];

        for (format, expected) in cases {
            assert_eq!(
                expected,
                tree.to_formatted_newick(format.clone()).unwrap(),
                "Failed to write newick for format: {format:?}"
            )
        }
    }

    // test cases from https://github.com/ila/Newick-validator
    #[test]
    fn read_newick() {
        let newick_strings = vec![
            "((D,E)B,(F,G)C)A;",
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
            "(dog:20,(elephant:30,horse:60):20):50;",
            "(A,B,(C,D));",
            "(A,B,(C,D)E)F;",
            "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0;",
            "(:0.1,:0.2,(:0.3,:0.4):0.5);",
            "(:0.1,:0.2,(:0.3,:0.4):0.5):0;",
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;",
            "(,,(,));",
            "('hungarian dog':20,('indian elephant':30,'swedish horse':60):20):50;",
            "('hungarian dog':20[Comment_1],('indian elephant':30,'swedish horse':60[Another interesting comment]):20):50;",
        ];
        for newick in newick_strings {
            let tree = Tree::from_newick(newick).unwrap();
            assert_eq!(newick, tree.to_newick().unwrap());
        }
    }

    #[test]
    fn reads_branch_support_with_explicit_options() {
        let options = NewickParseOptions {
            internal_node_labels: crate::tree::InternalNodeLabelMode::NumericSupport,
            support_comments: crate::tree::SupportCommentMode::Nhx,
        };
        let tree = Tree::from_newick_with_options("((A,B)95,C)[&&NHX:B=99];", options).unwrap();
        let root_id = tree.root().unwrap();
        let root = tree.get(root_id).unwrap();
        let inner = tree.get(root.children[0]).unwrap();

        assert_eq!(root.support, Some(99.0));
        assert_eq!(inner.support, Some(95.0));
    }

    #[ignore]
    #[test]
    fn to_nexus() {
        let tree = crate::generate_tree(10, true, crate::distr::Distr::Uniform).unwrap();
        println!("{}", tree.to_nexus().unwrap());
        panic!()
    }

    #[test]
    fn read_newick_fails() {
        assert!(matches!(
            Tree::from_newick("((D,E)B,(F,G,C)A;"),
            Err(NewickParseError::UnexpectedToken { .. })
        ));
        assert!(matches!(
            Tree::from_newick("((D,E)B,(F,G)C)A"),
            Err(NewickParseError::UnexpectedEOF)
        ));
    }

    #[test]
    fn test_subtree_leaves() {
        let test_cases = vec![
            (
                "((T0,T1)I1,(T2,T3)I2,((T4,T5)I4,(T6,T7)I4)I3)I0;",
                NodeId(7),
                vec!["T4", "T5", "T6", "T7"],
            ),
            (
                "((T0,T1)I1,(T2,T3)I2,((T4,T5)I4,(T6,T7)I4)I3)I0;",
                NodeId(1),
                vec!["T0", "T1"],
            ),
            (
                "((((((((T9,T8)I7,T7)I6,T6)I5,T5)I4,T4)I3,T3)I2,T2)I1,T0,T1)I0;",
                NodeId(5),
                vec!["T9", "T8", "T7", "T6"],
            ),
        ];

        for (newick, root, leaves) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let sub_leaves: Vec<String> = get_values(&tree.subtree_leaves(root).unwrap(), &tree)
                .iter()
                .map(|v| v.clone().unwrap())
                .collect();

            assert_eq!(sub_leaves, leaves);
        }
    }

    #[test]
    fn test_height() {
        // heights computed with ete3
        let test_cases = vec![
            ("((A:0.1,B:0.2)G:0.1,(C:0.3,D:0.4)E:0.5)F;", 0.9),
            ("((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1,G:0.1)F;", 1.0),
            ("((A,B)G,(C,D)E)F;", 2.0),
            (
                "((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1);",
                8.0,
            ),
        ];

        for (newick, height) in test_cases {
            assert_eq!(Tree::from_newick(newick).unwrap().height().unwrap(), height)
        }
    }

    #[test]
    fn test_diam() {
        let test_cases = vec![
            ("((D,E)B,(F,G)C)A;", 4.0),
            ("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;", 1.1),
            ("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);", 1.1),
            ("(A,B,(C,D));", 3.0),
            ("(A,B,(C,D)E)F;", 3.0),
        ];

        for (newick, diameter) in test_cases {
            assert_eq!(
                Tree::from_newick(newick).unwrap().diameter().unwrap(),
                diameter
            )
        }
    }

    #[test]
    fn test_cherries() {
        // Number of cherries computed with gotree
        let test_cases = vec![
            ("(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);", 1),
            ("(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);", 5),
            ("((a:0.2,b:0.2):0.2,((c:0.2,d:0.2):0.2,((e:0.2,f:0.2):0.2,((g:0.2,h:0.2):0.2,(i:0.2,j:0.2):0.2):0.2):0.2):0.2);", 5),
            ("(((d:0.3,e:0.3):0.3,((f:0.3,g:0.3):0.3,(h:0.3,(i:0.3,j:0.3):0.3):0.3):0.3):0.3,(a:0.3,(b:0.3,c:0.3):0.3):0.3);", 4),
        ];

        for (newick, true_cherries) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let cherries = tree.n_cherries();
            assert!(cherries.is_ok());
            assert_eq!(cherries.unwrap(), true_cherries);
        }
    }

    #[test]
    fn manual_colless() {
        let newick = "(((((((((T8,T9)I8,T7)I7,T6)I6,T5)I5,T4)I4,T3)I3,T2)I2,T1)I1,T0)I0;";
        let tree = Tree::from_newick(newick).unwrap();

        let mut colless = 0;

        dbg!(&tree);

        for node in tree.nodes.iter().filter(|node| !node.is_tip()) {
            let left = tree.subtree_leaves(node.children[0]).unwrap();
            let right = tree.subtree_leaves(node.children[1]).unwrap();
            eprintln!("Node {:?}:", node.name.clone().unwrap());
            eprintln!(
                "\tLeft: {:?}",
                get_values(&left, &tree)
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>()
            );
            eprintln!(
                "\tRight: {:?}",
                get_values(&right, &tree)
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>()
            );
            colless += left.len().abs_diff(right.len());
            eprintln!("Colless: {colless}\n")
        }

        assert_eq!(colless, 36)
    }

    #[test]
    fn test_colless_rooted() {
        // Colless index computed with gotree
        let test_cases = vec![
            ("(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);", 36),
            ("(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);", 4),
            ("((a:0.2,b:0.2):0.2,((c:0.2,d:0.2):0.2,((e:0.2,f:0.2):0.2,((g:0.2,h:0.2):0.2,(i:0.2,j:0.2):0.2):0.2):0.2):0.2);", 12),
            ("(((d:0.3,e:0.3):0.3,((f:0.3,g:0.3):0.3,(h:0.3,(i:0.3,j:0.3):0.3):0.3):0.3):0.3,(a:0.3,(b:0.3,c:0.3):0.3):0.3);", 10),
        ];

        for (newick, true_colless) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let colless = tree.colless();
            assert!(colless.is_ok());
            if tree.colless().unwrap() != true_colless {
                panic!(
                    "Computed colless {} not equal to true {true_colless}",
                    tree.colless().unwrap()
                )
            };
        }
    }

    #[test]
    fn test_sackin_rooted() {
        // Sackin index computed with gotree
        let test_cases = vec![
            ("(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);", 54),
            ("(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);", 34),
            ("((a:0.2,b:0.2):0.2,((c:0.2,d:0.2):0.2,((e:0.2,f:0.2):0.2,((g:0.2,h:0.2):0.2,(i:0.2,j:0.2):0.2):0.2):0.2):0.2);", 38),
            ("(((d:0.3,e:0.3):0.3,((f:0.3,g:0.3):0.3,(h:0.3,(i:0.3,j:0.3):0.3):0.3):0.3):0.3,(a:0.3,(b:0.3,c:0.3):0.3):0.3);", 36),
        ];

        for (newick, true_sackin) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let sackin = tree.sackin();
            assert!(sackin.is_ok());
            assert_eq!(tree.sackin().unwrap(), true_sackin);
        }
    }

    #[test]
    fn test_sackin_unrooted() {
        let test_cases = vec![
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);",
            "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1);",
            "(A,B,(C,D)E);",
            "((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip0,Tip1);",
        ];

        for newick in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            assert!(tree.sackin().is_err())
        }
    }

    #[test]
    fn test_rescale() {
        let test_cases = [
            ("((D:0.05307533041908017723553570021977,(C:0.08550401213833067060043902074540,(B:0.27463239708134284944307523801399,A:0.37113575171985613287972682883264)1:0.18134330279626256765546088445262)1:0.08033066840794983454188127325324)1:0.13864016688124142229199264875206,E:0.05060148260657528623829293223935);",
            "((D:0.04212872094323715649322181775460,(C:0.06786909546224775824363462106703,(B:0.21799038323938338401752901063446,A:0.29459024358034957558061250892933)1:0.14394185279875840177687962295749)1:0.06376273658252405718283029045779)1:0.11004609591585229333432494058798,E:0.04016509597234880352134567260691);",
            0.6525060248498331),
            ("(E:0.01699652764738122934229380689430,(D:0.00408169520164380558724381842239,(C:0.19713461567160570075962766622979,(B:0.12068059163592816107613003850929,A:0.45190753170439451613660253315174)1:0.03279750996120785189180679708443)1:0.21625179801434316062547225101298)1:0.03998705111996220251668887613050);",
            "(E:0.01986870266959113798255209815125,(D:0.00477144449924469995355513773916,(C:0.23044760352958004734347241537762,(B:0.14107392068250154681940955470054,A:0.52827357257097584675165080625447)1:0.03833982959587604877338407050047)1:0.25279532182407132845369801543711)1:0.04674430247278672095889717752470);",
            0.8860217291333011),
            ("((C:0.20738366520293352590620372666308,(B:0.19695170474498663315543467433599,A:0.02188551422116874478618342436675)1:0.05940680521299050026451382677806)1:0.13029006694844610936279138968530,(E:0.17189347707484656235799036494427,D:0.05867747522240193691622778260353)1:0.08673941227771603257323818070290);",
            "((C:0.18371634870356487456710681271943,(B:0.17447491841406459478491797199240,A:0.01938786624432843955223582099734)1:0.05262710219338979922287791168856)1:0.11542092936147484161235610145013,(E:0.15227641937588842768747099398752,D:0.05198100577716616849111019860175)1:0.07684042085359836515845444182560);",
            0.571639790198416),
        ];

        for (orig, rescaled, scale) in test_cases {
            let mut tree = Tree::from_newick(orig).unwrap();
            let rescaled = Tree::from_newick(rescaled).unwrap();

            let diam = tree.diameter().unwrap();

            tree.rescale(scale / diam);

            println!("Dealing with tree: {} and scale {}", orig, scale);
            for (n1, n2) in zip(tree.nodes, rescaled.nodes) {
                assert_eq!(n1, n2)
            }
        }
    }

    #[test]
    fn test_unique_tip_names() {
        let test_cases = vec![
            ("(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);",true),
            ("(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);", true),
            ("(((((((((Tip8,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);",false),
        ];

        for (newick, is_unique) in test_cases {
            assert_eq!(
                Tree::from_newick(newick)
                    .unwrap()
                    .has_unique_tip_names()
                    .unwrap(),
                is_unique,
                "Failed on: {newick}"
            )
        }

        assert!(Tree::from_newick("(((((((((,),),),),),),),),);")
            .unwrap()
            .has_unique_tip_names()
            .is_err())
    }

    #[test]
    fn test_descendants() {
        let tree = build_simple_tree().unwrap();
        eprintln!("{}", tree.to_newick().unwrap());
        dbg!(&tree);
        let descendants_b: Vec<_> = get_values(&tree.descendants(NodeId(1)).unwrap(), &tree)
            .into_iter()
            .flatten()
            .sorted()
            .collect();
        let descendants_g: Vec<_> = get_values(&tree.descendants(NodeId(2)).unwrap(), &tree)
            .into_iter()
            .flatten()
            .sorted()
            .collect();

        assert_eq!(descendants_b, vec!["A", "C", "D", "E"]);
        assert_eq!(descendants_g, vec!["H", "I"]);
    }

    #[test]
    fn test_compress() {
        let mut tree = Tree::new();
        tree.add(Node::new_named("root"));
        tree.add_child(Node::new_named("tip_A"), NodeId(0), Some(1.0))
            .unwrap();
        tree.add_child(Node::new_named("in_B"), NodeId(0), Some(1.0))
            .unwrap();
        tree.add_child(Node::new_named("in_C"), NodeId(2), Some(1.0))
            .unwrap();
        tree.add_child(Node::new_named("tip_D"), NodeId(3), Some(1.0))
            .unwrap();

        tree.compress().unwrap();

        assert_eq!(tree.to_newick().unwrap(), "(tip_A:1,tip_D:3)root;");
    }

    #[test]
    fn test_get_partitions() {
        let test_cases = vec![
            (
                "(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);",
                "(Tip0,((Tip2,(Tip3,(Tip4,(Tip5,(Tip6,((Tip8,Tip9),Tip7)))))),Tip1));",
            ),
            (
                "(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);",
                "(((c:0.1,d:0.1):0.1,((g:0.1,h:0.1):0.1,(f:0.1,e:0.1):0.1):0.1):0.1,((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1);",
            ),
            (
                "((a:0.2,b:0.2):0.2,((c:0.2,d:0.2):0.2,((e:0.2,f:0.2):0.2,((g:0.2,h:0.2):0.2,(i:0.2,j:0.2):0.2):0.2):0.2):0.2);",
                "((((e:0.2,f:0.2):0.2,((i:0.2,j:0.2):0.2,(g:0.2,h:0.2):0.2):0.2):0.2,(d:0.2,c:0.2):0.2):0.2,(b:0.2,a:0.2):0.2);",
            ),
            (
                "(((d:0.3,e:0.3):0.3,((f:0.3,g:0.3):0.3,(h:0.3,(i:0.3,j:0.3):0.3):0.3):0.3):0.3,(a:0.3,(b:0.3,c:0.3):0.3):0.3);",
                "((((g:0.3,f:0.3):0.3,((i:0.3,j:0.3):0.3,h:0.3):0.3):0.3,(d:0.3,e:0.3):0.3):0.3,((b:0.3,c:0.3):0.3,a:0.3):0.3);",
            ),
        ];

        for (newick, rot_newick) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let rota = Tree::from_newick(rot_newick).unwrap();

            let ps_orig: HashSet<_> = HashSet::from_iter(tree.get_partitions().unwrap());
            let ps_rota: HashSet<_> = HashSet::from_iter(rota.get_partitions().unwrap());

            assert_eq!(ps_orig, ps_rota);
        }
    }

    #[test]
    fn self_rf() {
        let test_cases = vec![
            (
                "(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);",
                "(Tip0,((Tip2,(Tip3,(Tip4,(Tip5,(Tip6,((Tip8,Tip9),Tip7)))))),Tip1));",
            ),
            (
                "(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);",
                "(((c:0.1,d:0.1):0.1,((g:0.1,h:0.1):0.1,(f:0.1,e:0.1):0.1):0.1):0.1,((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1);",
            ),
            (
                "((a:0.2,b:0.2):0.2,((c:0.2,d:0.2):0.2,((e:0.2,f:0.2):0.2,((g:0.2,h:0.2):0.2,(i:0.2,j:0.2):0.2):0.2):0.2):0.2);",
                "((((e:0.2,f:0.2):0.2,((i:0.2,j:0.2):0.2,(g:0.2,h:0.2):0.2):0.2):0.2,(d:0.2,c:0.2):0.2):0.2,(b:0.2,a:0.2):0.2);",
            ),
            (
                "(((d:0.3,e:0.3):0.3,((f:0.3,g:0.3):0.3,(h:0.3,(i:0.3,j:0.3):0.3):0.3):0.3):0.3,(a:0.3,(b:0.3,c:0.3):0.3):0.3);",
                "((((g:0.3,f:0.3):0.3,((i:0.3,j:0.3):0.3,h:0.3):0.3):0.3,(d:0.3,e:0.3):0.3):0.3,((b:0.3,c:0.3):0.3,a:0.3):0.3);",
            ),
        ];

        for (newick, rot_newick) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let rota = Tree::from_newick(rot_newick).unwrap();

            let tree_index = tree.taxon_index.clone();
            let rota_index = rota.taxon_index.clone();
            assert_eq!(
                tree.unrooted_robinson_foulds(&rota).unwrap(),
                0,
                "Ref{:#?}\nRot:{:#?}",
                tree.taxon_index(),
                rota.taxon_index()
            );
            assert_eq!(tree.taxon_index, tree_index);
            assert_eq!(rota.taxon_index, rota_index);
        }
    }

    #[test]
    // Robinson foulds distances according to
    // https://evolution.genetics.washington.edu/phylip/doc/treedist.html
    fn robinson_foulds_treedist() {
        let trees = [
            "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        ];
        let rfs = [
            vec![0, 4, 2, 10, 10, 10, 10, 10, 10, 10, 2, 10],
            vec![4, 0, 2, 10, 8, 10, 8, 10, 8, 10, 2, 10],
            vec![2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
            vec![10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
            vec![10, 8, 10, 2, 0, 4, 2, 4, 2, 2, 10, 4],
            vec![10, 10, 10, 2, 4, 0, 2, 2, 4, 2, 10, 2],
            vec![10, 8, 10, 4, 2, 2, 0, 4, 2, 4, 10, 4],
            vec![10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
            vec![10, 8, 10, 4, 2, 4, 2, 2, 0, 4, 10, 2],
            vec![10, 10, 10, 0, 2, 2, 4, 2, 4, 0, 10, 2],
            vec![2, 2, 0, 10, 10, 10, 10, 10, 10, 10, 0, 10],
            vec![10, 10, 10, 2, 4, 2, 4, 0, 2, 2, 10, 0],
        ];

        for indices in (0..trees.len()).combinations(2) {
            let (i0, i1) = (indices[0], indices[1]);

            let t0 = Tree::from_newick(trees[i0]).unwrap();
            let t1 = Tree::from_newick(trees[i1]).unwrap();

            assert_eq!(t0.unrooted_robinson_foulds(&t1).unwrap(), rfs[i0][i1])
        }
    }

    #[test]
    // Robinson foulds distances according to
    // https://evolution.genetics.washington.edu/phylip/doc/treedist.html
    fn weighted_robinson_foulds_treedist() {
        let trees = [
            "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        ];
        let rfs = [
            [
                0.,
                0.4,
                0.2,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.2,
                0.9999999999999999,
            ],
            [
                0.4,
                0.,
                0.2,
                0.9999999999999999,
                0.7999999999999999,
                0.9999999999999999,
                0.7999999999999999,
                0.9999999999999999,
                0.7999999999999999,
                0.9999999999999999,
                0.2,
                0.9999999999999999,
            ],
            [
                0.2,
                0.2,
                0.,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.,
                0.9999999999999999,
            ],
            [
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.,
                0.2,
                0.2,
                0.4,
                0.2,
                0.4,
                0.,
                0.9999999999999999,
                0.2,
            ],
            [
                0.9999999999999999,
                0.7999999999999999,
                0.9999999999999999,
                0.2,
                0.,
                0.4,
                0.2,
                0.4,
                0.2,
                0.2,
                0.9999999999999999,
                0.4,
            ],
            [
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.2,
                0.4,
                0.,
                0.2,
                0.2,
                0.4,
                0.2,
                0.9999999999999999,
                0.2,
            ],
            [
                0.9999999999999999,
                0.7999999999999999,
                0.9999999999999999,
                0.4,
                0.2,
                0.2,
                0.,
                0.4,
                0.2,
                0.4,
                0.9999999999999999,
                0.4,
            ],
            [
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.2,
                0.4,
                0.2,
                0.4,
                0.,
                0.2,
                0.2,
                0.9999999999999999,
                0.,
            ],
            [
                0.9999999999999999,
                0.7999999999999999,
                0.9999999999999999,
                0.4,
                0.2,
                0.4,
                0.2,
                0.2,
                0.,
                0.4,
                0.9999999999999999,
                0.2,
            ],
            [
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.,
                0.2,
                0.2,
                0.4,
                0.2,
                0.4,
                0.,
                0.9999999999999999,
                0.2,
            ],
            [
                0.2,
                0.2,
                0.,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.,
                0.9999999999999999,
            ],
            [
                0.9999999999999999,
                0.9999999999999999,
                0.9999999999999999,
                0.2,
                0.4,
                0.2,
                0.4,
                0.,
                0.2,
                0.2,
                0.9999999999999999,
                0.,
            ],
        ];

        for indices in (0..trees.len()).combinations(2) {
            let (i0, i1) = (indices[0], indices[1]);
            let t0 = Tree::from_newick(trees[i0]).unwrap();
            let t1 = Tree::from_newick(trees[i1]).unwrap();

            assert!((t0.weighted_robinson_foulds(&t1).unwrap() - rfs[i0][i1]).abs() <= f64::EPSILON)
        }
    }

    #[test]
    // Branch score distances according to
    // https://evolution.genetics.washington.edu/phylip/doc/treedist.html
    fn kuhner_felsenstein_treedist() {
        let trees = [
            "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        ];
        let rfs = [
            [
                0.,
                0.2,
                0.14142135623730953,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.14142135623730953,
                0.316227766016838,
            ],
            [
                0.2,
                0.,
                0.14142135623730953,
                0.316227766016838,
                0.28284271247461906,
                0.316227766016838,
                0.28284271247461906,
                0.316227766016838,
                0.28284271247461906,
                0.316227766016838,
                0.14142135623730953,
                0.316227766016838,
            ],
            [
                0.14142135623730953,
                0.14142135623730953,
                0.,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.,
                0.316227766016838,
            ],
            [
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.,
                0.14142135623730953,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.2,
                0.,
                0.316227766016838,
                0.14142135623730953,
            ],
            [
                0.316227766016838,
                0.28284271247461906,
                0.316227766016838,
                0.14142135623730953,
                0.,
                0.2,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.14142135623730953,
                0.316227766016838,
                0.2,
            ],
            [
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.14142135623730953,
                0.2,
                0.,
                0.14142135623730953,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.316227766016838,
                0.14142135623730953,
            ],
            [
                0.316227766016838,
                0.28284271247461906,
                0.316227766016838,
                0.2,
                0.14142135623730953,
                0.14142135623730953,
                0.,
                0.2,
                0.14142135623730953,
                0.2,
                0.316227766016838,
                0.2,
            ],
            [
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.2,
                0.,
                0.14142135623730953,
                0.14142135623730953,
                0.316227766016838,
                0.,
            ],
            [
                0.316227766016838,
                0.28284271247461906,
                0.316227766016838,
                0.2,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.14142135623730953,
                0.,
                0.2,
                0.316227766016838,
                0.14142135623730953,
            ],
            [
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.,
                0.14142135623730953,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.2,
                0.,
                0.316227766016838,
                0.14142135623730953,
            ],
            [
                0.14142135623730953,
                0.14142135623730953,
                0.,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.,
                0.316227766016838,
            ],
            [
                0.316227766016838,
                0.316227766016838,
                0.316227766016838,
                0.14142135623730953,
                0.2,
                0.14142135623730953,
                0.2,
                0.,
                0.14142135623730953,
                0.14142135623730953,
                0.316227766016838,
                0.,
            ],
        ];

        for indices in (0..trees.len()).combinations(2) {
            let (i0, i1) = (indices[0], indices[1]);
            let t0 = Tree::from_newick(trees[i0]).unwrap();
            let t1 = Tree::from_newick(trees[i1]).unwrap();

            println!(
                "[{i0}, {i1}] c:{:?} ==? t:{}",
                t0.kuhner_felsenstein(&t1).unwrap(),
                rfs[i0][i1]
            );

            assert_eq!(t0.kuhner_felsenstein(&t1).unwrap(), rfs[i0][i1])
        }
    }

    #[test]
    fn rooted_vs_unrooted_partitions() {
        let rooted = Tree::from_newick("((Tip_3,Tip_4),(Tip_0,(Tip_1,Tip_2)));").unwrap();
        let unrooted = Tree::from_newick("(Tip_3,Tip_4,(Tip_0,(Tip_1,Tip_2)));").unwrap();

        let parts_rooted = rooted.get_partitions().unwrap();
        let parts_unrooted = unrooted.get_partitions().unwrap();

        assert_eq!(parts_rooted, parts_unrooted);
    }

    #[test]
    fn new_vs_old() {
        let trees = [
            "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((F:0.1,I:0.1):0.1,(G:0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,(((J:0.1,H:0.1):0.1,D:0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,(G:0.1,((F:0.1,I:0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);",
            "(A:0.1,(B:0.1,(E:0.1,((G:0.1,(F:0.1,I:0.1):0.1):0.1,((J:0.1,(H:0.1,D:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);",
        ];

        for newicks in trees.iter().combinations(2) {
            let t1 = Tree::from_newick(newicks[0]).unwrap();
            let t2 = Tree::from_newick(newicks[1]).unwrap();

            assert_eq!(
                t1.unrooted_robinson_foulds(&t2).unwrap(),
                t1.unrooted_robinson_foulds(&t2).unwrap()
            )
        }
    }

    // the reference distance matrix was computed with ete3
    #[test]
    fn compute_distance_matrix() {
        let tree = Tree::from_newick("((A:0.1,B:0.2)F:0.6,(C:0.3,D:0.4)E:0.5)G;").unwrap();
        let true_dists: HashMap<(String, String), f64> = HashMap::from_iter(vec![
            (("A".into(), "B".into()), 0.30000000000000004),
            (("A".into(), "C".into()), 1.5),
            (("A".into(), "D".into()), 1.6),
            (("B".into(), "C".into()), 1.6),
            (("B".into(), "D".into()), 1.7000000000000002),
            (("C".into(), "D".into()), 0.7),
        ]);

        let matrix = tree.distance_matrix().unwrap();

        for ((n1, n2), dist) in true_dists {
            assert!(
                (dist - matrix.get(&n1, &n2).unwrap()) <= f64::EPSILON,
                "d({n1},{n2}) want:{dist} got:{}",
                matrix.get(&n1, &n2).unwrap()
            )
        }
    }
}

#[cfg(test)]
#[allow(dead_code)]
// These tests are lifted from the ETE3 test suite to ensure that the shared functionnalities
// behave the same between the 2 libraries.
mod tests_ete3_operations {

    use super::*;
    use rand::seq::SliceRandom;

    #[ignore]
    #[test]
    // Do I want to support this ?
    fn export_ordered_features() {
        let tree = Tree::from_newick("((A,B),C);").unwrap();
        let expected_newick = "((A:1[&&NHX:dist=1.0:name=A:support=1.0],B:1[&&NHX:0=0:1=1:2=2:3=3:4=4:5=5:6=6:7=7:8=8:9=9:a=a:b=b:c=c:d=d:dist=1.0:e=e:f=f:g=g:h=h:i=i:j=j:k=k:l=l:m=m:n=n:name=B:o=o:p=p:q=q:r=r:s=s:support=1.0:t=t:u=u:v=v:w=w])1:1[&&NHX:dist=1.0:name=:support=1.0],C:1[&&NHX:dist=1.0:name=C:support=1.0]);";

        // Add features to tree here
        let mut features: Vec<_> = String::from("abcdefghijklmnopqrstuvw0123456789")
            .chars()
            .collect();
        features.shuffle(&mut rand::thread_rng());
        let _target = tree.get_by_name("B").unwrap();
        for feature in features {
            let _val = feature.to_string().as_str();
            // target.add_feature(val, val)
        }

        assert_eq!(
            expected_newick,
            tree.to_newick().expect("Could not write newick")
        );
    }

    #[test]
    fn resolve_polytomies() {
        let mut tree = Tree::from_newick("((a,a,a),(b,b,b,(c,c,c)));").unwrap();

        assert!(tree.nodes.iter().any(|n| n.children.len() > 2));
        tree.resolve().unwrap();
        assert!(!tree.nodes.iter().any(|n| n.children.len() > 2));
    }

    #[test]
    fn edge_distances() {
        // Modified the ete3 test tree since this library does not handle NHX comments
        let tree = Tree::from_newick("(((A:0.1, B:0.01):0.001, C:0.0001)I:1.0[&&NHX:name=I], (D:0.00001)J:0.000001[&&NHX:name=J])root:2.0[&&NHX:name=root];").unwrap();

        let a = tree.get_by_name("A").unwrap();
        let b = tree.get_by_name("B").unwrap();
        let d = tree.get_by_name("D").unwrap();
        let i = tree.get_by_name("I").unwrap();
        let root = tree.get_by_name("root").unwrap();

        assert_eq!(tree.get_common_ancestor(a.id, i.id).unwrap(), i.id);
        assert_eq!(tree.get_common_ancestor(a.id, d.id).unwrap(), root.id);

        assert_eq!(tree.get_distance(a.id, i.id).unwrap().0.unwrap(), 0.101);
        assert_eq!(tree.get_distance(a.id, b.id).unwrap().0.unwrap(), 0.11);
        assert_eq!(tree.get_distance(a.id, a.id).unwrap().0.unwrap(), 0.0);
        assert_eq!(tree.get_distance(i.id, i.id).unwrap().0.unwrap(), 0.0);
        assert_eq!(
            tree.get_distance(a.id, root.id).unwrap(),
            tree.get_distance(root.id, a.id).unwrap()
        );
    }

    #[test]
    fn topological_distances() {
        let tree =
            Tree::from_newick("(((A:0.5, B:1.0):1.0, C:5.0):1, (D:10.0, F:1.0):2.0)root:20;")
                .unwrap();
        let a = tree.get_by_name("A").unwrap();
        let c = tree.get_by_name("C").unwrap();
        let d = tree.get_by_name("D").unwrap();
        let root = tree.get_by_name("root").unwrap();

        assert_eq!(tree.get_distance(root.id, a.id).unwrap().1, 3);
        assert_eq!(tree.get_distance(root.id, c.id).unwrap().1, 2);
        assert_eq!(tree.get_distance(root.id, root.id).unwrap().1, 0);
        assert_eq!(tree.get_distance(d.id, d.id).unwrap().1, 0);
    }

    #[test]
    fn traversals() {
        let tree = Tree::from_newick("((3,4)2,(6,7)5)1;").unwrap();
        let root = tree.root().unwrap();

        let postorder = "3426751";
        let preorder = "1234567";
        let levelorder = "1253467";

        fn get_str(iter: &[NodeId], tree: &Tree) -> String {
            iter.iter()
                .map(|&id| tree.get(id).unwrap().name.clone().unwrap())
                .collect()
        }

        assert_eq!(get_str(&tree.postorder(root).unwrap(), &tree), postorder);
        assert_eq!(get_str(&tree.preorder(root).unwrap(), &tree), preorder);
        assert_eq!(get_str(&tree.levelorder(root).unwrap(), &tree), levelorder);
    }

    #[ignore]
    #[test]
    fn compare() {
        let cases = [
            (
                28,
                true,
                "(((z,y),(x,(w,v))),(u,t),((s,r),((q,(p,o)),((n,(m,(l,(k,j)))),(i,(h,g))))));",
                "(((k,(j,(i,(h,g)))),z),(y,x),((w,v),((u,(t,(s,(r,q)))),(p,(o,(n,(m,l)))))));",
            ),
            (
                28,
                false,
                "(((t,s),((r,(q,p)),(o,n))),(((m,(l,(k,j))),(i,(h,g))),(z,(y,(x,(w,(v,u)))))));",
                "((((k,(j,i)),((h,g),z)),((y,(x,w)),((v,(u,t)),(s,(r,(q,p)))))),((o,n),(m,l)));",
            ),
            (
                18,
                true,
                "(((v,(u,(t,s))),((r,(q,(p,o))),((n,m),(l,k)))),(j,(i,(h,g))),(z,(y,(x,w))));",
                "(((z,(y,(x,w))),(v,(u,(t,s)))),((r,(q,p)),(o,(n,m))),((l,(k,(j,i))),(h,g)));",
            ),
            (
                26,
                true,
                "(((l,k),(j,i)),((h,g),(z,(y,(x,w)))),((v,(u,(t,(s,(r,q))))),((p,o),(n,m))));",
                "(((p,o),((n,(m,l)),(k,j))),((i,(h,g)),(z,y)),((x,(w,v)),((u,(t,s)),(r,q))));",
            ),
            (
                24,
                true,
                "(((o,(n,m)),(l,(k,(j,(i,(h,g)))))),(z,(y,x)),((w,v),((u,(t,(s,r))),(q,p))));",
                "(((t,(s,(r,(q,(p,o))))),(n,m)),((l,k),(j,(i,(h,g)))),((z,y),((x,w),(v,u))));",
            ),
            (
                24,
                true,
                "(((y,(x,(w,v))),(u,t)),((s,(r,(q,(p,o)))),(n,m)),((l,k),((j,(i,(h,g))),z)));",
                "(((z,(y,(x,w))),(v,(u,t))),(s,(r,(q,(p,(o,(n,(m,(l,k)))))))),(j,(i,(h,g))));",
            ),
            (
                28,
                false,
                "(((p,(o,(n,(m,l)))),((k,(j,i)),(h,g))),((z,y),((x,(w,(v,u))),(t,(s,(r,q))))));",
                "((((t,(s,r)),(q,p)),((o,n),(m,(l,(k,(j,i)))))),(((h,g),(z,(y,(x,w)))),(v,u)));",
            ),
            (
                28,
                true,
                "((((i,(h,g)),z),(y,x)),((w,v),((u,(t,(s,r))),(q,p))),((o,n),(m,(l,(k,j)))));",
                "((((h,g),z),(y,x)),(w,(v,u)),((t,s),((r,(q,p)),((o,(n,m)),(l,(k,(j,i)))))));",
            ),
            (
                28,
                true,
                "(((x,(w,(v,(u,(t,(s,(r,(q,(p,o))))))))),((n,(m,l)),(k,(j,i)))),(h,g),(z,y));",
                "(((u,t),(s,r)),((q,p),(o,(n,m))),(((l,(k,(j,i))),((h,g),(z,(y,x)))),(w,v)));",
            ),
            (
                22,
                false,
                "(((x,(w,(v,u))),((t,(s,r)),(q,p))),((o,(n,(m,l))),((k,j),((i,(h,g)),(z,y)))));",
                "(((z,(y,(x,(w,(v,u))))),(t,(s,r))),((q,(p,(o,(n,m)))),((l,k),(j,(i,(h,g))))));",
            ),
            (
                26,
                true,
                "((z,(y,(x,w))),(v,(u,(t,s))),((r,(q,(p,(o,(n,m))))),((l,k),(j,(i,(h,g))))));",
                "(((v,(u,t)),((s,r),((q,(p,o)),(n,(m,l))))),((k,j),((i,(h,g)),z)),(y,(x,w)));",
            ),
            (
                34,
                false,
                "((((i,(h,g)),(z,(y,x))),(w,v)),((u,t),((s,r),((q,(p,(o,n))),(m,(l,(k,j)))))));",
                "(((p,(o,(n,(m,(l,k))))),((j,i),(h,g))),(z,(y,(x,(w,(v,(u,(t,(s,(r,q))))))))));",
            ),
            (
                30,
                false,
                "(((i,(h,g)),(z,y)),((x,w),((v,(u,(t,(s,(r,q))))),(p,(o,(n,(m,(l,(k,j)))))))));",
                "((((l,k),(j,(i,(h,g)))),(z,(y,(x,w)))),((v,u),((t,s),((r,(q,p)),(o,(n,m))))));",
            ),
            (
                26,
                false,
                "(((v,(u,t)),((s,(r,q)),((p,o),((n,m),((l,k),(j,i)))))),((h,g),(z,(y,(x,w)))));",
                "(((y,(x,(w,v))),(u,(t,s))),(((r,q),((p,o),(n,(m,(l,k))))),((j,i),((h,g),z))));",
            ),
            (
                20,
                false,
                "(((u,(t,s)),(r,q)),(((p,o),((n,m),((l,k),((j,i),((h,g),z))))),(y,(x,(w,v)))));",
                "((((u,t),(s,r)),(((q,p),(o,(n,m))),(((l,k),(j,i)),((h,g),z)))),((y,x),(w,v)));",
            ),
            (
                20,
                true,
                "(((y,x),(w,v)),((u,(t,s)),((r,q),(p,(o,(n,(m,(l,k))))))),((j,(i,(h,g))),z));",
                "(((r,q),((p,o),(n,(m,(l,(k,j)))))),((i,(h,g)),(z,(y,(x,(w,v))))),(u,(t,s)));",
            ),
            (
                24,
                true,
                "((((k,(j,i)),(h,g)),((z,(y,(x,w))),((v,(u,t)),(s,r)))),(q,(p,(o,n))),(m,l));",
                "((((s,r),((q,p),(o,(n,m)))),((l,k),((j,i),((h,g),z)))),(y,x),(w,(v,(u,t))));",
            ),
            (
                18,
                true,
                "((w,(v,(u,(t,s)))),(r,q),((p,(o,n)),((m,(l,k)),((j,(i,(h,g))),(z,(y,x))))));",
                "(((y,x),((w,v),(u,(t,s)))),((r,(q,(p,(o,n)))),(m,l)),((k,j),((i,(h,g)),z)));",
            ),
            (
                26,
                true,
                "(((j,(i,(h,g))),(z,(y,(x,(w,(v,(u,t))))))),(s,r),((q,p),((o,(n,m)),(l,k))));",
                "(((s,(r,(q,(p,(o,(n,(m,l))))))),(k,j)),((i,(h,g)),(z,y)),((x,(w,v)),(u,t)));",
            ),
            (
                30,
                true,
                "((((r,(q,(p,(o,n)))),((m,l),(k,(j,i)))),((h,g),z)),(y,(x,(w,v))),(u,(t,s)));",
                "(((u,t),(s,r)),((q,p),(o,(n,(m,(l,(k,j)))))),(((i,(h,g)),(z,(y,x))),(w,v)));",
            ),
            (
                30,
                false,
                "((((m,(l,k)),(j,i)),(((h,g),(z,y)),(x,w))),((v,u),(t,(s,(r,(q,(p,(o,n))))))));",
                "(((u,t),((s,(r,q)),(p,(o,(n,(m,(l,k))))))),((j,(i,(h,g))),(z,(y,(x,(w,v))))));",
            ),
            (
                22,
                false,
                "(((k,(j,i)),(h,g)),((z,(y,x)),((w,(v,(u,(t,(s,r))))),((q,(p,(o,n))),(m,l)))));",
                "(((w,(v,u)),((t,(s,r)),((q,p),((o,(n,(m,l))),((k,(j,i)),((h,g),z)))))),(y,x));",
            ),
            (
                26,
                false,
                "(((x,(w,(v,(u,(t,s))))),(r,q)),((p,(o,(n,(m,l)))),((k,j),((i,(h,g)),(z,y)))));",
                "(((o,(n,m)),(l,(k,j))),(((i,(h,g)),(z,y)),((x,w),((v,u),((t,(s,r)),(q,p))))));",
            ),
            (
                28,
                true,
                "(((x,(w,v)),(u,(t,s))),((r,(q,(p,(o,(n,m))))),(l,(k,(j,(i,(h,g)))))),(z,y));",
                "((((i,(h,g)),(z,(y,x))),((w,v),((u,t),(s,(r,(q,p)))))),(o,n),((m,l),(k,j)));",
            ),
            (
                20,
                false,
                "((((m,l),(k,(j,(i,(h,g))))),(z,y)),((x,(w,(v,(u,(t,s))))),(r,(q,(p,(o,n))))));",
                "((((m,l),((k,(j,i)),(h,g))),(z,(y,(x,(w,v))))),((u,t),(s,(r,(q,(p,(o,n)))))));",
            ),
            (
                26,
                true,
                "(((o,(n,(m,(l,k)))),(j,i)),((h,g),(z,y)),((x,(w,(v,(u,(t,s))))),(r,(q,p))));",
                "((((t,(s,(r,(q,(p,(o,n)))))),(m,(l,k))),((j,i),(h,g))),(z,(y,x)),(w,(v,u)));",
            ),
            (
                22,
                false,
                "((((p,o),((n,m),((l,k),(j,i)))),((h,g),(z,y))),((x,(w,(v,u))),((t,s),(r,q))));",
                "((((v,(u,(t,s))),(r,q)),((p,o),((n,m),(l,k)))),(((j,i),(h,g)),(z,(y,(x,w)))));",
            ),
            (
                28,
                false,
                "((((r,(q,(p,(o,n)))),(m,(l,k))),(((j,i),(h,g)),((z,y),(x,w)))),((v,u),(t,s)));",
                "((((k,j),((i,(h,g)),(z,y))),(x,w)),(((v,(u,t)),(s,r)),((q,p),((o,n),(m,l)))));",
            ),
            (
                20,
                true,
                "((((q,(p,o)),(n,m)),((l,k),((j,i),(h,g)))),(z,(y,x)),((w,v),(u,(t,(s,r)))));",
                "((((l,(k,(j,i))),(h,g)),((z,y),(x,(w,v)))),(u,t),((s,(r,(q,(p,o)))),(n,m)));",
            ),
            (
                28,
                false,
                "(((t,(s,r)),(q,(p,o))),(((n,(m,(l,k))),(j,(i,(h,g)))),((z,y),(x,(w,(v,u))))));",
                "(((w,(v,u)),(t,s)),(((r,(q,p)),(o,n)),(((m,l),((k,j),((i,(h,g)),z))),(y,x))));",
            ),
            (
                24,
                true,
                "((((h,g),(z,y)),((x,(w,(v,u))),(t,(s,(r,q))))),(p,o),((n,m),((l,k),(j,i))));",
                "(((t,s),((r,(q,p)),((o,(n,(m,l))),((k,j),(i,(h,g)))))),(z,y),(x,(w,(v,u))));",
            ),
            (
                20,
                true,
                "(((p,o),(n,(m,(l,(k,(j,i)))))),((h,g),z),((y,(x,w)),((v,u),(t,(s,(r,q))))));",
                "(((y,(x,w)),(v,(u,t))),((s,r),(q,p)),((o,(n,m)),((l,(k,(j,i))),((h,g),z))));",
            ),
            (
                32,
                true,
                "((((s,(r,q)),((p,(o,n)),(m,(l,k)))),((j,(i,(h,g))),(z,y))),(x,w),(v,(u,t)));",
                "(((u,(t,(s,r))),((q,(p,o)),((n,(m,l)),(k,(j,i))))),((h,g),(z,(y,x))),(w,v));",
            ),
            (
                26,
                true,
                "(((z,(y,x)),(w,(v,(u,t)))),(s,(r,(q,(p,(o,n))))),((m,l),(k,(j,(i,(h,g))))));",
                "(((u,t),((s,r),((q,p),((o,n),((m,(l,k)),((j,i),((h,g),z))))))),(y,x),(w,v));",
            ),
            (
                10,
                true,
                "(((p,o),((n,m),((l,(k,(j,i))),((h,g),(z,y))))),(x,(w,(v,u))),((t,s),(r,q)));",
                "((((n,m),((l,(k,(j,i))),((h,g),(z,y)))),(x,w)),(v,(u,(t,(s,(r,q))))),(p,o));",
            ),
            (
                30,
                true,
                "((((h,g),z),((y,x),((w,v),(u,t)))),(s,r),((q,p),((o,n),((m,l),(k,(j,i))))));",
                "((((v,(u,(t,(s,r)))),(q,(p,o))),((n,m),((l,k),(j,(i,(h,g)))))),(z,y),(x,w));",
            ),
            (
                30,
                false,
                "(((q,(p,o)),((n,m),((l,(k,(j,(i,(h,g))))),(z,y)))),((x,(w,v)),(u,(t,(s,r)))));",
                "((((t,s),((r,q),((p,o),(n,m)))),((l,k),(j,i))),(((h,g),z),((y,(x,w)),(v,u))));",
            ),
            (
                // This is the only test case that does not work. I'm still not sure why RF should
                // 24...
                24,
                false,
                "(((p,o),(n,m)),(((l,(k,(j,i))),(h,g)),((z,y),((x,w),((v,u),(t,(s,(r,q))))))));",
                "((x,(w,v)),((u,(t,(s,(r,q)))),((p,(o,(n,(m,(l,(k,(j,(i,(h,g))))))))),(z,y))));",
            ),
            (
                28,
                false,
                "(((z,y),((x,w),((v,u),(t,s)))),((r,(q,(p,(o,(n,m))))),((l,k),((j,i),(h,g)))));",
                "((((s,(r,q)),((p,o),((n,(m,l)),(k,(j,(i,(h,g))))))),(z,y)),((x,w),(v,(u,t))));",
            ),
            (
                24,
                false,
                "((((o,n),((m,l),((k,(j,i)),(h,g)))),(z,(y,x))),((w,(v,(u,(t,(s,r))))),(q,p)));",
                "(((q,(p,(o,(n,m)))),((l,(k,j)),(i,(h,g)))),(z,(y,(x,(w,(v,(u,(t,(s,r)))))))));",
            ),
            (
                22,
                true,
                "(((p,(o,(n,m))),((l,k),((j,i),((h,g),(z,y))))),(x,w),((v,u),((t,s),(r,q))));",
                "(((u,(t,(s,(r,(q,(p,(o,(n,m)))))))),((l,k),((j,i),((h,g),(z,(y,x)))))),w,v);",
            ),
            (
                28,
                false,
                "((((r,q),((p,o),(n,(m,l)))),((k,(j,i)),(h,g))),((z,y),((x,(w,v)),(u,(t,s)))));",
                "(((h,g),z),((y,x),((w,v),((u,t),((s,(r,(q,(p,(o,(n,m)))))),(l,(k,(j,i))))))));",
            ),
            (
                30,
                true,
                "((((h,g),z),((y,(x,(w,(v,u)))),((t,s),((r,(q,(p,o))),(n,m))))),(l,k),(j,i));",
                "((((o,n),((m,(l,(k,j))),((i,(h,g)),z))),(y,(x,(w,v)))),(u,(t,s)),(r,(q,p)));",
            ),
            (
                30,
                true,
                "(((v,u),(t,(s,(r,(q,p))))),((o,(n,m)),((l,(k,j)),((i,(h,g)),z))),(y,(x,w)));",
                "((((m,(l,k)),((j,i),(h,g))),(z,y)),(x,w),((v,(u,(t,(s,(r,q))))),(p,(o,n))));",
            ),
            (
                26,
                true,
                "(((q,p),((o,(n,(m,l))),(k,(j,i)))),((h,g),z),((y,x),((w,(v,(u,t))),(s,r))));",
                "((((j,(i,(h,g))),(z,(y,x))),((w,v),(u,t))),(s,(r,q)),((p,o),(n,(m,(l,k)))));",
            ),
            (
                20,
                false,
                "((((o,(n,m)),((l,k),((j,i),((h,g),z)))),(y,x)),(((w,v),(u,t)),((s,r),(q,p))));",
                "((((j,i),((h,g),z)),((y,x),(w,(v,(u,(t,(s,r))))))),((q,p),((o,n),(m,(l,k)))));",
            ),
            (
                30,
                false,
                "(((x,w),(v,(u,(t,(s,(r,(q,(p,(o,(n,m)))))))))),((l,k),((j,(i,(h,g))),(z,y))));",
                "(((m,l),((k,(j,(i,(h,g)))),z)),((y,(x,(w,(v,(u,t))))),((s,r),((q,p),(o,n)))));",
            ),
            (
                32,
                true,
                "((((y,x),(w,v)),((u,(t,(s,r))),(q,(p,o)))),((n,m),(l,(k,j))),((i,(h,g)),z));",
                "(((m,l),(k,(j,i))),((h,g),z),((y,(x,w)),((v,u),((t,s),(r,(q,(p,(o,n))))))));",
            ),
            (
                28,
                true,
                "(((v,u),((t,(s,(r,(q,p)))),((o,n),((m,l),(k,(j,(i,(h,g)))))))),(z,y),(x,w));",
                "((((n,m),((l,k),((j,i),((h,g),(z,(y,(x,(w,(v,u))))))))),(t,s)),(r,q),(p,o));",
            ),
            (
                32,
                false,
                "(((r,(q,p)),(o,n)),(((m,(l,k)),(j,i)),(((h,g),(z,y)),((x,w),((v,u),(t,s))))));",
                "(((y,x),((w,v),(u,(t,(s,r))))),(((q,(p,(o,n))),(m,l)),((k,(j,(i,(h,g)))),z)));",
            ),
            (
                20,
                true,
                "(((w,v),((u,(t,(s,r))),((q,p),((o,(n,(m,l))),((k,j),((i,(h,g)),z)))))),y,x);",
                "(((w,v),((u,t),(s,(r,q)))),((p,o),((n,(m,l)),(k,j))),((i,(h,g)),(z,(y,x))));",
            ),
            (
                24,
                false,
                "(((x,(w,v)),((u,(t,s)),(r,q))),(((p,o),((n,(m,l)),(k,j))),((i,(h,g)),(z,y))));",
                "((((i,(h,g)),z),((y,x),(w,v))),((u,(t,s)),((r,(q,(p,(o,(n,m))))),(l,(k,j)))));",
            ),
            (
                22,
                false,
                "((((k,(j,(i,(h,g)))),(z,(y,x))),((w,v),(u,t))),((s,(r,(q,(p,o)))),(n,(m,l))));",
                "(((w,v),(u,(t,(s,(r,(q,(p,o))))))),(((n,m),((l,(k,(j,i))),((h,g),z))),(y,x)));",
            ),
            (
                28,
                true,
                "(((x,w),((v,u),((t,s),(r,(q,p))))),((o,n),(m,l)),((k,(j,i)),((h,g),(z,y))));",
                "((((p,o),(n,m)),((l,(k,(j,i))),((h,g),z))),(y,(x,(w,v))),((u,t),(s,(r,q))));",
            ),
            (
                30,
                false,
                "(((q,p),((o,(n,(m,l))),((k,(j,(i,(h,g)))),z))),((y,x),((w,(v,u)),(t,(s,r)))));",
                "((((m,(l,k)),((j,(i,(h,g))),z)),(y,(x,w))),((v,(u,(t,(s,(r,q))))),(p,(o,n))));",
            ),
            (
                30,
                false,
                "(((y,x),((w,(v,(u,(t,(s,r))))),(q,p))),((o,(n,(m,(l,(k,(j,i)))))),((h,g),z)));",
                "((((t,(s,(r,q))),((p,(o,(n,(m,l)))),((k,(j,i)),(h,g)))),(z,y)),((x,w),(v,u)));",
            ),
            (
                20,
                false,
                "(((u,(t,s)),(r,(q,(p,(o,(n,(m,(l,(k,j))))))))),(((i,(h,g)),z),(y,(x,(w,v)))));",
                "(((o,n),(m,(l,(k,j)))),(((i,(h,g)),(z,y)),((x,(w,v)),((u,(t,(s,r))),(q,p)))));",
            ),
            (
                26,
                false,
                "(((t,s),((r,(q,(p,(o,n)))),(m,(l,k)))),(((j,i),((h,g),z)),((y,(x,w)),(v,u))));",
                "(((r,(q,(p,o))),((n,(m,(l,k))),((j,i),(h,g)))),((z,(y,(x,(w,v)))),(u,(t,s))));",
            ),
            (
                28,
                true,
                "((((r,q),((p,(o,(n,(m,l)))),((k,(j,i)),(h,g)))),(z,(y,(x,w)))),(v,u),(t,s));",
                "(((x,(w,(v,(u,(t,s))))),(r,(q,(p,o)))),(n,m),((l,k),((j,(i,(h,g))),(z,y))));",
            ),
            (
                28,
                false,
                "(((t,s),((r,(q,p)),((o,n),(m,(l,(k,(j,i))))))),(((h,g),(z,y)),(x,(w,(v,u)))));",
                "((((h,g),(z,(y,(x,(w,v))))),(u,(t,(s,r)))),((q,(p,(o,(n,m)))),(l,(k,(j,i)))));",
            ),
            (
                26,
                true,
                "((((q,(p,o)),((n,m),((l,(k,(j,i))),(h,g)))),(z,(y,x))),(w,v),(u,(t,(s,r))));",
                "(((y,x),(w,(v,u))),((t,(s,r)),((q,p),(o,n))),((m,(l,k)),((j,(i,(h,g))),z)));",
            ),
            (
                28,
                false,
                "((((q,(p,(o,n))),((m,(l,k)),((j,(i,(h,g))),z))),(y,x)),((w,(v,(u,t))),(s,r)));",
                "(((z,(y,x)),(w,v)),(((u,t),((s,(r,(q,p))),((o,n),(m,l)))),((k,(j,i)),(h,g))));",
            ),
            (
                22,
                true,
                "(((x,w),((v,(u,(t,s))),(r,q))),((p,(o,n)),((m,(l,k)),(j,(i,(h,g))))),(z,y));",
                "((((j,(i,(h,g))),(z,(y,x))),(w,(v,u))),((t,s),((r,q),(p,o))),((n,m),(l,k)));",
            ),
            (
                26,
                false,
                "((((n,(m,l)),(k,j)),(((i,(h,g)),(z,y)),((x,w),((v,u),(t,s))))),((r,q),(p,o)));",
                "(((v,u),(t,s)),(((r,(q,(p,(o,n)))),((m,(l,k)),(j,i))),((h,g),(z,(y,(x,w))))));",
            ),
            (
                32,
                false,
                "((((n,(m,(l,(k,j)))),((i,(h,g)),z)),(y,x)),((w,v),((u,(t,(s,r))),(q,(p,o)))));",
                "((((v,u),(t,(s,(r,(q,p))))),((o,(n,(m,(l,k)))),(j,(i,(h,g))))),((z,y),(x,w)));",
            ),
            (
                20,
                false,
                "((((q,(p,(o,n))),(m,l)),((k,(j,(i,(h,g)))),z)),((y,(x,(w,(v,(u,t))))),(s,r)));",
                "(((w,(v,(u,t))),(s,r)),(((q,p),(o,n)),(((m,l),(k,(j,i))),((h,g),(z,(y,x))))));",
            ),
            (
                20,
                true,
                "(((z,(y,(x,w))),(v,u)),((t,(s,r)),(q,(p,o))),((n,(m,l)),((k,(j,i)),(h,g))));",
                "((((q,(p,(o,n))),(m,l)),((k,j),(i,(h,g)))),(z,y),((x,w),((v,u),(t,(s,r)))));",
            ),
            (
                34,
                false,
                "(((w,(v,(u,(t,(s,(r,q)))))),(p,o)),(((n,m),(l,(k,j))),((i,(h,g)),(z,(y,x)))));",
                "(((y,(x,(w,(v,u)))),(t,(s,r))),(((q,(p,(o,(n,(m,(l,k)))))),(j,i)),((h,g),z)));",
            ),
            (
                26,
                false,
                "(((y,x),(w,(v,(u,t)))),(((s,r),((q,(p,o)),(n,(m,l)))),((k,(j,(i,(h,g)))),z)));",
                "(((s,(r,(q,(p,o)))),(n,m)),(((l,k),((j,i),((h,g),(z,(y,(x,w)))))),(v,(u,t))));",
            ),
            (
                30,
                false,
                "(((v,(u,t)),((s,r),((q,p),((o,(n,(m,(l,k)))),(j,i))))),(((h,g),z),(y,(x,w))));",
                "(((y,(x,(w,v))),((u,(t,s)),(r,(q,(p,o))))),((n,(m,l)),((k,(j,i)),((h,g),z))));",
            ),
            (
                26,
                false,
                "(((y,x),(w,v)),(((u,t),((s,(r,(q,p))),(o,n))),((m,(l,k)),((j,i),((h,g),z)))));",
                "((((s,(r,q)),((p,(o,n)),((m,l),(k,(j,i))))),((h,g),z)),((y,(x,w)),(v,(u,t))));",
            ),
            (
                22,
                true,
                "(((w,v),(u,t)),((s,r),((q,p),((o,(n,m)),((l,k),((j,i),(h,g)))))),(z,(y,x)));",
                "(((z,y),(x,(w,(v,u)))),(t,(s,r)),((q,(p,o)),((n,m),((l,(k,(j,i))),(h,g)))));",
            ),
            (
                28,
                false,
                "(((y,x),(w,(v,(u,t)))),(((s,(r,q)),((p,o),(n,(m,(l,k))))),((j,i),((h,g),z))));",
                "((((i,(h,g)),(z,(y,x))),((w,(v,u)),(t,s))),((r,q),((p,o),((n,m),(l,(k,j))))));",
            ),
            (
                26,
                false,
                "(((v,(u,(t,s))),(r,(q,p))),(((o,n),((m,(l,(k,j))),((i,(h,g)),(z,y)))),(x,w)));",
                "(((q,p),((o,n),((m,l),((k,j),((i,(h,g)),z))))),(y,(x,(w,(v,(u,(t,(s,r))))))));",
            ),
            (
                26,
                true,
                "(((t,(s,(r,q))),((p,o),((n,(m,l)),((k,j),((i,(h,g)),z))))),(y,x),(w,(v,u)));",
                "(((z,y),(x,w)),(v,u),((t,(s,r)),((q,(p,(o,(n,(m,l))))),((k,(j,i)),(h,g)))));",
            ),
            (
                30,
                true,
                "(((w,(v,(u,(t,(s,r))))),(q,p)),((o,(n,m)),((l,k),(j,i))),(((h,g),z),(y,x)));",
                "((((p,o),(n,(m,(l,(k,(j,(i,(h,g)))))))),(z,(y,x))),(w,(v,u)),((t,s),(r,q)));",
            ),
            (
                26,
                true,
                "((((i,(h,g)),(z,y)),(x,w)),((v,u),((t,(s,r)),(q,p))),((o,n),(m,(l,(k,j)))));",
                "(((l,k),((j,i),((h,g),(z,y)))),(x,w),((v,u),((t,s),((r,(q,(p,o))),(n,m)))));",
            ),
            (
                26,
                false,
                "(((x,w),((v,(u,(t,s))),((r,(q,p)),((o,(n,(m,(l,k)))),((j,i),(h,g)))))),(z,y));",
                "(((p,(o,(n,m))),(l,k)),(((j,i),(h,g)),((z,y),((x,(w,v)),((u,t),(s,(r,q)))))));",
            ),
            (
                24,
                true,
                "(((x,w),((v,(u,t)),(s,r))),((q,p),(o,(n,(m,(l,k))))),((j,i),((h,g),(z,y))));",
                "(((h,g),(z,y)),(x,(w,(v,u))),((t,(s,r)),(q,(p,(o,(n,(m,(l,(k,(j,i))))))))));",
            ),
            (
                24,
                true,
                "(((y,x),(w,v)),((u,t),((s,r),((q,p),((o,n),(m,(l,k)))))),((j,(i,(h,g))),z));",
                "((((r,(q,p)),(o,(n,(m,(l,(k,(j,(i,(h,g))))))))),(z,y)),(x,(w,v)),(u,(t,s)));",
            ),
            (
                28,
                false,
                "(((y,(x,(w,v))),((u,t),((s,(r,q)),((p,(o,n)),((m,l),(k,(j,i))))))),((h,g),z));",
                "(((v,u),(t,(s,(r,(q,(p,(o,n))))))),(((m,l),((k,j),((i,(h,g)),z))),(y,(x,w))));",
            ),
            (
                26,
                true,
                "((((h,g),z),((y,x),((w,(v,u)),((t,(s,(r,q))),(p,(o,n)))))),(m,(l,k)),(j,i));",
                "((z,y),(x,(w,(v,(u,t)))),((s,r),((q,p),((o,n),((m,(l,k)),(j,(i,(h,g))))))));",
            ),
            (
                24,
                true,
                "(((u,t),(s,r)),((q,p),((o,n),((m,(l,(k,(j,(i,(h,g)))))),z))),(y,(x,(w,v))));",
                "((((j,(i,(h,g))),z),(y,x)),(w,(v,(u,t))),((s,(r,(q,p))),((o,(n,m)),(l,k))));",
            ),
            (
                30,
                true,
                "(((t,(s,r)),((q,p),((o,n),(m,(l,(k,j)))))),((i,(h,g)),z),((y,x),(w,(v,u))));",
                "((((w,(v,(u,t))),(s,(r,q))),((p,(o,(n,m))),(l,k))),((j,i),(h,g)),(z,(y,x)));",
            ),
            (
                30,
                false,
                "((((x,(w,v)),(u,t)),((s,(r,q)),(p,o))),(((n,m),((l,k),((j,i),(h,g)))),(z,y)));",
                "((r,q),((p,(o,n)),((m,(l,(k,(j,i)))),((h,g),(z,(y,(x,(w,(v,(u,(t,s)))))))))));",
            ),
            (
                28,
                true,
                "((((k,j),((i,(h,g)),(z,(y,x)))),(w,v)),(u,t),((s,(r,q)),(p,(o,(n,(m,l))))));",
                "(((z,y),(x,w)),(v,(u,(t,(s,(r,q))))),((p,o),((n,(m,(l,(k,(j,i))))),(h,g))));",
            ),
            (
                18,
                true,
                "(((t,s),((r,(q,(p,o))),(n,m))),((l,(k,j)),((i,(h,g)),(z,y))),((x,w),(v,u)));",
                "((((l,k),(j,i)),(((h,g),(z,y)),(x,w))),((v,u),(t,s)),((r,q),((p,o),(n,m))));",
            ),
            (
                26,
                true,
                "(((h,g),z),(y,(x,w)),((v,(u,(t,s))),((r,(q,p)),((o,(n,(m,l))),(k,(j,i))))));",
                "(((s,r),(q,p)),((o,n),(m,l)),(((k,j),((i,(h,g)),(z,(y,x)))),(w,(v,(u,t)))));",
            ),
            (
                30,
                true,
                "(((x,w),((v,(u,(t,(s,(r,(q,(p,(o,n)))))))),((m,(l,k)),((j,i),(h,g))))),z,y);",
                "((((h,g),z),(y,x)),((w,v),((u,(t,s)),(r,q))),((p,(o,(n,(m,l)))),(k,(j,i))));",
            ),
            (
                30,
                false,
                "(((v,(u,(t,(s,(r,q))))),((p,(o,(n,m))),((l,(k,(j,i))),(h,g)))),((z,y),(x,w)));",
                "(((v,u),((t,(s,(r,(q,(p,o))))),(n,(m,(l,(k,j)))))),((i,(h,g)),(z,(y,(x,w)))));",
            ),
            (
                22,
                true,
                "(((z,y),((x,(w,v)),((u,(t,(s,r))),(q,(p,o))))),(n,m),((l,k),(j,(i,(h,g)))));",
                "(((r,q),(p,(o,(n,m)))),((l,(k,(j,(i,(h,g))))),(z,y)),((x,w),(v,(u,(t,s)))));",
            ),
            (
                30,
                true,
                "(((x,w),((v,(u,(t,(s,r)))),(q,p))),((o,n),(m,l)),((k,j),((i,(h,g)),(z,y))));",
                "((((p,o),((n,(m,(l,k))),((j,i),(h,g)))),((z,y),(x,(w,v)))),(u,t),(s,(r,q)));",
            ),
            (
                32,
                false,
                "(((r,(q,p)),(o,(n,m))),(((l,(k,(j,i))),(h,g)),((z,(y,(x,(w,(v,u))))),(t,s))));",
                "((((j,(i,(h,g))),(z,y)),(x,(w,(v,(u,t))))),(((s,r),(q,(p,o))),((n,m),(l,k))));",
            ),
            (
                30,
                false,
                "((((q,p),((o,(n,(m,(l,k)))),((j,(i,(h,g))),(z,y)))),(x,w)),((v,u),(t,(s,r))));",
                "((((o,(n,m)),((l,(k,(j,i))),((h,g),z))),(y,x)),((w,v),((u,t),((s,r),(q,p)))));",
            ),
            (
                28,
                false,
                "((((s,r),((q,(p,o)),(n,(m,l)))),((k,(j,i)),(h,g))),((z,(y,x)),(w,(v,(u,t)))));",
                "(((m,l),(k,j)),(((i,(h,g)),z),((y,x),((w,(v,(u,(t,(s,r))))),((q,p),(o,n))))));",
            ),
            (
                20,
                true,
                "((((z,y),(x,(w,(v,u)))),((t,s),(r,q))),((p,o),(n,(m,l))),((k,(j,i)),(h,g)));",
                "(((j,i),(h,g)),(z,(y,x)),((w,(v,u)),((t,(s,(r,q))),((p,o),((n,m),(l,k))))));",
            ),
            (
                20,
                false,
                "(((v,u),((t,s),(r,q))),(((p,o),(n,(m,l))),(((k,(j,i)),((h,g),z)),(y,(x,w)))));",
                "((((s,(r,q)),(p,o)),(((n,(m,l)),(k,(j,i))),((h,g),z))),((y,x),((w,v),(u,t))));",
            ),
            (
                28,
                true,
                "((z,y),(x,w),((v,u),((t,(s,(r,q))),((p,(o,(n,m))),(l,(k,(j,(i,(h,g)))))))));",
                "((((r,q),((p,o),((n,m),((l,k),(j,i))))),((h,g),(z,(y,x)))),(w,v),(u,(t,s)));",
            ),
            (
                24,
                false,
                "((((k,(j,(i,(h,g)))),(z,y)),(x,(w,v))),(((u,t),(s,(r,q))),((p,o),(n,(m,l)))));",
                "(((w,v),(u,(t,s))),(((r,(q,(p,o))),((n,m),(l,(k,(j,(i,(h,g))))))),(z,(y,x))));",
            ),
            (
                24,
                true,
                "((((n,m),((l,(k,j)),(i,(h,g)))),(z,y)),(x,(w,v)),((u,(t,(s,(r,q)))),(p,o)));",
                "(((r,q),(p,o)),((n,(m,l)),((k,j),((i,(h,g)),z))),((y,x),(w,(v,(u,(t,s))))));",
            ),
        ];

        let mut failed = vec![];

        for (i, (expected, unrooted, nw1, nw2)) in cases.into_iter().enumerate() {
            let t1 = Tree::from_newick(nw1).unwrap();
            let t2 = Tree::from_newick(nw2).unwrap();

            let rf = t1.unrooted_robinson_foulds(&t2).unwrap();

            // Add information to failure output
            if expected != rf {
                let p1 = t1.get_partitions().unwrap();
                let p2 = t2.get_partitions().unwrap();

                let p1_s: Result<Vec<_>, _> =
                    p1.iter().map(|part| t1.partition_to_leaves(part)).collect();
                let p2_s: Result<Vec<_>, _> =
                    p2.iter().map(|part| t2.partition_to_leaves(part)).collect();

                let mut p1_s = p1_s.unwrap();
                let mut p2_s = p2_s.unwrap();

                p1_s.sort();
                p2_s.sort();

                failed.push((i, rf, expected, nw1, nw2, unrooted, p1_s, p2_s));
            }
        }

        assert!(failed.is_empty(), "Failed cases:\n{failed:#?}")
    }
}

#[cfg(test)]
mod rescale_cache_tests {
    use super::Tree;

    #[test]
    fn rescaling_invalidates_cached_subtree_distances() {
        for target_diameter in [false, true] {
            let mut tree = Tree::from_newick("((A:1,B:2):3,C:4);").unwrap();
            tree.distance_matrix().unwrap();
            assert!(tree
                .nodes
                .iter()
                .any(|n| n.subtree_distances.borrow().is_some()));
            if target_diameter {
                tree.rescale_diameter(18.0).unwrap();
            } else {
                tree.rescale(2.0);
            }
            assert!(tree
                .nodes
                .iter()
                .all(|n| n.subtree_distances.borrow().is_none()));
            assert_eq!(tree.diameter().unwrap(), 18.0);
        }
    }
}
