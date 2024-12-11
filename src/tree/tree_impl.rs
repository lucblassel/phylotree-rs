use accurate::sum::NaiveSum;
use accurate::traits::*;
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use ptree::{print_tree, TreeBuilder};
use rand::seq::SliceRandom;
use std::collections::VecDeque;
use std::iter::zip;
use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    fs,
    path::Path,
};

use thiserror::Error;

use super::node::{Node, NodeError};
use super::{EdgeDepth, EdgeLength, NewickFormat, NodeId};

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
    /// General error
    #[error("Encountered an error: {0}")]
    GeneralError(&'static str),
}

/// Errors that can occur when parsing newick files.
#[derive(Error, Debug)]
pub enum NewickParseError {
    /// There is whitespace in one of the branch lengths
    #[error("Cannot have whitespace in number field.")]
    WhiteSpaceInNumber,
    /// There is an unclosed bracket in the newick String
    #[error("Missing a closing bracket.")]
    UnclosedBracket,
    /// The newick string is missing a final semi-colon
    #[error("The tree is missin a semi colon at the end.")]
    NoClosingSemicolon,
    /// We are trying to close a subtre but have no parent node.
    #[error("Parent node of subtree not found")]
    NoSubtreeParent,
    /// There was a [`TreeError`] when building a tree fromthe newick string
    #[error("Problem with building the tree.")]
    TreeError(#[from] TreeError),
    /// There was a [`std::num::ParseFloatError`] when parsing branch lengths
    #[error("Could not parse a branch length")]
    FloatError(#[from] std::num::ParseFloatError),
    /// There was a [`std::io::Error`] when reading a newick file
    #[error("Problem reading file")]
    IoError(#[from] std::io::Error),
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
    /// Khuner Felsenstein Branch Score
    pub branch_score: f64,
}

/// Used to hold compared tree edges
type EdgeCompare = (
    Vec<(EdgeDepth, EdgeLength)>,
    Vec<(EdgeDepth, EdgeLength)>,
    Vec<((EdgeDepth, EdgeLength), (EdgeDepth, EdgeLength))>,
);

type Partition = FixedBitSet;
type WrappedPartitionMap = HashMap<Partition, (usize, Option<EdgeLength>)>;
type PartitionMap = HashMap<Partition, (EdgeDepth, EdgeLength)>;
type PartitionSet = HashSet<Partition>;

/// A Phylogenetic tree
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    leaf_index: RefCell<Option<Vec<String>>>,
    partitions: RefCell<Option<WrappedPartitionMap>>,
}

/// Base methods to add and get [`Node`] objects to and from the [`Tree`].
///   
/// ----
/// ----
impl Tree {
    /// Create a new empty Tree object
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            leaf_index: RefCell::new(None),
            partitions: RefCell::new(None),
        }
    }

    // ############################
    // # ADDING AND GETTING NODES #
    // ############################

    /// Add a new node to the tree.
    pub fn add(&mut self, node: Node) -> NodeId {
        let idx = self.nodes.len();
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
    /// assert_eq!(tree.get(&root_id).unwrap().children.len(), 2);
    ///
    /// // The depths of child nodes are derived from the parent node
    /// assert_eq!(tree.get(&left).unwrap().get_depth(), 1);
    /// assert_eq!(tree.get(&right).unwrap().get_depth(), 1);
    ///
    /// // If an edge length is specified then it is set in both child and parent
    /// assert_eq!(tree.get(&right).unwrap().parent_edge, Some(0.1));
    /// assert_eq!(tree.get(&root_id).unwrap().get_child_edge(&right), Some(0.1));
    /// ```
    pub fn add_child(
        &mut self,
        node: Node,
        parent: NodeId,
        edge: Option<EdgeLength>,
    ) -> Result<NodeId, TreeError> {
        if parent >= self.nodes.len() {
            return Err(TreeError::NodeNotFound(parent));
        }

        let mut node = node;

        node.set_parent(parent, edge);
        node.set_depth(self.get(&parent)?.depth + 1);

        let id = self.add(node);

        self.get_mut(&id)?.set_id(id);
        self.get_mut(&parent)?.add_child(id, edge);

        Ok(id)
    }

    /// Get a reference to a specific Node of the tree
    pub fn get(&self, id: &NodeId) -> Result<&Node, TreeError> {
        if *id >= self.nodes.len() {
            return Err(TreeError::NodeNotFound(*id));
        }
        let node = &self.nodes[*id];
        if node.deleted {
            return Err(TreeError::NodeNotFound(*id));
        }

        Ok(node)
    }

    /// Get a mutable reference to a specific Node of the tree
    pub fn get_mut(&mut self, id: &NodeId) -> Result<&mut Node, TreeError> {
        if *id >= self.nodes.len() {
            return Err(TreeError::NodeNotFound(*id));
        }
        let node = &mut self.nodes[*id];
        if node.deleted {
            return Err(TreeError::NodeNotFound(*id));
        }

        Ok(node)
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
    /// assert_eq!(tree.get_by_name("child"), Some(tree.get(&child_idx).unwrap()));
    /// ```
    pub fn get_by_name(&self, name: &str) -> Option<&Node> {
        self.nodes
            .iter()
            .find(|node| node.name.is_some() && node.name == Some(String::from(name)))
    }

    /// Get a mutable reference to a node in the tree by name
    pub fn get_by_name_mut(&mut self, name: &str) -> Option<&mut Node> {
        self.nodes
            .iter_mut()
            .find(|node| node.name.is_some() && node.name == Some(String::from(name)))
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
    pub fn get_root(&self) -> Result<NodeId, TreeError> {
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
    /// assert_eq!(tree.get_leaves(), vec![left, right]);
    /// ```
    pub fn get_leaves(&self) -> Vec<NodeId> {
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
    /// let names: Vec<_> = tree.get_leaf_names()
    ///     .into_iter()
    ///     .flatten()
    ///     .collect();
    /// assert_eq!(names, vec!["left", "right"]);
    /// ```
    pub fn get_leaf_names(&self) -> Vec<Option<String>> {
        self.get_leaves()
            .iter()
            .map(|leaf_id| self.get(leaf_id).unwrap().name.clone())
            .collect()
    }

    /// Gets the node ids of all the nodes in the subtree rooted at the specified node
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let subtree: Vec<_> = tree.get_subtree(&sub_root.id)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(subtree, vec!["E", "C", "D"])
    /// ```
    pub fn get_subtree(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![*root];

        for child in self.get(root)?.children.iter() {
            indices.extend(self.get_subtree(child)?);
        }

        Ok(indices)
    }

    /// Gets the node ids of all the nodes in the subtree rooted at the specified node, except the root
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let subtree: Vec<_> = tree.get_descendants(&sub_root.id)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(subtree, vec!["C", "D"])
    /// ```
    pub fn get_descendants(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];

        for child in self.get(root)?.children.iter() {
            indices.extend(self.get_subtree(child)?);
        }

        Ok(indices)
    }

    /// Gets the node ids of all the leaves in the subtree rooted at the specified node
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let sub_leaves: Vec<_> = tree.get_subtree_leaves(&sub_root.id)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(sub_leaves, vec!["C", "D"])
    /// ```
    pub fn get_subtree_leaves(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        Ok(self
            .get_subtree(root)?
            .into_iter()
            .filter(|id| self.get(id).unwrap().is_tip())
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
    /// let preorder: Vec<_> = tree.preorder(&tree.get_root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(preorder, vec!["F", "B", "A", "D", "C", "E", "G", "I", "H"])
    /// ```
    pub fn preorder(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![*root];
        for child in self.get(root)?.children.iter() {
            indices.extend(self.preorder(child)?)
        }

        Ok(indices)
    }

    /// Returns a vector containing node ids in the same order as the
    /// [postorder](https://en.wikipedia.org/wiki/Tree_traversal#Post-order,_LRN ) tree traversal
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let postorder: Vec<_> = tree.postorder(&tree.get_root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(postorder, vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"])
    /// ```
    pub fn postorder(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];
        for child in self.get(root)?.children.iter() {
            indices.extend(self.postorder(child)?)
        }
        indices.push(*root);

        Ok(indices)
    }

    /// Returns a vector containing node ids in the same order as the
    /// [inorder](https://en.wikipedia.org/wiki/Tree_traversal#In-order,_LNR) tree traversal.
    /// This assumes that the tree is binary.
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let inorder: Vec<_> = tree.inorder(&tree.get_root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(inorder, vec!["A", "B", "C", "D", "E", "F", "H", "I", "G"])
    /// ```
    /// *N.B.: the resulting traversal in the example is different than the one in the
    /// [wikipedia](https://en.wikipedia.org/wiki/Tree_traversal) page. That is because
    /// the [`Node`] struct cannot have a right child only, so in our tree I is the left
    /// child of G.*
    pub fn inorder(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];
        let children = &self.get(root)?.children;

        // Tree is not binary
        if children.len() > 2 {
            return Err(TreeError::IsNotBinary);
        }

        // There is a left child
        if !children.is_empty() {
            indices.extend(self.inorder(&children[0])?)
        }

        indices.push(*root);

        // There is a right child
        if children.len() > 1 {
            indices.extend(self.inorder(&children[1])?)
        }

        Ok(indices)
    }

    /// Returns a vector containing node ids in the same order as the
    /// [levelorder](https://en.wikipedia.org/wiki/Tree_traversal#Breadth-first_search) tree traversal
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let levelorder: Vec<_> = tree.levelorder(&tree.get_root().unwrap())
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(levelorder, vec!["F", "B", "G", "A", "D", "I", "C", "E", "H"])
    /// ```
    pub fn levelorder(&self, root: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut indices = vec![];
        let mut queue = VecDeque::new();
        queue.push_back(root);
        while !queue.is_empty() {
            let root = queue.pop_front().unwrap();
            indices.push(*root);
            for child in self.get(root)?.children.iter() {
                queue.push_back(child)
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
    pub fn is_rooted(&self) -> Result<bool, TreeError> {
        let root_id = self.get_root()?;

        Ok(!self.nodes.is_empty() && self.get(&root_id)?.children.len() == 2)
    }

    /// Checks if all the tips have unique names (This check assumes that all tips have a name)
    pub fn has_unique_tip_names(&self) -> Result<bool, TreeError> {
        let mut names = HashSet::new();
        for name in self.get_leaf_names() {
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

        let root = self.get_root()?;

        self.get_leaves()
            .iter()
            .map(|leaf| {
                let (edge_sum, num_edges) = self.get_distance(&root, leaf).unwrap();
                match edge_sum {
                    Some(height) => height,
                    None => num_edges as f64,
                }
            })
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .ok_or(TreeError::IsEmpty)
    }

    /// Returns the diameter of the tree
    /// (i.e. longest tip to tip distance)
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
        self.get_leaves()
            .iter()
            .combinations(2)
            .map(|pair| {
                let (edge_sum, num_edges) = self.get_distance(pair[0], pair[1]).unwrap();
                match edge_sum {
                    Some(height) => height,
                    None => num_edges as f64,
                }
            })
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .ok_or(TreeError::IsEmpty)
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
    pub fn cherries(&self) -> Result<usize, TreeError> {
        if !self.is_binary()? {
            return Err(TreeError::IsNotBinary);
        }
        if !self.nodes.is_empty() {
            let mut n = 0;
            for node in self.nodes.iter() {
                if node.children.len() == 2
                    && self.get(&node.children[0])?.is_tip()
                    && self.get(&node.children[1])?.is_tip()
                {
                    n += 1;
                }
            }
            Ok(n)
        } else {
            Err(TreeError::IsEmpty)
        }
    }

    /// Computes the Colless index for the tree.
    /// The colless index, $I_c$, measures the imbalance of a phylogenetic tree:  
    /// $$
    /// I_c = \sum_{i \in nodes} |L_i - R_i|
    /// $$
    ///
    /// Where $L_i$ is the number of leaves in the left subtree of node $i$ and
    /// $R_i$ the number of leaves in the right subtree of $i$.
    ///
    pub fn colless(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;

        let mut colless = 0;

        for node in self.nodes.iter().filter(|node| !node.is_tip()) {
            let left = self.get_subtree_leaves(&node.children[0])?.len();
            let right = if node.children.len() > 1 {
                self.get_subtree_leaves(&node.children[1])?.len()
            } else {
                0
            };

            colless += left.abs_diff(right);
        }

        Ok(colless)
    }

    /// Computes the normalized colless statistic with a Yule null model:  
    /// $$
    /// I_{yule} = \frac{I_c - n\cdot\log(n) - n(\gamma-1-\log(2))}{n}
    /// $$
    /// Where $I_c$ is the unnormalized colless index *(computed with [`Tree::colless()`])*,
    /// $n$ the number of leaves
    /// and $\gamma$ the [Euler constant](https://en.wikipedia.org/wiki/Euler%27s_constant).  
    /// *([see also apTreeshape](https://search.r-project.org/CRAN/refmans/apTreeshape/html/colless.html))*
    pub fn colless_yule(&self) -> Result<f64, TreeError> {
        self.colless().map(|i_c| {
            let n = self.n_leaves() as f64;
            let e_i_c = n * n.ln() + (0.57721566 - 1. - f64::ln(2.0)) * n;

            (i_c as f64 - e_i_c) / n
        })
    }

    /// Computes the normalized colless statistic with a PDA null model:  
    /// $$
    /// I_{PDA} = \frac{I_c}{n^{3/2}}
    /// $$
    /// Where $I_c$ is the unnormalized colless index *(computed with [`Tree::colless()`])*
    /// and $n$ the number of leaves.  
    /// *([see also apTreeshape](https://search.r-project.org/CRAN/refmans/apTreeshape/html/colless.html))*
    pub fn colless_pda(&self) -> Result<f64, TreeError> {
        self.colless()
            .map(|i_c| i_c as f64 / f64::powf(self.n_leaves() as f64, 3.0 / 2.0))
    }

    /// Computes the Sackin index. The Sackin index, $I_s$, is computed by taking the
    /// sum over all internal nodes of the number of leaves descending from that node.
    /// A smaller Sackin index means a more balanced tree.
    pub fn sackin(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;

        Ok(self
            .get_leaves()
            .iter()
            .map(|tip_idx| self.get(tip_idx).unwrap().depth)
            .sum())
    }

    /// Computes the normalized Sackin index with a Yule null model:
    /// $$
    /// I_{yule} = \frac{I_s - 2n\cdot \sum_{j=2}^n \frac{1}{j}}{n}
    /// $$
    /// With $I_s$ the unnormalized Sackin index *(computed with [`Tree::sackin()`])*
    /// and $n$ the number of leaves in the tree.  
    /// *([see also apTreeshape](https://search.r-project.org/CRAN/refmans/apTreeshape/html/sackin.html))*
    pub fn sackin_yule(&self) -> Result<f64, TreeError> {
        self.sackin().map(|i_n| {
            let n = self.n_leaves();
            let sum: f64 = (2..=n).map(|i| 1.0 / (i as f64)).sum();

            (i_n as f64 - 2.0 * (n as f64) * sum) / n as f64
        })
    }

    /// Computes the normalized sackin statistic with a PDA null model:
    /// $$
    /// I_{PDA} = \frac{I_s}{n^{3/2}}
    /// $$
    /// With $I_s$ the unnormalized Sackin index *(computed with [`Tree::sackin()`])*
    /// and $n$ the number of leaves in the tree.  
    /// *([see also apTreeshape](https://search.r-project.org/CRAN/refmans/apTreeshape/html/sackin.html))*
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

    /// Initializes the leaf index
    fn init_leaf_index(&self) -> Result<(), TreeError> {
        if self.nodes.is_empty() {
            return Err(TreeError::IsEmpty);
        }
        if self.leaf_index.borrow().is_some() {
            return Ok(());
        }

        let names = self.get_leaf_names();
        if names.len() != self.n_leaves() {
            return Err(TreeError::UnnamedLeaves);
        }

        if !self.has_unique_tip_names()? {
            return Err(TreeError::DuplicateLeafNames);
        }

        (*self.leaf_index.borrow_mut()) = Some(names.into_iter().flatten().sorted().collect());

        Ok(())
    }

    /// Get the partition corresponding to the branch associated to the node at index
    fn get_partition(&self, index: &NodeId) -> Result<Partition, TreeError> {
        self.init_leaf_index()?;

        let subtree_leaves = self.get_subtree_leaves(index)?;
        let l_index = self.leaf_index.borrow();
        let indices = subtree_leaves
            .iter()
            .filter_map(|index| self.get(index).unwrap().name.as_ref())
            .map(|name| l_index.iter().flatten().position(|n| n == name).unwrap());

        let mut bitset = FixedBitSet::with_capacity(self.n_leaves());
        for index in indices {
            bitset.insert(index);
        }

        let mut toggled = bitset.clone();
        toggled.toggle_range(..);

        Ok(toggled.min(bitset))
    }

    /// Helper function to view a partition as
    pub fn partition_to_leaves(&self, partition: &Partition) -> Result<String, TreeError> {
        self.init_leaf_index()?;

        let v = self.leaf_index.borrow().clone().unwrap();
        Ok(partition.ones().map(|i| v[i].clone()).collect())
    }

    /// Caches partitions for distance computation
    fn init_partitions(&self) -> Result<(), TreeError> {
        self.init_leaf_index()?;

        if self.partitions.borrow().is_some() {
            return Ok(());
        }

        let mut partitions: WrappedPartitionMap = HashMap::new();

        for node in self
            .nodes
            .iter()
            .filter(|n| !(n.deleted || n.parent.is_none() || n.is_tip()))
        {
            let part = self.get_partition(&node.id)?;

            if part.count_ones(..) == 1 {
                continue;
            }

            let new_len = node.parent_edge;
            let old_value = partitions.get(&part);

            let len = match (new_len, old_value) {
                (None, None) => None,
                (Some(new_len), Some((_, old_len))) => old_len.map(|v| v + new_len),
                (Some(new_len), None) => Some(new_len),
                (None, Some((_, old_len))) => *old_len,
            };

            partitions.insert(part, (node.depth, len));
        }

        (*self.partitions.borrow_mut()) = Some(partitions);

        Ok(())
    }

    /// Get all partitions of a tree
    pub fn get_partitions(&self) -> Result<PartitionSet, TreeError> {
        self.init_leaf_index()?;
        self.init_partitions()?;

        Ok(HashSet::from_iter(
            self.partitions
                .borrow()
                .as_ref()
                .unwrap()
                .iter()
                .map(|(k, _)| k.clone()),
        ))
    }

    /// Get all partitions of a tree along with corresponding branch lengths and branch depths
    pub(crate) fn get_partitions_with_lengths(&self) -> Result<PartitionMap, TreeError> {
        self.init_leaf_index()?;
        self.init_partitions()?;

        let mut partitions = HashMap::new();
        for (bitset, (depth, len)) in self.partitions.borrow().as_ref().unwrap().iter() {
            let len = len.ok_or(TreeError::MissingBranchLengths)?;
            partitions.insert(bitset.clone(), (*depth, len));
        }

        Ok(partitions)
    }

    /// Empties the partitions cache
    fn reset_partitions(&mut self) {
        (*self.partitions.borrow_mut()) = None;
    }

    /// Empties the leaf index
    fn reset_leaf_index(&mut self) {
        (*self.leaf_index.borrow_mut()) = None
    }

    /// Resets the caches used when computing bipartitions
    /// *(i.e. with [`Tree::compare_topologies()`])*.
    /// You should call this if you have computed bipartitions in the tree
    /// and then changed the tree.
    pub fn reset_bipartition_cache(&mut self) {
        self.reset_leaf_index();
        self.reset_partitions();
    }

    // #################
    // # COMPARE TREES #
    // #################

    /// Computes the [Robinson Foulds distance](https://en.wikipedia.org/wiki/Robinson–Foulds_metric)
    /// [(Robinson & Foulds, 1981)](https://doi.org/10.1016/0025-5564(81)90043-2)
    /// between two trees. The RF distance is defined as the number of unique bipartitions for each tree:
    /// $$
    /// RF = |A\cup B| - |A\cap B|
    /// $$
    /// Where $A$ and $B$ are the sets of bipartitions of the first and second trees.  
    /// See also [Tree::compare_topologies()]
    pub fn robinson_foulds(&self, other: &Self) -> Result<usize, TreeError> {
        let partitions_s = self.get_partitions()?;
        let partitions_o = other.get_partitions()?;

        if *(self.leaf_index.borrow()) != *(other.leaf_index.borrow()) {
            return Err(TreeError::DifferentTipIndices);
        }

        let mut root_s = HashSet::new();
        for i in self.get(&self.get_root()?)?.children.iter() {
            root_s.insert(self.get_partition(i)?);
        }
        let mut root_o = HashSet::new();
        for i in other.get(&other.get_root()?)?.children.iter() {
            root_o.insert(other.get_partition(i)?);
        }

        let same_root = root_s == root_o;

        let i = partitions_o.intersection(&partitions_s).count();
        let rf = partitions_o.len() + partitions_s.len() - 2 * i;

        // Hacky...
        if self.is_rooted()? && rf != 0 && !same_root {
            Ok(rf + 2)
        } else {
            Ok(rf)
        }
    }

    /// Computes the normalized Robinson Foulds distance between two trees
    /// [(Robinson & Foulds, 1981)](https://doi.org/10.1016/0025-5564(81)90043-2).
    /// The RF distance is normalized by the maximum possible RF distance for both trees
    /// *(i.e the number of bipartitions in both trees)* so that the resulting distance
    /// is contained within [0, 1]:  
    /// $$
    /// RF_{norm} = \frac{RF}{|A| + |B|}
    /// $$
    /// Where $A$ and $B$ are the sets of bipartitions of the first and second trees.  
    /// See also [Tree::compare_topologies()]
    pub fn robinson_foulds_norm(&self, other: &Self) -> Result<f64, TreeError> {
        let rf = self.robinson_foulds(other)?;

        let partitions_s = self.get_partitions()?;
        let partitions_o = other.get_partitions()?;

        let tot = partitions_o.len() + partitions_s.len();

        Ok((rf as f64) / (tot as f64))
    }

    /// Computes the weighted Robinson Foulds distance between two trees
    /// [(Robinson & Foulds, 1979)](https://doi.org/10.1007/BFb0102690).
    /// This distance is equal to the absolute difference of branch lengths for
    /// matched bipartitions between the two trees, plus branch lenghts for unique bipartitions:  
    /// $$
    /// RF_{weighted} = \sum_{e \in A\cap B} |d_{(e,A)} - d_{(e,B)}| +
    /// \sum_{e \in A\setminus B}d_{(e,A)} +
    /// \sum_{e \in B\setminus A}d_{(e,B)}
    /// $$
    /// Where $A$ and $B$ are the sets of bipartitions of the first and second trees,
    /// and $d_{(e,A)}$ the branch length of bipartition $e$ in the first tree ($A$).  
    /// See also [Tree::compare_topologies()]
    pub fn weighted_robinson_foulds(&self, other: &Self) -> Result<f64, TreeError> {
        let partitions_s = self.get_partitions_with_lengths()?;
        let partitions_o = other.get_partitions_with_lengths()?;

        let mut dist = 0.;

        for (edge, (_, len_s)) in partitions_s.iter() {
            if let Some((_, len_o)) = partitions_o.get(edge) {
                dist += (len_s - len_o).abs()
            } else {
                dist += len_s
            }
        }

        for (edge, (_, len_o)) in partitions_o.iter() {
            if !partitions_s.contains_key(edge) {
                dist += len_o
            }
        }

        Ok(dist)
    }

    /// Computes the khuner felsenstein branch score between two trees,
    /// [(Khuner & Felsenstein, 1994)](https://doi.org/10.1093/oxfordjournals.molbev.a040126).
    /// The distance is computed by taking the squared difference of branch lengths for
    /// matched bipartitions between the two trees, plus squared branch lenghts for unique bipartitions.
    /// The branch score is then derived by taking the square root of that total sum:  
    /// $$
    /// KF = \sqrt{
    ///     \sum_{e \in A\cap B} (d_{(e,A)} - d_{(e,B)})^2 +
    ///     \sum_{e \in A\setminus B}d_{(e,A)}^2 +
    ///     \sum_{e \in B\setminus A}d_{(e,B)}^2
    /// }
    /// $$
    /// See also [Tree::compare_topologies()]
    pub fn khuner_felsenstein(&self, other: &Self) -> Result<f64, TreeError> {
        let partitions_s = self.get_partitions_with_lengths()?;
        let partitions_o = other.get_partitions_with_lengths()?;

        let mut dist = 0.;

        for (edge, (_, len_s)) in partitions_s.iter() {
            if let Some((_, len_o)) = partitions_o.get(edge) {
                dist += f64::powi(len_s - len_o, 2)
            } else {
                dist += f64::powi(*len_s, 2)
            }
        }

        for (edge, (_, len_o)) in partitions_o.iter() {
            if !partitions_s.contains_key(edge) {
                dist += f64::powi(*len_o, 2)
            }
        }

        Ok(dist.sqrt())
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
    /// let rf = tree1.robinson_foulds(&tree2).unwrap() as f64;
    /// let norm_rf = tree1.robinson_foulds_norm(&tree2).unwrap();
    /// let weighted_rf = tree1.weighted_robinson_foulds(&tree2).unwrap();
    /// let branch_score = tree1.khuner_felsenstein(&tree2).unwrap();
    ///
    /// let comparison = tree1.compare_topologies(&tree2).unwrap();
    ///
    /// assert_eq!(rf, comparison.rf);
    /// assert_eq!(norm_rf, comparison.norm_rf);
    /// assert_eq!(weighted_rf, comparison.weighted_rf);
    /// assert_eq!(branch_score, comparison.branch_score);
    /// ```
    pub fn compare_topologies(&self, other: &Self) -> Result<Comparison, TreeError> {
        let partitions_s = self.get_partitions_with_lengths()?;
        let partitions_o = other.get_partitions_with_lengths()?;

        let tot = partitions_o.len() + partitions_s.len();

        let mut intersection = 0.;
        let mut rf_weight = 0.;
        let mut kf = 0.;

        for (edge, (_, len_s)) in partitions_s.iter() {
            if let Some((_, len_o)) = partitions_o.get(edge) {
                rf_weight += (len_s - len_o).abs();
                kf += f64::powi(len_s - len_o, 2);
                intersection += 1.0;
            } else {
                rf_weight += len_s;
                kf += f64::powi(*len_s, 2);
            }
        }

        for (edge, (_, len_o)) in partitions_o.iter() {
            if !partitions_s.contains_key(edge) {
                rf_weight += len_o;
                kf += f64::powi(*len_o, 2);
            }
        }

        let mut rf = tot as f64 - 2.0 * intersection;

        // Hacky...
        let mut root_s = HashSet::new();
        for i in self.get(&self.get_root()?)?.children.iter() {
            root_s.insert(self.get_partition(i)?);
        }
        let mut root_o = HashSet::new();
        for i in other.get(&other.get_root()?)?.children.iter() {
            root_o.insert(other.get_partition(i)?);
        }
        let same_root = root_s == root_o;
        if self.is_rooted()? && other.is_rooted()? && rf != 0.0 && !same_root {
            rf += 2.0;
        }

        Ok(Comparison {
            rf,
            norm_rf: rf / tot as f64,
            weighted_rf: rf_weight,
            branch_score: kf.sqrt(),
        })
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

        Ok(HashMap::from_iter(self.get_leaves().iter().map(|idx| {
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
    /// let path: Vec<_> = tree.get_path_from_root(&5)
    ///     .unwrap()
    ///     .iter()
    ///     .map(|id| tree.get(id).unwrap().name.clone())
    ///     .flatten()
    ///     .collect();
    ///
    /// assert_eq!(path, vec!["F", "B", "D", "E"])
    /// ```
    pub fn get_path_from_root(&self, node: &NodeId) -> Result<Vec<NodeId>, TreeError> {
        let mut path = vec![];
        let mut current_node = *node;
        loop {
            path.push(current_node);
            match self.get(&current_node)?.parent {
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
    ///     &tree.get_by_name("A").unwrap().id,
    ///     &tree.get_by_name("D").unwrap().id,
    /// ).unwrap();
    ///
    /// assert_eq!(tree.get(&ancestor).unwrap().name, Some("B".to_owned()))
    /// ```
    pub fn get_common_ancestor(
        &self,
        source: &NodeId,
        target: &NodeId,
    ) -> Result<usize, TreeError> {
        if source == target {
            return Ok(*source);
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
    ///     &tree.get_by_name("A").unwrap().id,
    ///     &tree.get_by_name("I").unwrap().id,
    /// ).unwrap();
    ///
    /// assert_eq!(num_edges, 4);
    /// assert!(sum_edge_lengths.is_none());
    /// ```
    pub fn get_distance(
        &self,
        source: &NodeId,
        target: &NodeId,
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
            for node in list.iter().skip(cursor) {
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
        current: &NodeId,
        prev: Option<&NodeId>,
        lengths: &mut [f64],
        currlength: f64,
    ) -> Result<(), TreeError> {
        if prev.is_some() && self.get(current)?.is_tip() {
            lengths[*current] = currlength;
            return Ok(());
        }

        let children = self.get(current)?.children.clone();
        let mut neighbors: Vec<_> = children
            .iter()
            .map(|idx| (*idx, self.get(idx).unwrap().parent_edge))
            .collect();

        if let Some(parent) = self.get(current)?.parent {
            neighbors.push((parent, self.get(current).unwrap().parent_edge))
        }

        for (neighbor, brlen) in neighbors {
            if Some(&neighbor) != prev {
                if let Some(brlen) = brlen {
                    self.distance_matrix_recursive_impl(
                        &neighbor,
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
        let mut matrix = DistanceMatrix::new_with_size(self.n_leaves());
        let mut cache: Vec<Vec<_>> = vec![vec![f64::INFINITY; size]; size];

        self.init_leaf_index()?;
        let taxa = self.leaf_index.borrow().as_ref().unwrap().clone();
        matrix.set_taxa(taxa)?;

        for tip in self.get_leaves().iter() {
            self.distance_matrix_recursive_impl(tip, None, &mut cache[*tip], 0.0)?
        }

        for pair in self.get_leaves().iter().combinations(2) {
            let (i1, i2) = (pair[0], pair[1]);
            let d = cache[*i1][*i2];
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
        let mut leaf_order = self.get_leaves();
        leaf_order.sort_by(|a, b| self.get(a).unwrap().name.cmp(&self.get(b).unwrap().name));

        let n = self.n_leaves();
        let mut pairwise_vec = vec![NaiveSum::zero(); n * (n - 1) / 2];

        let leaf_idx_to_leaf_order = self
            .nodes
            .iter()
            .enumerate()
            .map(|(idx, _)| leaf_order.iter().position(|&v| v == idx))
            .collect_vec();

        // Converts the node index of a leaf to its index in the leaf_order array
        let get_leaf_index = |leaf: usize| -> Result<usize, TreeError> {
            leaf_idx_to_leaf_order[leaf].ok_or(TreeError::NodeNotFound(leaf))
        };

        for current_node in self.levelorder(&self.get_root()?)?.iter().rev() {
            let mut node_cache: HashMap<_, _, BuildIdentityHasher> = HashMap::default();

            let parent = self.get(current_node)?;
            if parent.is_tip() {
                node_cache.insert(*current_node, 0.);
            }

            // Compute distances from current node to descendant leaves
            for child in parent.children.iter() {
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
                    node_cache.insert(*leaf, len);
                }
            }

            // Compute distances between leaves
            for subtree_roots in parent.children.iter().combinations(2) {
                let subtree1 = self.get(subtree_roots[0])?;
                let subtree2 = self.get(subtree_roots[1])?;

                for (leaf1, _) in subtree1
                    .subtree_distances
                    .borrow()
                    .as_ref()
                    .ok_or(TreeError::MissingBranchLengths)?
                    .iter()
                {
                    for (leaf2, _) in subtree2
                        .subtree_distances
                        .borrow()
                        .as_ref()
                        .ok_or(TreeError::MissingBranchLengths)?
                        .iter()
                    {
                        let distance1 = node_cache.get(leaf1).unwrap();
                        let distance2 = node_cache.get(leaf2).unwrap();

                        let mut i = get_leaf_index(*leaf1)?;
                        let mut j = get_leaf_index(*leaf2)?;
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
                .map(|i| self.get(i).unwrap().clone().name.unwrap())
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
    /// tree.prune(&root_idx);
    ///
    /// assert_eq!(tree.to_newick().unwrap(), "((A,(C,E)D)B)F;")
    /// ```
    pub fn prune(&mut self, root: &NodeId) -> Result<(), TreeError> {
        for child in self.get(root)?.children.clone() {
            self.prune(&child)?
        }

        if let Some(parent) = self.get(root)?.parent {
            self.get_mut(&parent)?.remove_child(root)?;
        }

        self.get_mut(root)?.delete();

        Ok(())
    }

    // Removes a single node
    fn compress_node(&mut self, id: &NodeId) -> Result<(), TreeError> {
        let node = self.get(id)?;

        if node.parent.is_none() || node.children.len() != 1 {
            return Err(TreeError::CouldNotCompressNode(*id));
        }

        let parent = node.parent.unwrap();
        let child = node.children[0];
        let to_remove = node.id;

        let parent_edge = node.parent_edge;
        let child_edge = node.get_child_edge(&child);

        let new_edge = match (parent_edge, child_edge) {
            (Some(p), Some(c)) => Some(p + c),
            (None, None) => None,
            _ => return Err(TreeError::MissingBranchLengths),
        };

        self.get_mut(&child)?.set_parent(parent, new_edge);
        self.get_mut(&parent)?.add_child(child, new_edge);
        self.get_mut(&parent)?.remove_child(&to_remove)?;

        self.get_mut(&to_remove)?.delete();

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
            .cloned()
            .map(|node| node.id)
            .collect();

        for id in to_compress {
            self.compress_node(&id)?;
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
            node.rescale_edges(factor)
        }
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
                let mut children = self.get(&node_id)?.children.clone();
                children.shuffle(rng);

                let parent = self.add_child(Node::new(), node_id, Some(0.0))?;

                for _ in 0..2 {
                    let child = children.pop().unwrap();
                    let edge = self.get(&child)?.parent_edge;
                    self.get_mut(&parent)?.add_child(child, edge);
                    self.get_mut(&child)?.set_parent(parent, edge);
                    self.get_mut(&node_id)?.remove_child(&child)?;
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
        let root = self.get_root()?;
        // Go from tips to root
        for node_id in self.levelorder(&root)?.into_iter().rev() {
            let node = self.get_mut(&node_id)?;
            for child in node.children.iter() {
                descendant_counter[node_id] += descendant_counter[*child] + 1;
            }
            node.children.sort_by_key(|v| descendant_counter[*v]);
        }

        Ok(())
    }

    // recusrive implementation of depth recomputation
    fn reset_depth_impl(&mut self, root: &NodeId, depth: usize) -> Result<(), TreeError> {
        let root = self.get_mut(root)?;
        root.set_depth(depth);

        for &child in root.children.clone().iter() {
            self.reset_depth_impl(&child, depth + 1)?
        }

        Ok(())
    }

    /// Recompute node depths and set them correctly.
    pub fn reset_depths(&mut self) -> Result<(), TreeError> {
        let root = self.get_root()?;
        self.reset_depth_impl(&root, 0)
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
    /// tree.merge_children(&a, &b, None, None, None, Some("D".into()));
    ///
    /// let expected = Tree::from_newick("((A,B)D, C);").unwrap();
    /// assert_eq!(tree.robinson_foulds(&expected).unwrap(), 0);
    /// ```
    pub fn merge_children(
        &mut self,
        child1: &NodeId,
        child2: &NodeId,
        edge1: Option<EdgeLength>,
        edge2: Option<EdgeLength>,
        parent_edge: Option<EdgeLength>,
        parent_name: Option<String>,
    ) -> Result<NodeId, TreeError> {
        // Check that nodes are siblings
        let parent = self.get(child1)?.parent;
        if parent != self.get(child2)?.parent {
            return Err(TreeError::MergingNonSiblingNodes(*child1, *child2));
        }

        // Add new parent node as child of current parent
        let parent = match parent {
            Some(parent_id) => {
                // Remove merged nodes as children of current parent
                let parent_node = self.get_mut(&parent_id)?;
                parent_node.remove_child(child1)?;
                parent_node.remove_child(child2)?;
                // Add new parent
                self.add_child(Node::new(), parent_id, parent_edge)?
            }
            None => self.add(Node::new()),
        };

        // Set parent/child relationships between merged nodes and new parent node
        let p = self.get_mut(&parent)?;
        p.add_child(*child1, edge1);
        p.add_child(*child2, edge2);
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

    /// Generate newick representation of tree
    fn to_newick_impl(&self, root: &NodeId, format: NewickFormat) -> Result<String, TreeError> {
        let root = self.get(root)?;
        if root.children.is_empty() {
            Ok(root.to_newick(format))
        } else {
            Ok("(".to_string()
                + &(root
                    .children
                    .iter()
                    .map(|child_idx| self.to_newick_impl(child_idx, format).unwrap()))
                .collect::<Vec<String>>()
                .join(",")
                + ")"
                + &(root.to_newick(format)))
        }
    }

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
        let root = self.get_root()?;
        Ok(self.to_newick_impl(&root, NewickFormat::AllFields)? + ";")
    }

    /// Writes the tree as a newick formatted string with a specified
    /// output format from [`NewickFormat`].
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
    pub fn to_formatted_newick(&self, format: NewickFormat) -> Result<String, TreeError> {
        let root = self.get_root()?;
        Ok(self.to_newick_impl(&root, format)? + ";")
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
        #[derive(Debug, PartialEq)]
        enum Field {
            Name,
            Length,
            Comment,
        }

        let mut tree = Tree::new();

        let mut parsing = Field::Name;
        let mut current_name: Option<String> = None;
        let mut current_length: Option<String> = None;
        let mut current_comment: Option<String> = None;
        let mut current_index: Option<NodeId> = None;
        let mut parent_stack: Vec<NodeId> = Vec::new();

        let mut open_delimiters = Vec::new();
        let mut within_quotes = false;

        for c in newick.chars() {
            // Add character in quotes to name
            if within_quotes && parsing == Field::Name && c != '"' {
                if let Some(name) = current_name.as_mut() {
                    name.push(c)
                } else {
                    current_name = Some(c.into())
                }
                continue;
            }

            // Add current character to comment
            if parsing == Field::Comment && c != ']' {
                if let Some(comment) = current_comment.as_mut() {
                    comment.push(c)
                } else {
                    current_comment = Some(c.into())
                }
                continue;
            }

            // Skip unquoted whitespace
            if c.is_whitespace() && !within_quotes {
                continue;
            }

            match c {
                '"' => {
                    // Enter or close quoted section (name)
                    // TODO: handle escaped quotes
                    within_quotes = !within_quotes;
                    if parsing == Field::Name {
                        if let Some(name) = current_name.as_mut() {
                            name.push(c)
                        } else {
                            current_name = Some(c.into())
                        }
                    }
                }
                '[' => {
                    parsing = Field::Comment;
                }
                ']' => {
                    parsing = Field::Name;
                }
                '(' => {
                    // Start subtree
                    match parent_stack.last() {
                        None => parent_stack.push(tree.add(Node::new())),
                        Some(parent) => {
                            parent_stack.push(tree.add_child(Node::new(), *parent, None)?)
                        }
                    };
                    open_delimiters.push(0);
                }
                ':' => {
                    // Start parsing length
                    parsing = Field::Length;
                }
                ',' => {
                    // Add sibling
                    let node = if let Some(index) = current_index {
                        tree.get_mut(&index)?
                    } else {
                        if let Some(parent) = parent_stack.last() {
                            current_index = Some(tree.add_child(Node::new(), *parent, None)?);
                        } else {
                            unreachable!("Sould not be possible to have named child with no parent")
                        };
                        tree.get_mut(current_index.as_ref().unwrap())?
                    };

                    if let Some(name) = current_name {
                        node.set_name(name);
                    }

                    let edge = if let Some(length) = current_length {
                        Some(length.parse()?)
                    } else {
                        None
                    };
                    if let Some(parent) = node.parent {
                        node.set_parent(parent, edge);
                    }

                    node.comment = current_comment;

                    current_name = None;
                    current_comment = None;
                    current_length = None;
                    current_index = None;

                    parsing = Field::Name;
                }
                ')' => {
                    // Close subtree
                    open_delimiters.pop();
                    let node = if let Some(index) = current_index {
                        tree.get_mut(&index)?
                    } else {
                        if let Some(parent) = parent_stack.last() {
                            current_index = Some(tree.add_child(Node::new(), *parent, None)?);
                        } else {
                            unreachable!("Sould not be possible to have named child with no parent")
                        };
                        tree.get_mut(current_index.as_ref().unwrap())?
                    };

                    if let Some(name) = current_name {
                        node.set_name(name);
                    }

                    let edge = if let Some(length) = current_length {
                        Some(length.parse()?)
                    } else {
                        None
                    };
                    if let Some(parent) = node.parent {
                        node.set_parent(parent, edge);
                    }

                    node.comment = current_comment;

                    current_name = None;
                    current_comment = None;
                    current_length = None;

                    parsing = Field::Name;

                    if let Some(parent) = parent_stack.pop() {
                        current_index = Some(parent)
                    } else {
                        return Err(NewickParseError::NoSubtreeParent);
                    }
                }
                ';' => {
                    // Finish parsing the Tree
                    if !open_delimiters.is_empty() {
                        return Err(NewickParseError::UnclosedBracket);
                    }
                    let node = tree.get_mut(current_index.as_ref().unwrap())?;
                    node.name = current_name;
                    node.comment = current_comment;
                    if let Some(length) = current_length {
                        node.parent_edge = Some(length.parse()?);
                    }

                    // Finishing pass to make sure that branch lenghts are set in both children and parents
                    let ids: Vec<_> = tree.nodes.iter().map(|node| node.id).collect();
                    for node_id in ids {
                        if let Some(edge) = tree.get(&node_id)?.parent_edge {
                            if let Some(parent) = tree.get(&node_id)?.parent {
                                tree.get_mut(&parent)?.set_child_edge(&node_id, Some(edge));
                            }
                        }
                    }

                    return Ok(tree);
                }
                _ => {
                    // Parse characters in fields
                    match parsing {
                        Field::Name => {
                            if let Some(name) = current_name.as_mut() {
                                name.push(c)
                            } else {
                                current_name = Some(c.into())
                            }
                        }
                        Field::Length => {
                            if c.is_whitespace() {
                                return Err(NewickParseError::WhiteSpaceInNumber);
                            }
                            if let Some(length) = current_length.as_mut() {
                                length.push(c)
                            } else {
                                current_length = Some(c.into())
                            }
                        }
                        Field::Comment => unimplemented!(),
                    };
                }
            }
        }

        Err(NewickParseError::NoClosingSemicolon)
    }

    /// Writes the tree to a newick file
    pub fn to_file(&self, path: &Path) -> Result<(), TreeError> {
        match fs::write(path, self.to_newick()?) {
            Ok(_) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    /// Creates a tree from a newick file
    pub fn from_file(path: &Path) -> Result<Self, NewickParseError> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_newick(&newick_string)
    }

    /// Outputs a Nexus formatted string of the tree
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
    fn print_nodes(
        &self,
        root_idx: &NodeId,
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
            for child_idx in root.children.iter() {
                self.print_nodes(child_idx, output_tree, debug)?;
            }
            output_tree.end_child();
        }

        Ok(())
    }

    /// Print a debug view of the tree to the console
    pub fn print_debug(&self) -> Result<(), TreeError> {
        let root = self.get_root()?;
        let mut builder = TreeBuilder::new(format!("{:?}", root));
        for child_idx in self.get(&root)?.children.iter() {
            self.print_nodes(child_idx, &mut builder, true)?;
        }
        let tree = builder.build();
        print_tree(&tree)?;
        Ok(())
    }

    /// Print the tree to the console
    pub fn print(&self) -> Result<(), TreeError> {
        let root = self.get_root()?;
        let mut builder = TreeBuilder::new(format!("{:?}", root));
        for child_idx in self.get(&root)?.children.iter() {
            self.print_nodes(child_idx, &mut builder, false)?;
        }
        let tree = builder.build();
        print_tree(&tree)?;
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct IdentityHasher(usize);

impl core::hash::Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.0 as u64
    }

    fn write(&mut self, _bytes: &[u8]) {
        unimplemented!("IdentityHasher only supports usize keys")
    }

    fn write_usize(&mut self, i: usize) {
        self.0 = i;
    }
}

type BuildIdentityHasher = core::hash::BuildHasherDefault<IdentityHasher>;

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
        tree.add_child(Node::new_named("B"), 0, None)?; // 1
        tree.add_child(Node::new_named("G"), 0, None)?; // 2
        tree.add_child(Node::new_named("A"), 1, None)?; // 3
        tree.add_child(Node::new_named("D"), 1, None)?; // 4
        tree.add_child(Node::new_named("I"), 2, None)?; // 5
        tree.add_child(Node::new_named("C"), 4, None)?; // 6
        tree.add_child(Node::new_named("E"), 4, None)?; // 7
        tree.add_child(Node::new_named("H"), 5, None)?; // 8

        Ok(tree)
    }

    /// Generates example tree from the newick format wikipedia page
    /// https://en.wikipedia.org/wiki/Newick_format#Examples
    fn build_tree_with_lengths() -> Result<Tree, TreeError> {
        let mut tree = Tree::new();
        tree.add(Node::new_named("F")); // 0
        tree.add_child(Node::new_named("A"), 0, Some(0.1))?; // 1
        tree.add_child(Node::new_named("B"), 0, Some(0.2))?; // 2
        tree.add_child(Node::new_named("E"), 0, Some(0.5))?; // 3
        tree.add_child(Node::new_named("C"), 3, Some(0.3))?; // 4
        tree.add_child(Node::new_named("D"), 3, Some(0.4))?; // 5

        Ok(tree)
    }

    fn build_tree_without_lengths() -> Result<Tree, TreeError> {
        let mut tree = Tree::new();
        tree.add(Node::new_named("F")); // 0
        tree.add_child(Node::new_named("A"), 0, None)?; // 1
        tree.add_child(Node::new_named("B"), 0, None)?; // 2
        tree.add_child(Node::new_named("E"), 0, None)?; // 3
        tree.add_child(Node::new_named("C"), 3, None)?; // 4
        tree.add_child(Node::new_named("D"), 3, None)?; // 5

        Ok(tree)
    }
    fn get_values(indices: &[usize], tree: &Tree) -> Vec<Option<String>> {
        indices
            .iter()
            .map(|idx| tree.get(idx).unwrap().name.clone())
            .collect()
    }

    #[test]
    fn test_tips() {
        let mut tree = Tree::new();
        tree.add(Node::new_named("root"));
        assert_eq!(tree.get_leaves(), vec![0]);

        tree.add_child(Node::new_named("A"), 0, Some(0.1)).unwrap(); // 1
        tree.add_child(Node::new_named("B"), 0, Some(0.2)).unwrap(); // 2
        tree.add_child(Node::new_named("E"), 0, Some(0.5)).unwrap(); // 3

        assert_eq!(tree.get_leaves(), vec![1, 2, 3]);

        tree.add_child(Node::new_named("C"), 3, Some(0.3)).unwrap(); // 4
        tree.add_child(Node::new_named("D"), 3, Some(0.4)).unwrap(); // 5

        assert_eq!(tree.get_leaves(), vec![1, 2, 4, 5]);
    }

    #[test]
    fn test_binary() {
        let mut tree = Tree::new();
        tree.add(Node::new_named("root"));

        tree.add_child(Node::new_named("0L"), 0, None).unwrap(); //1
        tree.add_child(Node::new_named("0R"), 0, None).unwrap(); //2

        assert!(tree.is_binary().unwrap());

        tree.add_child(Node::new_named("1L"), 1, None).unwrap(); //3
        tree.add_child(Node::new_named("1R"), 1, None).unwrap(); //4

        assert!(tree.is_binary().unwrap());

        tree.add_child(Node::new_named("3L"), 3, None).unwrap(); //5
        tree.add_child(Node::new_named("3R"), 3, None).unwrap(); //6
        assert!(tree.is_binary().unwrap());

        tree.add_child(Node::new_named("3?"), 3, None).unwrap(); //7
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
        tree.prune(&4).unwrap(); // prune D subtree

        assert_eq!(tree.to_newick().unwrap(), "((A)B,((H)I)G)F;");
    }

    #[test]
    fn path_from_root() {
        let tree = build_simple_tree().unwrap();
        let values: Vec<_> = get_values(&(tree.get_path_from_root(&7).unwrap()), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["F", "B", "D", "E"])
    }

    #[test]
    fn last_common_ancestor() {
        let test_cases = vec![
            ((3, 7), 1), // (A,E) -> B
            ((6, 8), 0), // (C,H) -> F
            ((3, 3), 3), // (A,A) -> A
            ((8, 5), 5), // (H,I) -> I
            ((4, 7), 4), // (D,E) -> D
        ];
        let tree = build_simple_tree().unwrap();
        for ((source, target), ancestor) in test_cases {
            println!(
                "Testing: ({:?}, {:?}) -> {:?}",
                tree.get(&source).unwrap().name,
                tree.get(&target).unwrap().name,
                tree.get(&ancestor).unwrap().name
            );
            assert_eq!(
                ancestor,
                tree.get_common_ancestor(&source, &target).unwrap()
            );
        }
    }

    #[test]
    fn get_distances_lengths() {
        let test_cases = vec![
            ((1, 3), (0.6, 2)), // (A,E)
            ((1, 4), (0.9, 3)), // (A,C)
            ((4, 5), (0.7, 2)), // (C,D)
            ((5, 2), (1.1, 3)), // (D,B)
            ((2, 5), (1.1, 3)), // (B,D)
            ((0, 2), (0.2, 1)), // (F,B)
            ((1, 1), (0.0, 0)), // (A,A)
        ];
        let tree = build_tree_with_lengths().unwrap();
        let tree2 = build_tree_without_lengths().unwrap();

        for ((idx_s, idx_t), (dist, branches)) in test_cases {
            let (d_pred, b_pred) = tree.get_distance(&idx_s, &idx_t).unwrap();
            let (d_pred2, b_pred2) = tree2.get_distance(&idx_s, &idx_t).unwrap();
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
        let values: Vec<_> = get_values(&(tree.get_leaves()), &tree)
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
                tree.to_formatted_newick(format).unwrap(),
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
            "(\"hungarian dog\":20,(\"indian elephant\":30,\"swedish horse\":60):20):50;",
            "(\"hungarian dog\":20[Comment_1],(\"indian elephant\":30,\"swedish horse\":60[Another interesting comment]):20):50;",
        ];
        for newick in newick_strings {
            let tree = Tree::from_newick(newick).unwrap();
            assert_eq!(newick, tree.to_newick().unwrap());
        }
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
        let newick_strings = vec![
            ("((D,E)B,(F,G,C)A;", NewickParseError::UnclosedBracket),
            ("((D,E)B,(F,G)C)A", NewickParseError::NoClosingSemicolon),
        ];
        for (newick, _error) in newick_strings {
            let tree = Tree::from_newick(newick);
            assert!(tree.is_err());
        }
    }

    #[test]
    fn test_subtree_leaves() {
        let test_cases = vec![
            (
                "((T0,T1)I1,(T2,T3)I2,((T4,T5)I4,(T6,T7)I4)I3)I0;",
                7,
                vec!["T4", "T5", "T6", "T7"],
            ),
            (
                "((T0,T1)I1,(T2,T3)I2,((T4,T5)I4,(T6,T7)I4)I3)I0;",
                1,
                vec!["T0", "T1"],
            ),
            (
                "((((((((T9,T8)I7,T7)I6,T6)I5,T5)I4,T4)I3,T3)I2,T2)I1,T0,T1)I0;",
                5,
                vec!["T9", "T8", "T7", "T6"],
            ),
        ];

        for (newick, root, leaves) in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            let sub_leaves: Vec<String> =
                get_values(&tree.get_subtree_leaves(&root).unwrap(), &tree)
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
            let cherries = tree.cherries();
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
            let left = tree.get_subtree_leaves(&node.children[0]).unwrap();
            let right = tree.get_subtree_leaves(&node.children[1]).unwrap();
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
        let descendants_b: Vec<_> = get_values(&tree.get_descendants(&1).unwrap(), &tree)
            .into_iter()
            .flatten()
            .sorted()
            .collect();
        let descendants_g: Vec<_> = get_values(&tree.get_descendants(&2).unwrap(), &tree)
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
        tree.add_child(Node::new_named("tip_A"), 0, Some(1.0))
            .unwrap();
        tree.add_child(Node::new_named("in_B"), 0, Some(1.0))
            .unwrap();
        tree.add_child(Node::new_named("in_C"), 2, Some(1.0))
            .unwrap();
        tree.add_child(Node::new_named("tip_D"), 3, Some(1.0))
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

            tree.init_leaf_index().unwrap();
            rota.init_leaf_index().unwrap();

            assert_eq!(
                tree.robinson_foulds(&rota).unwrap(),
                0,
                "Ref{:#?}\nRot:{:#?}",
                tree.leaf_index,
                rota.leaf_index
            );
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

            assert_eq!(t0.robinson_foulds(&t1).unwrap(), rfs[i0][i1])
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
    fn khuner_felsenstein_treedist() {
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
                t0.khuner_felsenstein(&t1).unwrap(),
                rfs[i0][i1]
            );

            assert_eq!(t0.khuner_felsenstein(&t1).unwrap(), rfs[i0][i1])
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
                t1.robinson_foulds(&t2).unwrap(),
                t1.robinson_foulds(&t2).unwrap()
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
mod tests_ete3 {

    use super::*;
    use rand::seq::SliceRandom;

    // These newick trees come from the ETE3 test suite
    // https://github.com/etetoolkit/ete/blob/439c7efd8668408fb813e78c3febce1dc627fb43/ete3/test/datasets.py

    const NW_SIMPLE2: &str = "((D, (B,C)),A);";
    const NW_SIMPLE3: &str = "((B, (A,C)),D);";
    const NW_SIMPLE4: &str = "((A, (B,C)),D);";
    const NW_SIMPLE5: &str = "(H,(A,(B,C,D)),D,T,S,(U,Y));";
    const NW_SIMPLE6: &str = "(H,(A,(B,(C),(T))),D);";

    const NW_FULL: &str = "(Ddi0002240:1.45747,Dme0014628:1.23513,(Aga0007658:1.75256,(Cin0011239:0.72821,((Fru0004507:0.184484,((Dre0008391:0,Dre0008390:0)1:0.002729,Dre0008392:0.010931)1:0.12242)1:0.14253,((Xtr0044988:0.422481,(Gga0000982:0,Gga0000981:0)1:0.109228)1:0.035488,(Mdo0014718:0.129337,((Mms0024821:0.027982,Rno0030248:0.029287)1:0.074667,((Cfa0016700:0.031643,Bta0018700:0.047366)1:0.007962,(Ptr0000001:0.005433,((Hsa0010730:0.018311[&&NHX:flag=Red],Hsa0000001:0.003656)1:0.007173,Hsa0010711:0.00273)1:0.014995)1:0.052577)1:0.01007)1:0.035671)1:0.074417)1:0.034385)1:0.18438[&&NHX:flag=Black])1:0.219467)1:0.317782[&&NHX:flag=White]);";
    const NW_DFLT: &str = "(Ddi0002240:1.45747,Dme0014628:1.23513,(Aga0007658:1.75256,(Cin0011239:0.72821,((Fru0004507:0.184484,((Dre0008391:0,Dre0008390:0)1:0.002729,Dre0008392:0.010931)1:0.12242)1:0.14253,((Xtr0044988:0.422481,(Gga0000982:0,Gga0000981:0)1:0.109228)1:0.035488,(Mdo0014718:0.129337,((Mms0024821:0.027982,Rno0030248:0.029287)1:0.074667,((Cfa0016700:0.031643,Bta0018700:0.047366)1:0.007962,(Ptr0000001:0.005433,((Hsa0010730:0.018311,Hsa0000001:0.003656)1:0.007173,Hsa0010711:0.00273)1:0.014995)1:0.052577)1:0.01007)1:0.035671)1:0.074417)1:0.034385)1:0.18438)1:0.219467)1:0.317782);";
    const NW_TOPO: &str = "(Ddi0002240,Dme0014628,(Aga0007658,(Cin0011239,((Fru0004507,((Dre0008391,Dre0008390),Dre0008392)),((Xtr0044988,(Gga0000982,Gga0000981)),(Mdo0014718,((Mms0024821,Rno0030248),((Cfa0016700,Bta0018700),(Ptr0000001,((Hsa0010730,Hsa0000001),Hsa0010711))))))))));";
    const NW_DIST: &str = "(Ddi0002240:1.45747,Dme0014628:1.23513,(Aga0007658:1.75256,(Cin0011239:0.72821,((Fru0004507:0.184484,((Dre0008391:0,Dre0008390:0):0.002729,Dre0008392:0.010931):0.12242):0.14253,((Xtr0044988:0.422481,(Gga0000982:0,Gga0000981:0):0.109228):0.035488,(Mdo0014718:0.129337,((Mms0024821:0.027982,Rno0030248:0.029287):0.074667,((Cfa0016700:0.031643,Bta0018700:0.047366):0.007962,(Ptr0000001:0.005433,((Hsa0010730:0.018311,Hsa0000001:0.003656):0.007173,Hsa0010711:0.00273):0.014995):0.052577):0.01007):0.035671):0.074417):0.034385):0.18438):0.219467):0.317782);";
    const NW2_FULL: &str = "((((((((YGR138C:0.038472,YPR156C:0.033397)1:0.050097,YOR230W:0.220261)1:0.114267,(YAL018C:0.03583,YBR287W:0.024912,YCL075W:0.030383,YDR055w:0.044474,YOR358W:0.045875)1:0.030762)1:0.061899,(YBR006W:0.012167,YBR241C:0.045848,YCR021c:0.011703,YCR061W:0.034785,YDL024c:0.030633,YDR298C:0.011081,YER141w:0.052598,YER158c:0.088493,YGR028W:0.015282,YGR149W:0.025905,YIR038C:0.021086,YJL155C:0.011531,YLR297W:0.052342,YLR423C:0.023426,YOL083W:0.046727,YOR049C:0.02909,YPL087W:0.015005)1:0.010338)1:0.02253,(YBR052C:0.044257,YBR054W:0.022772,YBR056W:0.033843,YBR183W:0.027453,YCL040w:0.076683,YCL042W:0.186523,YDL021W:0.007351,YDR032c:0.034443,YDR342C:0.024799,YDR343C:0.077931,YER053c:0.042805,YER150w:0.028763,YGR194C:0.019262,YGR244C:0.019371,YHR092C:0.064983,YIL111W:0.085602,YIR039C:0.061777,YJL079C:0.068661,YJL164C:0.06653,YJR073C:0.00477,YKL035W:0.114557,YMR250W:0.026205,YMR297W:0.057263,YNL160W:0.02005,YOR136W:0.031351,YOR273C:0.045478,YOR347C:0.017452,Ydr021w:0.029241,YPL154C:0.026716)1:0.014156)1:0.03636,((YBL078C:0.018348,YBR072W:0.008034,YBR139W:0.030263,YBR149W:0.017972,YBR169C:0.00192,YBR204C:0.017318,YDL022w:0.025509,YDL023c:0.042978,YDL091c:0.013934,YDR001C:0.014137,YDR077W:0.044424,YDR171W:0.011951,YDR178W:0.013526,YDR231C:0.005223,YDR258C:0.019661,YDR513W:0.016566,YDR529C:0.019541,YEL024w:0.020349,YFL014W:0.019925,YFR015C:0.027107,YFR033C:0.019186,YGL006W:0.019943,YGL187C:0.01129,YGL191W:0.03019,YGR008C:0.017479,YGR019W:0.017015,YGR088W:0.01558,YGR174C:0.01118,YHL021C:0.023666,YIL087C:0.011636,YIL107C:0.005697,YKL085W:0.003164,YKL103C:0.039112,YKL148C:0.022439,YKL150W?1:0.015616,YKR016W:0.017934,YKR067W:0.024008,YLL026w:0.019043,YLL041c:0.014446,YLR178C:0.029435,YLR258W:0.052186,YLR294C:0.024676,YLR299W:0.032933,YLR304C:0.045758,YLR327C:0.007435,YLR345W:0.005768,YLR395C:0.045549,YML100W:0.04311,YML120C:0.018385,YML128C:0.016054,YMR105C:0.024179,YMR110C:0.017258,YMR133W:0.023707,YMR196W:0.036526,YNL015W:0.019308,YNL055C:0.025667,YNL100W:0.004518,YNL144C:0.004901,YOL048C:0.02664,YOL053C:0.019814,YOR215C:0.010501,YOR317W:0.018442,YOR374W:0.011387,YPL078C:0.070312,YPL230W:0.020924,YPR149W:0.010801,YDR258C:0.014695)1:0.008367,(YAL060W:0.015064,YBL038W:0.014364,YBL064C:0.028571,YBR269C:0.016747,YCL035C:0.025807,YDL004W:0.007085,YDR031w:0.003761,YDR074W:0.069217,YDR272W:0.011113,YDR272W:0.013582,YDR277C:0.032421,YDR516C:0.040793,YEL039c:0.016632,YER035w:0.0248,YER067w:0.038474,YER182w:0.011733,YFR053C:0.252371,YGL037C:0.047136,YGL121C:0.022573,YGL199C:0.011657,YGR248W:0.012191,YHR104W:0.015812,YIL113W:0.013517,YIL136W:0.039322,YIL162W:0.067884,YJL144W:0.022571,YJL151C:0.017277,YJR096W:0.019693,YKL026C:0.016814,YKL142W:0.033154,YKL151C:0.02722,YKL193C:0.017688,YLL023C:0.013887,YLR080W?1:0.010144,YLR270W:0.017368,YLR356W?1:0.007843,YML004C:0.047041,YMR030W:0.012315,YMR090W:0.048551,YMR170C:0.010852,YMR195W:0.03576,YNL045W:0.020669,YNL115C:0.048069,YNL173C:0.025198,YNL200C:0.026012,YNL274C:0.005603,YNL305C:0.034282,YNR073C:0.022308,YOL071W:0.026411,YOR052C:0.024637,YOR161C:0.021991,YOR220W:0.040045,YPL004C:0.024521)1:0.013616)1:0.006005)1:0.021741,((((YAL034C:0.018415,YBL043W:0.024425,YBL049W:0.022513,YBR046C:0.03177,YCR039c:0.006232,YCR091w:0.026004,YDL079C:0.094125,YDL085w:0.032467,YDL204w:0.022855,YDL218w:0.016136,YDL223c:0.07437,YDR003W:0.07325,YDR018c:0.053591,YDR313C:0.025321,YIL097W:0.012089,YJL141C:0.035049,YKL093W:0.014752,YLR164W?1:0.011162,YLR168C:0.01707,YMR139W:0.035376,YNL093W:0.025831)1:0.015975,((YBL048W:0.007604,YBL075C:0.014956,YBR067C:0.01658,YBR101C:0.007598,YBR147W:0.00615,YDR070c:0.010301,YGR236C:0.008648,YGR238C:0.010611,YGR243W:0.004473,YJL067W:0.00774,YKL016C:0.009956,YKL217W:0.009919,YLR216C:0.041156,YLR217W:0.018315,YLR219W:0.006671,YLR295C:0.017458,YMR107W:0.003236,YNL037C:0.011856,YNL134C:0.001122,YNL194C:0.003001,YNL202W:0.013719,YOR031W:0.014812)1:0.006651,(YBR116C:0.033659,YDL020C:0.022419,YDR043C:0.056132,YDR262w:0.008193,YEL065w:0.266647,YGR146C:0.007376,YGR256W:0.035848,YJR155W:0.008535,YKL163W:0.036797,YNL036W:0.057108,YOR027W:0.016952,YPL135W:0.014251)1:0.024637)1:0.005255)1:0.007999,(YOL109W:0.104113,YAL067C:0.023929,YBR053C:0.029405,YBR203W:0.013784,YBR214W:0.010883,YBR280C:0.01757,YDL169C:0.023945,YDL245C:0.04012,YDR148C:0.019495,YDR255C:0.031336,YDR329C:0.018032,YDR533C:0.00125,YER098w:0.018484,YER142c:0.030633,YER175c:0.004409,YGL045W:0.032391,YGL059W:0.050804,YGR043C:0.016315,YGR052W:0.044183,YGR066C:0.008994,YGR130C:0.022106,YGR201C:0.026083,YGR289C:0.01392,YHR195W:0.038417,YJL153C:0.033018,YJL161W:0.02582,YJR008W:0.027243,YJR080C:0.030356,YKR058W:0.01074,YLL001w:0.008319,YLR271W:0.02894,YMR271C:0.020724,YMR311C:0.031845,YMR322C:0.021606,YNR007C:0.0205,YOL032W:0.017985,YOL053W:0.04006,YOR178C:0.004028,YOR289W:0.017451,YOR386W:0.025749,YPL185W:0.017283,YPL186C:0.00914,YPR150W:0.025491)1:0.01111)1:0.012866,((YBL015W:0.018473,YBL099W:0.046953,YBL100C:0.046467,YBR039W:0.02728,YBR132C:0.057252,YDL181W:0.052427,YDR059C:0.028415,YDR085C:0.055567,YEL011w:0.02023,YER079w:0.013469,YFL030W:0.019404,YGL259W:0.019404,YGR070W:0.028021,YGR142W:0.002517,YGR250C:0.017097,YHL024W:0.016718,YHR051W:0.020669,YIL101C:0.040969,YJL103C:0.015318,YJR048W:0.013049,YJR121W:0.018256,YKL109W:0.013086,YKL141W:0.01828,YLR038c:0.046595,YLR149C:0.034217,YLR312C:0.025271,YML054C:0.014777,YMR031C:0.041817,YMR056C:0.024105,YMR081C:0.019086,YMR136W:0.078413,YMR181C:0.059683,YMR191W:0.007052,YNL052W:0.018461,YNR001C:0.008023,YOL117W:0.037811,YOL126C:0.058282,YOL153C:0.01988,YOR035C:0.045863,YOR065W:0.023359,YOR100C:0.021039,YPL165C:0.012048,YPL201C?1:0.009305,YPL223C:0.013347,YPR020W:0.014464,YPR098C:0.027269,YPR184W:0.035196)1:0.007269,((YAL054C:0.01077,YDR125C:0.021905,YDR216W:0.027122,YEL012w:0.014583,YER024w:0.016562,YGL062W:0.266984,YIL125W:0.016436,YJL089W:0.01055,YJL163C:0.007033,YKL187C:0.012379,YLR259C:0.059779,YLR311C:0.01144,YNL009W:0.025277,YOL084W?1:0.005872,YPL262W:0.002576,YPR030W:0.011058)1:0.010088,(YBL045C:0.015491,YBR001C:0.059987,YBR117C:0.033327,YBR298C:0.021627,YCR005c:0.008913,YDL199c:0.02288,YGL153W:0.006091,YGR110W:0.019326,YHL032C:0.008845,YHR096C:0.012855,YJL045W:0.007726,YJL137C:0.016,YJR095W:0.03513,YKR076W:0.013664,YLR174W:0.025886,YLR267W:0.009587,YML042W:0.009641,YMR114C:0.242923,YMR206W:0.019201,YOR019W:0.035474,YOR391C:0.009738)1:0.005262)1:0.013775)1:0.014916)1:0.021186)1:1.75393,((((((YCR054C:0.153934,YKL061W:0.133755,YLR109W:0.142228)1:0.198112,(YBL054W:0.050037,YBR034C:0.038518,YGR177C:0.059524,YKL078W:0.026939,YKL106W:0.046488,YKR024C:0.035324,YLR380W:0.028761,YMR310C:0.034731,YOR271C:0.006213,YPL256C:0.059639)1:0.012497)1:0.065443,(YBR247C:0.050191,YBR267W:0.034203,YCLX02C:0.011558,YCR058C:0.051217,YDL151c:0.037685,YDR384C:0.011912,YDR465C:0.012999,YFL022C:0.009532,YGR159C:0.012299,YGR280C:0.013983,YHR206W:0.043975,YIL019W:0.086266,YIL065C:0.027422,YJL194W:0.096096,YKR056W:0.007362,YLL012W:0.031657,YLR073C:0.016216,YLR214W:0.039115,YLR293C:0.015673,YMR093W:0.009965,YMR259C:0.016481,YNR038W:0.01569,YOL019W:0.028424,YOL101C:0.023597,YOR272W:0.033415,YOR294W:0.033731,YPL043W:0.024688,YPL126W:0.011288,YPR144C:0.013358)1:0.010101)1:0.01987,((YBL068W:0.038979,YBR158W:0.008382,YBR186W:0.116086,YCR055C:0.061205,YCR065w:0.003695,YDL037c:0.048015,YDL042C:0.022897,YDL182w:0.039913,YDR508C:0.030289,YER145c:0.024536,YER165w:0.005328,YGR245C:0.009444,YJL069C:0.017205,YLR397C:0.033973,YMR108W:0.024272,YMR309C:0.017863,YNL124W:0.016342,YNL189W:0.033389,YOL007C:0.062976,YOL130W?1:0.068678)1:0.017449,((YAL043C:0.085408,YDR006C:0.131879,YHR205W:0.03885,YJL195C:0.040006,YPL061W:0.070549)1:0.031474,(YDR101C:0.068051,YKL185W:0.029377,YNL061W:0.047111,YOR091W:0.057478)1:0.047015)1:0.102596)1:0.046055)1:0.039969,(YAL003W:0.027682,YAL025C:0.008438,YAL046C:0.01662,YAL059W:0.05127,YAR071W:0.027037,YAR074C:0.005698,YBL024W:0.00817,YBL039C:0.039903,YBR069C:0.243928,YBR079C:0.025916,YBR089W:0.031468,YBR093C:0.043355,YBR142W:0.05141,YBR238C:0.05107,YBR266C:0.043018,YCL046W:0.019256,YCL054W:0.041245,YCR053w:0.025386,YDL050c:0.008063,YDL051W:0.022365,YDL052C:0.053016,YDL063c:0.018176,YDL084w:0.016863,YDL111c:0.008781,YDL122W:0.043325,YDL131w:0.108001,YDL140C:0.024567,YDL148c:0.046571,YDL152w:0.055074,YDL153c:0.017418,YDL167C:0.009333,YDL213c:0.017884,YDR024w:0.148026,YDR144C:0.010083,YDR165W:0.015923,YDR206W:0.005912,YDR324C:0.042943,YDR341C:0.015004,YDR361C:0.015177,YDR398W:0.018739,YDR492W:0.019353,YDR502C:0.023435,YDR527W:0.025343,YEL033w:0.020575,YEL040w:0.048325,YEL046c:0.015069,YER002w:0.014535,YER006w:0.017578,YER025w:0.014587,YER036c:0.016587,YER043c:0.03011,YER049w:0.045735,YER052c:0.022026,YER070w:0.021812,YER110c:0.022076,YGR103W:0.034426,YGR123C:0.04964,YGR124W:0.006737,YGR145W:0.02528,YGR155W:0.054731,YGR160W:0.013018,YGR162W:0.020365,YGR264C:0.018057,YHR046C:0.180512,YHR052W:0.010195,YHR064C:0.016477,YHR070W:0.03684,YHR215W:0.030425,YIL053W:0.012927,YIL066C:0.008671,YIL091C:0.003393,YJL050W:0.016488,YJL109C:0.060257,YJL122W:0.017642,YJR003C:0.022002,YJR016C:0.014266,YJR041C:0.020306,YJR054W:0.037849,YKL076C:0.038156,YKL082C:0.012586,YKL181W:0.005029,YKL191W:0.099536,YLL008w:0.016875,YLL021w:0.024558,YLR068W:0.079458,YLR084C:0.035147,YLR180W:0.009062,YLR196W:0.046001,YLR221C:0.028516,YLR223C:0.057444,YLR243W:0.010906,YLR244C:0.037764,YLR280C:0.020274,YLR355C:0.025286,YLR409C:0.009258,YLR435W:0.092123,YLR449W:0.01871,YML123C:0.018546,YMR011W:0.032472,YMR014W:0.007539,YMR037C:0.063027,YMR049C:0.011294,YMR058W:0.010664,YMR080C:0.040422,YMR128W:0.033444,YMR129W:0.046421,YMR146C:0.021802,YMR217W:0.041486,YMR229C:0.004043,YMR239C:0.006788,YMR241W:0.012776,YMR246W:0.023201,YMR290C:0.003882,YNL002C:0.029441,YNL013C:0.037125,YNL060C:0.001669,YNL062C:0.014338,YNL075W:0.012303,YNL110C:0.007903,YNL112W:0.034459,YNL120C:0.018624,YNL123W:0.022098,YNL132W:0.00337,YNL141W:0.009783,YNL174W:0.101658,YNL182C:0.007662,YNL207W:0.014884,YNL216W:0.170324,YNL221C:0.015522,YNL256W:0.011304,YNL308C:0.045905,YNL313C:0.026778,YNL327W:0.008828,YNR025C:0.048977,YNR043W:0.020504,YNR046W:0.055271,YNR054C:0.00537,YNR075W:0.01383,YOL021C:0.010752,YOL080C:0.031182,YOR116C:0.029862,YOR146W:0.05367,YOR206W:0.051387,YOR233W:0.03687,YOR335C:0.009945,YOR341W:0.023363,YOR355W:0.015312,YOR361C:0.022609,YPL012W:0.003133,YPL019C:0.014463,YPL030W:0.011898,YPL032C:0.047324,YPL093W:0.009715,YPL183C:0.087818,YPL211W:0.00693,YPL226W:0.013388,YPR009W:0.025512,YPR034W:0.010374,YPR110C:0.00993,YPR136C:0.133899,YPR143W:0.021514,YPR145W:0.025839)1:0.00584)1:0.024107,((YHR007C:0.024321,YAL023C:0.037896,YAL036C:0.027457,YAR073W:0.007629,YBL076C:0.007966,YBL087C:0.016367,YBR048W:0.03756,YBR092C:0.02723,YBR121C:0.027699,YBR143C:0.022332,YBR156C:0.004017,YBR249C:0.000577,YCL053C:0.012758,YDL014W:0.011777,YDL060w:0.008162,YDL082w:0.006478,YDL145C:0.02523,YDL208W:0.003093,YDR023W:0.039428,YDR037W:0.023439,YDR060w:0.005655,YDR062W:0.01396,YDR064W:0.002168,YDR321W:0.006361,YDR365C:0.006326,YDR399W:0.003641,YDR447C:0.015148,YDR449C:0.001937,YDR450W:0.016031,YDR471W:0.015488,YEL026w:0.019446,YEL054c:0.008555,YER056c:0.022337,YER060w:0.036023,YER129w:0.011346,YFL045C:0.053502,YGL008C:0.012515,YGL029W:0.067475,YGL078C:0.045179,YGL092W:0.045425,YGL103W:0.00955,YGL225W:0.027201,YGR034W:0.012058,YGR060W:0.013105,YGR061C:0.009325,YGR214W:0.007962,YHL033C:0.014192,YHR019C:0.019476,YHR049W:0.036925,YHR089C:0.010167,YHR128W:0.008165,YHR203C:0.0076,YHR216W:0.002384,YIL052C:0.013691,YIL069C:0.012335,YIL133C:0.004129,YJL111W:0.017878,YJL138C:0.009234,YJL177W:0.002417,YJL183W:0.005059,YJR071W:0.008702,YJR105W:0.058597,YJR123W:0.003913,YKL009W:0.009543,YKL057C:0.034921,YKL081W:0.007191,YKL156W:0.011042,YKR013W:0.029496,YKR025W:0.021172,YKR094C:0.008062,YLL004w:0.026716,YLL044W:0.008058,YLR009W:0.012315,YLR048w:0.004366,YLR056w:0.028151,YLR061W:0.003944,YLR083c:0.016167,YLR129w:0.052379,YLR134w:0.018072,YLR175W:0.025815,YLR179C:0.014679,YLR185W:0.007868,YLR186W:0.055986,YLR212C:0.006866,YLR339C:0.013427,YLR384C:0.043396,YLR413W:0.007107,YLR432W:0.00039,YLR448W:0.013913,YML059C:0.068679,YMR199W:0.180049,YMR202W:0.01383,YMR205C:0.01129,YMR242C:0.002943,YMR308C:0.012121,YMR318C:0.009922,YMR321C:0.007602,YNL065W:0.04082,YNL087W:0.011444,YNL096C:0.016268,YNL111C:0.009092,YNL175C:0.021305,YNL235C:0.017678,YNL247W:0.010779,YNL301C:0.006244,YNL303W:0.026882,YNR009W:0.069305,YNR051C:0.005803,YNR053C:0.011495,YNR067C:0.038533,YOL010W:0.013445,YOL014W:0.050775,YOL120C:0.01174,YOL061W:0.014505,YOL077C:0.027831,YOR063W:0.004912,YOR096W:0.009238,YOR153W:0.010446,YOR234C:0.008998,YOR277C:0.023939,YOR293W:0.013614,YOR309C:0.007759,YPL034W:0.035088,YPL143W:0.030537,YPL160W:0.018401,YPL198W:0.002841,YPL243W:0.014395,YPL245W:0.004341,YPL266W:0.01551,YPR033C:0.028234,YPR044C:0.016689,YPR074C:0.006226,YPR112C:0.012512,YPR125W:0.03021,YHR098C:0.005812)1:0.005387,((YAL012W:0.037345,YAL038W:0.007617,YBL027W:0.013554,YBR106W:0.006914,YBR189W:0.010072,YBR191W:0.00515,YDL075W:0.007897,YDL083C:0.031836,YDL136w:0.006651,YDL191W:0.011001,YDL210W:0.01511,YDR366C:0.015185,YDR417C:0.003784,YDR418W:0.001455,YDR500C:0.007442,YER074w:0.008117,YER117w:0.008649,YER131w:0.002672,YFR031BC:0.012496,YGL030W:0.003864,YGL076C:0.013045,YGL123W:0.012111,YGL135W:0.018213,YGL147C:0.004678,YGR118W:0.013782,YGR148C:0.005154,YGR285C:0.024308,YHL001W:0.021713,YHL015W:0.009606,YHR010W:0.002032,YHR141C:0.007074,YHR208W:0.038459,YIL018W:0.000289,YIL039W:0.015756,YIL148W:0.016059,YJL080C:0.017726,YJL136C:0.009858,YJL148W:0.025592,YJL188C?1:0.009377,YJL189W:0.030581,YJR063W:0.036817,YJR145C:0.011576,YKL006W:0.005959,YKR043C:0.016532,YKR057W:0.009924,YKR059W:0.011284,YLL045c:0.014089,YLL047W:0.023354,YLR029c:0.010768,YLR044c:0.011219,YLR060w:0.020979,YLR062C:0.006945,YLR076C:0.010549,YLR198C:0.003174,YLR249W:0.022251,YLR325C:0.009715,YLR340W:0.008365,YLR344W:0.010317,YLR367W:0.003094,YLR388W:0.021736,YML063W:0.006051,YMR121C:0.017411,YMR131C:0.026306,YNL119W:0.014916,YNL162W:0.016264,YNL302C:0.004123,YNR050C:0.058206,YOL121C:0.009887,YOL127W:0.013974,YOL040C:0.007057,YOR224C:0.013626,YOR254C:0.017923,YOR310C:0.006634,YOR312C:0.00351,YPL079W:0.030247,YPL131W:0.022146,YPL142C:0.015716,YPL220W:0.009552,YPL244C:0.022966,YPR102C:0.008204,YPR137W:0.017679)1:0.002907,((YBR032W:0.016097,YBR181C:0.009522,YBR187W:0.018681,YDR012W:0.011279,YGL102C:0.009205,YGR027C:0.00369,YJL157C:0.024018,YJL190C:0.007303,YLR058c:0.073551,YLR264W:0.018569,YLR300W:0.020115,YLR342W:0.034792,YLR359W:0.022991,YLR372W:0.011058,YMR300C:0.0321,YMR305C:0.014402,YNL069C:0.028621,YOR133W:0.004537,YOR167C:0.014198,YOR182C:0.023853,YOR326W:0.024632,YPL080C:0.064427,YPL081W:0.012279,YPL090C:0.012152)1:0.006473,(YOR298W:0.233587,((YDR019C:0.040149,YKL096W:0.037806)1:0.064948,(YKL108W:0.033111,YPR113W:0.038607)1:0.063282)1:0.016535)1:0.229462)1:0.024472)1:0.008239)1:0.017711)1:0.082401);";

    #[test]
    fn parse_full_nhx() {
        let tree = Tree::from_newick(NW_FULL).expect("Could not parse NW_FULL");
        assert_eq!(NW_FULL, tree.to_newick().expect("Could not write NW_FULL"));
    }

    #[test]
    fn parse_complex_newick() {
        let tree = Tree::from_newick(NW2_FULL).expect("Could not parse NW2_FULL");
        assert_eq!(
            NW2_FULL,
            tree.to_newick().expect("Could not write NW2_FULL")
        );
    }

    #[test]
    fn parse_weird_topologies() {
        let cases = [(NW_SIMPLE5, "NW_SIMPLE5"), (NW_SIMPLE6, "NW_SIMPLE6")];
        for (newick, name) in cases {
            let tree =
                Tree::from_newick(newick).unwrap_or_else(|_| panic!("Could not parse {name}"));
            assert_eq!(
                newick,
                tree.to_newick()
                    .unwrap_or_else(|_| panic!("Could not write {name}"))
            )
        }
    }

    #[ignore]
    #[test]
    fn parse_single_node_trees() {
        for newick in ["(hola);", "hola;"] {
            let tree = Tree::from_newick(newick).expect("Could not read tree");
            assert_eq!(newick, tree.to_newick().expect("Could not write tree"));
            println!("Parsed '{newick}' succesfully");
        }
    }

    #[test]
    fn parse_with_features() {
        let newick = "(((A[&&NHX:name=A],B[&&NHX:name=B])[&&NHX:name=NoName],C[&&NHX:name=C])[&&NHX:name=I],(D[&&NHX:name=D],F[&&NHX:name=F])[&&NHX:name=J])[&&NHX:name=root];";
        let tree = Tree::from_newick(newick).expect("Could not parse tree");

        assert_eq!(newick, tree.to_newick().expect("Could not write tree"));
    }

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

        assert_eq!(tree.get_common_ancestor(&a.id, &i.id).unwrap(), i.id);
        assert_eq!(tree.get_common_ancestor(&a.id, &d.id).unwrap(), root.id);

        assert_eq!(tree.get_distance(&a.id, &i.id).unwrap().0.unwrap(), 0.101);
        assert_eq!(tree.get_distance(&a.id, &b.id).unwrap().0.unwrap(), 0.11);
        assert_eq!(tree.get_distance(&a.id, &a.id).unwrap().0.unwrap(), 0.0);
        assert_eq!(tree.get_distance(&i.id, &i.id).unwrap().0.unwrap(), 0.0);
        assert_eq!(
            tree.get_distance(&a.id, &root.id).unwrap(),
            tree.get_distance(&root.id, &a.id).unwrap()
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

        assert_eq!(tree.get_distance(&root.id, &a.id).unwrap().1, 3);
        assert_eq!(tree.get_distance(&root.id, &c.id).unwrap().1, 2);
        assert_eq!(tree.get_distance(&root.id, &root.id).unwrap().1, 0);
        assert_eq!(tree.get_distance(&d.id, &d.id).unwrap().1, 0);
    }

    #[test]
    fn traversals() {
        let tree = Tree::from_newick("((3,4)2,(6,7)5)1;").unwrap();
        let root = tree.get_root().unwrap();

        let postorder = "3426751";
        let preorder = "1234567";
        let levelorder = "1253467";

        fn get_str(iter: &[usize], tree: &Tree) -> String {
            iter.iter()
                .map(|id| tree.get(id).unwrap().name.clone().unwrap())
                .collect()
        }

        assert_eq!(get_str(&tree.postorder(&root).unwrap(), &tree), postorder);
        assert_eq!(get_str(&tree.preorder(&root).unwrap(), &tree), preorder);
        assert_eq!(get_str(&tree.levelorder(&root).unwrap(), &tree), levelorder);
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

            let rf = t1.robinson_foulds(&t2).unwrap();

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
