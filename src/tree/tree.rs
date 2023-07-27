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
use super::{Edge, NodeId};

use crate::distance::{DistanceMatrix, MatrixError};

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

/// A Phylogenetic tree
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    leaf_index: RefCell<Option<Vec<String>>>,
    partitions: RefCell<Option<HashMap<FixedBitSet, Option<Edge>>>>,
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
        edge: Option<Edge>,
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

    /// Get a reference to a node in the tree by name
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
    pub fn height(&self) -> Result<Edge, TreeError> {
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
            .map_or(Err(TreeError::IsEmpty), Ok)
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
    pub fn diameter(&self) -> Result<Edge, TreeError> {
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
            .map_or(Err(TreeError::IsEmpty), Ok)
    }

    /// Returns the lenght of the trees
    /// (i.e. the sum of branch lenghts)
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// assert_eq!(tree.length().unwrap(), 1.5);
    /// ```
    pub fn length(&self) -> Result<Edge, TreeError> {
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
    fn get_partition(&self, index: &NodeId) -> Result<FixedBitSet, TreeError> {
        self.init_leaf_index()?;

        let subtree_leaves = self.get_subtree_leaves(index)?;
        let indices = subtree_leaves
            .iter()
            .filter_map(|index| self.get(index).as_ref().unwrap().name.clone())
            .map(|name| {
                let v = self.leaf_index.borrow().clone();
                v.map(|v| v.iter().position(|n| *n == name).unwrap())
                    .unwrap()
            });

        let mut bitset = FixedBitSet::with_capacity(self.n_leaves());
        for index in indices {
            bitset.insert(index);
        }

        let mut toggled = bitset.clone();
        toggled.toggle_range(..);

        Ok(toggled.min(bitset))
    }

    /// Caches partitions for distance computation
    fn init_partitions(&self) -> Result<(), TreeError> {
        self.init_leaf_index()?;

        if self.partitions.borrow().is_some() {
            return Ok(());
        }

        let mut partitions: HashMap<FixedBitSet, Option<f64>> = HashMap::new();

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
            let old_len = partitions.get(&part);

            let len = match (new_len, old_len) {
                (None, None) => None,
                (Some(new_len), Some(old_len)) => old_len.map(|v| v + new_len),
                (Some(new_len), None) => Some(new_len),
                (None, Some(old_len)) => *old_len,
            };

            partitions.insert(part, len);
        }

        (*self.partitions.borrow_mut()) = Some(partitions);

        Ok(())
    }

    /// Get all partitions of a tree
    pub fn get_partitions(&self) -> Result<HashSet<FixedBitSet>, TreeError> {
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

    /// Get all partitions of a tree along with corresponding branch lengths
    pub(crate) fn get_partitions_with_lengths(
        &self,
    ) -> Result<HashMap<FixedBitSet, f64>, TreeError> {
        self.init_leaf_index()?;
        self.init_partitions()?;

        let mut partitions = HashMap::new();
        for (bitset, len) in self.partitions.borrow().as_ref().unwrap().iter() {
            if len.is_none() {
                return Err(TreeError::MissingBranchLengths);
            }
            partitions.insert(bitset.clone(), len.unwrap());
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

        for (edge, len_s) in partitions_s.iter() {
            if let Some(len_o) = partitions_o.get(edge) {
                dist += (len_s - len_o).abs()
            } else {
                dist += len_s
            }
        }

        for (edge, len_o) in partitions_o.iter() {
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

        for (edge, len_s) in partitions_s.iter() {
            if let Some(len_o) = partitions_o.get(edge) {
                dist += f64::powi(len_s - len_o, 2)
            } else {
                dist += f64::powi(*len_s, 2)
            }
        }

        for (edge, len_o) in partitions_o.iter() {
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

        for (edge, len_s) in partitions_s.iter() {
            if let Some(len_o) = partitions_o.get(edge) {
                rf_weight += (len_s - len_o).abs();
                kf += f64::powi(len_s - len_o, 2);
                intersection += 1.0;
            } else {
                rf_weight += len_s;
                kf += f64::powi(*len_s, 2);
            }
        }

        for (edge, len_o) in partitions_o.iter() {
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
            return Ok((None, 0));
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

        for list in vec![root_to_source, root_to_target] {
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
    pub fn distance_matrix_recursive(&self) -> Result<DistanceMatrix<Edge>, TreeError> {
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
        let mut pairwise_vec = vec![NaiveSum::zero(); n * n / 2];

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
                        let vec_idx = ((2 * n - 3 - i) * i) / 2 + j - 1;

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

        Ok(matrix)
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
            .cloned()
            .filter(|node| !node.deleted && node.parent.is_some() && node.children.len() == 1)
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
    ///
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

        for node_id in to_binarize.iter() {
            loop {
                let mut children = self.get(node_id)?.children.clone();
                children.shuffle(rng);

                let parent = self.add_child(Node::new(), *node_id, Some(0.0))?;

                for _ in 0..2 {
                    let child = children.pop().unwrap();
                    let edge = self.get(&child)?.parent_edge.clone();
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
    fn to_newick_impl(&self, root: &NodeId) -> Result<String, TreeError> {
        let root = self.get(root)?;
        if root.children.is_empty() {
            Ok(root.to_newick())
        } else {
            Ok("(".to_string()
                + &(root
                    .children
                    .iter()
                    .map(|child_idx| self.to_newick_impl(child_idx).unwrap()))
                .collect::<Vec<String>>()
                .join(",")
                + ")"
                + &(root.to_newick()))
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
    /// dbg!(&tree);
    ///
    /// assert_eq!(tree.to_newick().unwrap(), newick);
    /// ```
    pub fn to_newick(&self) -> Result<String, TreeError> {
        let root = self.get_root()?;
        Ok(self.to_newick_impl(&root)? + ";")
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

/// Methods to infer [Tree] ojects from [DistanceMatrix] objects
///
/// ----
/// ----
// impl Tree {
// /// #########################
// /// # INFER TREES FROM DATA #
// /// #########################
//
// fn neighbour_joining(matrix: DistanceMatrix<f64>) -> Self {
//     let mut tree = Tree::new();
//     let mut matrix = matrix;
//     let mut q = DistanceMatrix::new(matrix.size);
//
//     for _ in 0..(matrix.size - 3) {
//         for pair in matrix.ids.iter().enumerate().combinations(2) {
//
//         }
//     }
//
//     tree
// }
// }

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
// #[allow(clippy::excessive_precision)]
mod tests {

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
            ((1, 3), (Some(0.6), 2)), // (A,E)
            ((1, 4), (Some(0.9), 3)), // (A,C)
            ((4, 5), (Some(0.7), 2)), // (C,D)
            ((5, 2), (Some(1.1), 3)), // (D,B)
            ((2, 5), (Some(1.1), 3)), // (B,D)
            ((0, 2), (Some(0.2), 1)), // (F,B)
            ((1, 1), (None, 0)),      // (A,A)
        ];
        let tree = build_tree_with_lengths().unwrap();

        for ((idx_s, idx_t), (dist, branches)) in test_cases {
            let (d_pred, b_pred) = tree.get_distance(&idx_s, &idx_t).unwrap();
            assert_eq!(branches, b_pred);
            match dist {
                None => assert!(d_pred.is_none()),
                Some(d) => {
                    assert!(d_pred.is_some());
                    assert!((d_pred.unwrap() - d).abs() < f64::EPSILON);
                }
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
            "(\"hungarian dog\"[Comment_1]:20,(\"indian elephant\":30,\"swedish horse[Another interesting comment]\":60):20):50;",
        ];
        for newick in newick_strings {
            let tree = Tree::from_newick(newick).unwrap();
            assert_eq!(newick, tree.to_newick().unwrap());
        }
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
        let trees = vec![
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
        let rfs = vec![
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
        let trees = vec![
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
        let rfs = vec![
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
        let trees = vec![
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
        let rfs: Vec<Vec<f64>> = vec![
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
            vec![
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
        let trees = vec![
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
