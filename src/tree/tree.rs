use fixedbitset::FixedBitSet;
use std::iter::zip;
use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    fs,
    path::Path,
};

use thiserror::Error;

use super::node::Node;
use super::{Edge, NodeId};

#[derive(Error, Debug)]
pub enum TreeError {
    #[error("This tree is not Binary.")]
    IsNotBinary,
    #[error("This tree is not rooted.")]
    IsNotRooted,
    #[error("This tree is empty.")]
    IsEmpty,
    #[error("All your leaf nodes must ne named.")]
    UnnamedLeaves,
    #[error("Your leaf names must be unique.")]
    DuplicateLeafNames,
    #[error("The leaf index of the tree is not initialized.")]
    LeafIndexNotInitialized,
    #[error("The tree must have all branch lengths.")]
    MissingBranchLengths,
    #[error("The trees have different tips indices.")]
    DifferentTipIndices,
    #[error("There is no node with index: {0}")]
    NodeNotFound(NodeId),
    #[error("No root node found")]
    RootNotFound,
    #[error("Error writing tree to file")]
    IoError(#[from] std::io::Error),
}

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("Cannot have whitespace in number field.")]
    WhiteSpaceInNumber,
    #[error("Missing a closing bracket.")]
    UnclosedBracket,
    #[error("The tree is missin a semi colon at the end.")]
    NoClosingSemicolon,
    #[error("Problem with building the tree.")]
    TreeError(#[from] TreeError),
    #[error("Could not parse a branch length")]
    FloatError(#[from] std::num::ParseFloatError),
    #[error("Parent node of subtree not found")]
    NoSubtreeParent,
    #[error("Problem reading file")]
    IoError(#[from] std::io::Error),
}

/// A Vector backed Tree structure
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    leaf_index: RefCell<Option<Vec<String>>>,
    partitions: RefCell<Option<HashMap<FixedBitSet, Option<Edge>>>>,
}

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
    // # adding and getting nodes #
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
    /// assert_eq!(tree.get(&root_id).children.len(), 2);
    ///
    /// // The depths of child nodes are derived from the parent node
    /// assert_eq!(tree.get(&left).depth, 1);
    /// assert_eq!(tree.get(&right).depth, 1);
    ///
    /// // If an edge length is specified then it is set in both child and parent
    /// assert_eq!(tree.get(&right).parent_edge, Some(0.1));
    /// assert_eq!(tree.get(&root_id).get_child_edge(&right), Some(0.1));
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
        node.set_depth(self.get(&parent).depth + 1);

        let id = self.add(node);

        self.get_mut(&id).set_id(id);
        self.get_mut(&parent).add_child(id, edge);

        Ok(id)
    }

    /// Get a reference to a specific Node of the tree
    pub fn get(&self, id: &NodeId) -> &Node {
        &self.nodes[*id]
    }

    /// Get a mutable reference to a specific Node of the tree
    pub fn get_mut(&mut self, id: &NodeId) -> &mut Node {
        &mut self.nodes[*id]
    }

    /// Get a reference to a node in the tree by name
    /// ```
    /// use phylotree::tree::{Tree, Node};
    /// 
    /// let mut tree = Tree::new();
    /// let root_idx = tree.add(Node::new_named("root"));
    /// let child_idx = tree.add_child(Node::new_named("child"), root_idx, None).unwrap();
    /// 
    /// assert_eq!(tree.get_by_name("child"), Some(tree.get(&child_idx)));
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
            .filter(|&node| node.is_tip())
            .map(|node| node.id)
            .collect()
    }

    /// Gets the leaf indices in the subtree rooted at node of a specified index
    /// ```
    /// use phylotree::tree::Tree;
    /// 
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    /// let sub_root = tree.get_by_name("E").unwrap();
    /// let sub_leaves: Vec<_> = tree.get_subtree_leaves(&sub_root.id)
    ///     .iter()
    ///     .map(|id| tree.get(id).name.clone().unwrap())
    ///     .collect();
    /// 
    /// assert_eq!(
    ///     sub_leaves,
    ///     vec![String::from("C"), String::from("D")]
    /// )
    /// ```
    pub fn get_subtree_leaves(&self, index: &NodeId) -> Vec<NodeId> {
        let mut indices = vec![];
        if self.get(index).is_tip() {
            return vec![*index];
        }

        for &child_idx in self.get(index).children.iter() {
            indices.extend(self.get_subtree_leaves(&child_idx))
        }

        indices
    }

    // #######################################
    // # getting characteristics of the tree #
    // #######################################

    /// Check if the tree is Binary
    pub fn is_binary(&self) -> bool {
        for node in self.nodes.iter() {
            // The "root" node of an unrooted binary tree has 3 children
            if node.parent.is_none() && node.children.len() > 3 {
                return false;
            }
            if node.children.len() > 2 {
                return false;
            }
        }
        true
    }

    /// Checks if the tree is rooted (i.e. the root node exists and has exactly 2 children)
    pub fn is_rooted(&self) -> Result<bool, TreeError> {
        let root_id = self.get_root()?;

        Ok(!self.nodes.is_empty() && self.get(&root_id).children.len() == 2)
    }

    /// Returns the number of nodes in the tree
    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the number of leaves in the tree
    pub fn n_leaves(&self) -> usize {
        self.nodes.iter().filter(|&node| node.is_tip()).count()
    }

    // ##########################
    // # Find paths in the tree #
    // ##########################

    /// Returns the path from the node to the root
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let path: Vec<_> = tree.get_path_from_root(&5)
    ///     .iter()
    ///     .map(|id| tree.get(id).name.clone().unwrap())
    ///     .collect();
    ///
    /// assert_eq!(
    ///     path, 
    ///     vec![String::from("F"), String::from("B"), String::from("D"), String::from("E")]
    /// )
    /// ```
    pub fn get_path_from_root(&self, node: &NodeId) -> Vec<NodeId> {
        let mut path = vec![];
        let mut current_node = *node;
        loop {
            path.push(current_node);
            match self.get(&current_node).parent {
                Some(parent) => current_node = parent,
                None => break,
            }
        }

        path.into_iter().rev().collect()
    }

    /// Gets the most recent common ancestor between two tree nodes
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
    /// let ancestor = tree.get_common_ancestor(
    ///     &tree.get_by_name("A").unwrap().id,
    ///     &tree.get_by_name("D").unwrap().id,
    /// );
    /// 
    /// assert_eq!(tree.get(&ancestor).name, Some("B".to_owned()))
    /// ```
    pub fn get_common_ancestor(&self, source: &NodeId, target: &NodeId) -> usize {
        if source == target {
            return *source;
        }
        let root_to_source = self.get_path_from_root(source);
        let root_to_target = self.get_path_from_root(target);

        let cursor = zip(root_to_source.iter(), root_to_target.iter())
            .enumerate()
            .filter(|(_, (s, t))| s != t)
            .map(|(idx, _)| idx)
            .next()
            .unwrap_or_else(|| {
                // One node is a child of the other
                root_to_source.len().min(root_to_target.len())
            });

        root_to_source[cursor - 1]
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
    /// );
    /// 
    /// assert_eq!(num_edges, 4);
    /// assert!(sum_edge_lengths.is_none());
    /// ```
    pub fn get_distance(&self, source: &NodeId, target: &NodeId) -> (Option<f64>, usize) {
        let mut dist = 0.0;
        let mut branches = 0;
        let mut all_dists = true;

        if source == target {
            return (None, 0);
        }

        let root_to_source = self.get_path_from_root(source);
        let root_to_target = self.get_path_from_root(target);

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
                if let Some(d) = self.get(node).parent_edge {
                    dist += d;
                } else {
                    all_dists = false;
                }
                branches += 1;
            }
        }

        if all_dists {
            (Some(dist), branches)
        } else {
            (None, branches)
        }
    }

    /// Returns the height of the tree
    /// (i.e. the number of edges from the root to the deepest tip)
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
    ///
    /// assert_eq!(tree.height(), Some(2));
    /// ```
    pub fn height(&self) -> Option<usize> {
        self.nodes.iter().map(|node| node.depth).max()
    }

    /// Checks if the tree is rooted and binary
    fn check_rooted_binary(&self) -> Result<(), TreeError> {
        if !self.is_rooted()? {
            Err(TreeError::IsNotRooted)
        } else if !self.is_binary() {
            Err(TreeError::IsNotBinary)
        } else {
            Ok(())
        }
    }

    /// Computes the number of cherries in a tree
    pub fn cherries(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;
        if !self.nodes.is_empty() {
            let mut n = 0;
            for node in self.nodes.iter() {
                if node.children.len() == 2
                    && self.get(&node.children[0]).is_tip()
                    && self.get(&node.children[1]).is_tip()
                {
                    n += 1;
                }
            }
            Ok(n)
        } else {
            Err(TreeError::IsEmpty)
        }
    }

    /// Computes the Colless index for the tree
    pub fn colless(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;

        let colless = self
            .nodes
            .iter()
            .filter(|node| !node.is_tip())
            .map(|node| {
                if node.children.is_empty() {
                    return 0;
                }
                let left = self.get_subtree_leaves(&node.children[0]).len();
                let right = if node.children.len() > 1 {
                    self.get_subtree_leaves(&node.children[1]).len()
                } else {
                    0
                };

                left.abs_diff(right)
            })
            .sum();

        Ok(colless)
    }

    /// Computes the normalized colless statistic with a Yule null model
    pub fn colless_yule(&self) -> Result<f64, TreeError> {
        self.colless().map(|i_c| {
            let n = self.n_leaves() as f64;
            let e_i_c = n * n.ln() + (0.57721566 - 1. - f64::ln(2.0)) * n;

            (i_c as f64 - e_i_c) / n
        })
    }

    /// Computes the normalized colless statistic with a PDA null model
    pub fn colless_pda(&self) -> Result<f64, TreeError> {
        self.colless()
            .map(|i_c| i_c as f64 / f64::powf(self.n_leaves() as f64, 3.0 / 2.0))
    }

    /// Computes the sackin statistic for tree
    pub fn sackin(&self) -> Result<usize, TreeError> {
        self.check_rooted_binary()?;

        Ok(self
            .get_leaves()
            .iter()
            .map(|tip_idx| self.get(tip_idx).depth)
            .sum())
    }

    /// Computes the normalized sackin statistic with a Yule null model
    pub fn sackin_yule(&self) -> Result<f64, TreeError> {
        self.sackin().map(|i_n| {
            let n = self.n_leaves();
            let sum: f64 = (2..=n).map(|i| 1.0 / (i as f64)).sum();

            (i_n as f64 - 2.0 * (n as f64) * sum) / n as f64
        })
    }

    /// Computes the normalized sackin statistic with a PDA null model
    pub fn sackin_pda(&self) -> Result<f64, TreeError> {
        self.sackin()
            .map(|i_n| i_n as f64 / f64::powf(self.n_leaves() as f64, 3.0 / 2.0))
    }

    // ########################
    // # read and write trees #
    // ########################

    /// Generate newick representation of tree
    fn to_newick_impl(&self, root: &NodeId) -> String {
        let root = self.get(root);
        if root.children.is_empty() {
            root.to_newick()
        } else {
            "(".to_string()
                + &(root
                    .children
                    .iter()
                    .map(|child_idx| self.to_newick_impl(child_idx)))
                .collect::<Vec<String>>()
                .join(",")
                + ")"
                + &(root.to_newick())
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
        Ok(self.to_newick_impl(&root) + ";")
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
    pub fn from_newick(newick: &str) -> Result<Self, ParseError> {
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
                        tree.get_mut(&index)
                    } else {
                        if let Some(parent) = parent_stack.last() {
                            current_index = Some(tree.add_child(Node::new(), *parent, None)?);
                        } else {
                            unreachable!("Sould not be possible to have named child with no parent")
                        };
                        tree.get_mut(current_index.as_ref().unwrap())
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
                        tree.get_mut(&index)
                    } else {
                        if let Some(parent) = parent_stack.last() {
                            current_index = Some(tree.add_child(Node::new(), *parent, None)?);
                        } else {
                            unreachable!("Sould not be possible to have named child with no parent")
                        };
                        tree.get_mut(current_index.as_ref().unwrap())
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
                        return Err(ParseError::NoSubtreeParent);
                    }
                }
                ';' => {
                    // Finish parsing the Tree
                    if !open_delimiters.is_empty() {
                        return Err(ParseError::UnclosedBracket);
                    }
                    let node = tree.get_mut(current_index.as_ref().unwrap());
                    node.name = current_name;
                    node.comment = current_comment;
                    if let Some(length) = current_length {
                        node.parent_edge = Some(length.parse()?);
                    }

                    // Finishing pass to make sure that branch lenghts are set in both children and parents
                    let ids: Vec<_> = tree.nodes.iter().map(|node| node.id).collect();
                    for node_id in ids {
                        if let Some(edge) = tree.get(&node_id).parent_edge {
                            if let Some(parent) = tree.get(&node_id).parent {
                                tree.get_mut(&parent).set_child_edge(&node_id, Some(edge));
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
                                return Err(ParseError::WhiteSpaceInNumber);
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

        Err(ParseError::NoClosingSemicolon)
    }

    /// Writes the tree to a newick file
    pub fn to_file(&self, path: &Path) -> Result<(), TreeError> {
        match fs::write(path, self.to_newick()?) {
            Ok(_) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    /// Creates a tree from a newick file
    pub fn from_file(path: &Path) -> Result<Self, ParseError> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_newick(&newick_string)
    }
}

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}
