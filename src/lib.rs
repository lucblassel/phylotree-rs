#![feature(is_some_and)]

use std::{
    collections::{HashSet, VecDeque},
    fmt::{Debug, Display},
    fs,
    iter::zip,
    path::Path,
};

use itertools::Itertools;
use ptree::{print_tree, TreeBuilder};
use rand::prelude::*;

// pub mod data;
pub mod io;

type Error = Box<dyn std::error::Error>;
type Result<T> = std::result::Result<T, Error>;

/// A Vector backed Tree structure
#[derive(Debug)]
pub struct Tree {
    nodes: Vec<TreeNode>,
    tips: HashSet<usize>,
    is_binary: bool,
    height: Option<f32>,
    diameter: Option<f32>,
}

impl Tree {
    /// Creates a Tree with a single root node
    pub fn new(name: Option<&str>) -> Self {
        Self {
            nodes: vec![TreeNode::new(0, name.map(String::from), None)],
            tips: HashSet::from_iter(vec![0]),
            is_binary: true,
            height: None,
            diameter: None,
        }
    }

    /// Creates an empty Tree
    pub fn new_empty() -> Self {
        Self {
            nodes: vec![],
            tips: HashSet::new(),
            is_binary: true,
            height: None,
            diameter: None,
        }
    }

    /// Adds a root to an empty
    pub fn add_root(&mut self, name: Option<&str>) -> Result<usize> {
        if !self.nodes.is_empty() {
            return Err("Trying to add a root to a non empty tree".into());
        }
        self.nodes
            .push(TreeNode::new(0, name.map(String::from), None));

        self.tips.insert(0);

        Ok(0)
    }

    /// Reset cached values for metrics like tree diameter and height
    pub fn reset_cache(&mut self) {
        self.height = None;
        self.diameter = None;
    }

    /// Creates a node and appends it as a child of the specified parent
    pub fn add_child(&mut self, name: Option<&str>, parent: usize) -> usize {
        self.add_child_node(
            parent,
            TreeNode::new(self.nodes.len(), name.map(String::from), Some(parent)),
        )
    }

    /// Creates a node and appends it as a child of the specified parent
    pub fn add_child_with_len(
        &mut self,
        name: Option<&str>,
        parent: usize,
        len: Option<f32>,
    ) -> usize {
        self.add_child_node(
            parent,
            TreeNode::new_with_length(self.nodes.len(), name.map(String::from), Some(parent), len),
        )
    }

    fn add_child_node(&mut self, parent: usize, mut node: TreeNode) -> usize {
        // This insertion might change the height or diameter of the tree
        self.reset_cache();

        node.distance_to_root = self.get(parent).distance_to_root + 1;

        let idx = self.nodes.len();
        self.nodes.push(node);

        self.nodes[parent].children.push(idx);
        if self.nodes[parent].children.len() > 2 {
            self.is_binary = false
        }

        self.tips.remove(&parent);
        self.tips.insert(idx);

        idx
    }

    /// Returns the size of the tree, i.e. the number of tips
    pub fn size(&self) -> Option<usize> {
        if self.nodes.is_empty() {
            None
        } else {
            Some(self.nodes.len())
        }
    }

    /// Returns height of the tree (i.e. longest distance from tip to root)
    pub fn height(&mut self) -> Option<f32> {
        if self.height.is_some() {
            return self.height;
        }

        if self.nodes.is_empty() {
            None
        } else {
            let height = self
                .tips
                .iter()
                .map(|&tip_index| {
                    let (d, n) = self.get_distance(0, tip_index);
                    match d {
                        Some(d) => d,
                        None => n as f32,
                    }
                })
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            self.height = height;
            height
        }
    }

    /// Gets the diameter of the tree (i.e. the longest distance between tips)
    pub fn diameter(&mut self) -> Option<f32> {
        if self.diameter.is_some() {
            return self.diameter;
        }

        if self.nodes.is_empty() {
            None
        } else {
            let diameter = self
                .tips
                .iter()
                .combinations(2)
                .map(|ids| {
                    let (d, n) = self.get_distance(*ids[0], *ids[1]);
                    match d {
                        Some(d) => d,
                        None => n as f32,
                    }
                })
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            self.diameter = diameter;

            diameter
        }
    }

    /// Computes the sackin statistic for tree with Yule normalization by default
    pub fn sackin(&self) -> Option<usize> {
        if self.is_rooted() && self.is_binary() {
            Some(
                self.tips
                    .iter()
                    .map(|&tip_idx| self.get(tip_idx).distance_to_root)
                    .sum(),
            )
        } else {
            None
        }
    }

    /// Computes the normalized sackin statistic with a Yule null model
    pub fn sackin_yule(&self) -> Option<f64> {
        self.sackin().map(|i_n| {
            let n = self.tips.len();
            let sum: f64 = (2..=n).map(|i| 1.0 / (i as f64)).sum();

            (i_n as f64 - 2.0 * (n as f64) * sum) / n as f64
        })
    }

    /// Computes the normalized sackin statistic with a PDA null model
    pub fn sackin_pda(&self) -> Option<f64> {
        self.sackin()
            .map(|i_n| i_n as f64 / f64::powf(self.tips.len() as f64, 3.0 / 2.0))
    }

    /// Checks if the tree is binary or not
    pub fn is_binary(&self) -> bool {
        self.is_binary
    }

    /// Checks if the tree is rooted
    pub fn is_rooted(&self) -> bool {
        !self.nodes.is_empty() && self.nodes[0].children.len() == 2
    }

    /// Returns indices in the pre-order traversal
    pub fn preorder(&self, root: usize) -> Vec<usize> {
        if root >= self.nodes.len() {
            panic!("Leaf number {root} does not exist in this tree")
        }

        let mut indices = vec![root];
        for child in self.nodes[root].children.iter() {
            indices.extend(self.preorder(*child));
        }

        indices
    }

    /// Returns indices in the post-order traversal
    pub fn postorder(&self, root: usize) -> Vec<usize> {
        if root >= self.nodes.len() {
            panic!("Leaf number {root} does not exist in this tree")
        }

        let mut indices = vec![];
        for child in self.nodes[root].children.iter() {
            indices.extend(self.postorder(*child));
        }
        indices.push(root);

        indices
    }

    /// Returns the indices in the level-order traversal
    pub fn levelorder(&self, root: usize) -> Vec<usize> {
        if root >= self.nodes.len() {
            panic!("Leaf number {root} does not exist in this tree")
        }
        let mut indices = vec![];
        let mut queue = VecDeque::new();
        queue.push_back(root);
        while !queue.is_empty() {
            let root_idx = queue.pop_front().unwrap();
            indices.push(root_idx);
            for child_idx in self.get(root_idx).children.iter() {
                queue.push_back(*child_idx);
            }
        }

        indices
    }

    /// Gets reference to a specified node in the tree
    pub fn get(&self, node: usize) -> &TreeNode {
        &self.nodes[node]
    }

    /// Gets mutable reference to a specified node in the tree
    pub fn get_mut(&mut self, node: usize) -> &mut TreeNode {
        &mut self.nodes[node]
    }

    /// Prunes subtree at given node
    pub fn prune(&mut self, node: usize) {
        let children = self.get(node).children.clone();
        for child in children {
            self.prune(child)
        }
        if let Some(p_idx) = self.get(node).parent {
            self.get_mut(p_idx).children.retain(|&val| val != node);
        }
        let n = self.get_mut(node);
        n.length = None;
        n.parent = None;
        n.children = vec![];
    }

    /// Gets the index of leaf nodes in the tree
    pub fn get_leaves(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .filter(|node| node.children.is_empty())
            .map(|node| node.idx)
            .collect()
    }

    /// Returns the path from the node to the root
    pub fn get_path_from_root(&self, node: usize) -> Vec<usize> {
        let mut path = vec![];
        let mut current_node = node;
        loop {
            path.push(current_node);
            match self.get(current_node).parent {
                Some(parent) => current_node = parent,
                None => break,
            }
        }

        path.into_iter().rev().collect()
    }

    /// Gets the most recent common ancestor between two tree nodes
    pub fn get_common_ancestor(&self, source: usize, target: usize) -> usize {
        if source == target {
            return source;
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

    pub fn get_distance(&self, source: usize, target: usize) -> (Option<f32>, usize) {
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
                if let Some(d) = self.get(*node).length {
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

    /// Recursive function that adds node representation to a printable tree builder
    fn print_nodes(&self, root_idx: usize, output_tree: &mut TreeBuilder, debug: bool) {
        let root = self.get(root_idx);
        if root.children.is_empty() {
            if debug {
                output_tree.add_empty_child(format!("{root:?}"));
            } else {
                output_tree.add_empty_child(root.to_string());
            }
        } else {
            if debug {
                output_tree.begin_child(format!("{root:?}"));
            } else {
                output_tree.begin_child(root.to_string());
            }
            for child_idx in root.children.iter() {
                self.print_nodes(*child_idx, output_tree, debug);
            }
            output_tree.end_child();
        }
    }

    /// Print the tree to the cli
    pub fn print(&self) {
        let mut builder = TreeBuilder::new(self.get(0).to_string());
        for child_idx in self.get(0).children.iter() {
            self.print_nodes(*child_idx, &mut builder, false);
        }
        let tree = builder.build();
        print_tree(&tree).ok();
    }

    /// Print the tree to the cli
    pub fn print_debug(&self) {
        let mut builder = TreeBuilder::new(format!("{:?}", self.get(0)));
        for child_idx in self.get(0).children.iter() {
            self.print_nodes(*child_idx, &mut builder, true);
        }
        let tree = builder.build();
        print_tree(&tree).ok();
    }

    /// Generate newick representation of tree
    fn to_newick_impl(&self, root: usize) -> String {
        if self.get(root).children.is_empty() {
            match self.get(root).name.clone() {
                Some(v) => v,
                None => String::from(""),
            }
        } else {
            "(".to_string()
                + &self
                    .get(root)
                    .children
                    .iter()
                    .map(|child_idx| match self.get(*child_idx).length {
                        Some(l) => format!("{}:{l}", self.to_newick_impl(*child_idx)),
                        None => self.to_newick_impl(*child_idx),
                    })
                    .collect::<Vec<String>>()
                    .join(",")
                + ")"
                + &(match self.get(root).name.clone() {
                    Some(v) => v,
                    None => String::from(""),
                })
        }
    }

    /// Parses a newick string
    /// Heavily inspired (borrowed...) by https://github.com/RagnarGrootKoerkamp/rust-phylogeny
    fn from_newick_impl<'a>(
        mut s: &'a str,
        parent: Option<usize>,
        tree: &mut Tree,
    ) -> Result<&'a str> {
        let node = if let Some(index) = parent {
            // descendant
            tree.add_child(None, index)
        } else {
            // root
            0
        };

        if let Some(_s) = s.strip_prefix('(') {
            s = _s;
            loop {
                s = Self::from_newick_impl(s, Some(node), tree)?;

                let len = if let Some(_s) = s.strip_prefix(':') {
                    s = _s;
                    // The length of the branch ends with , or )
                    let end = s.find(|c| c == ',' || c == ')').ok_or("No , or ) found")?;
                    let len;
                    (len, s) = s.split_at(end);
                    Some(
                        len.parse()
                            .map_err(|e: std::num::ParseFloatError| -> String {
                                "Could not parse length: ".to_string() + &e.to_string()
                            })?,
                    )
                } else {
                    None
                };

                if let Some(child) = tree.get(node).children.last() {
                    tree.get_mut(*child).length = len;
                } else {
                    // Root node has length
                    tree.get_mut(node).length = len;
                }

                let done = s.starts_with(')');
                s = s.split_at(1).1;
                if done {
                    break;
                }
            }
        }

        let end = s
            .find(|c| c == ':' || c == ';' || c == ')' || c == ',')
            .ok_or("No : or ; found")?;
        let (name, s) = s.split_at(end);
        tree.get_mut(node).name = Some(name.to_string());

        Ok(s)
    }

    /// Generate newick representation of tree
    pub fn to_newick(&self) -> String {
        self.to_newick_impl(0) + ";"
    }

    /// Parses a newick string into to Tree
    pub fn from_newick(newick: &str) -> Result<Self> {
        let mut tree = Tree::new(None);
        Self::from_newick_impl(newick, None, &mut tree)?;
        Ok(tree)
    }

    /// Saves the tree to a newick file
    pub fn to_file(&self, path: &Path) -> Result<()> {
        match fs::write(path, self.to_newick()) {
            Ok(_) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    /// returns a preorder traversal iterator
    pub fn iter_preorder(&self) -> PreOrder<'_> {
        PreOrder {
            tree: self,
            indices: vec![0],
        }
    }

    /// returns a postorder traversal iterator
    pub fn iter_postorder(&self) -> PostOrder<'_> {
        PostOrder::new(self)
    }
}

/// A struct to implement the iterator trait on a pre-order  tree traversal
pub struct PreOrder<'a> {
    tree: &'a Tree,
    indices: Vec<usize>,
}

impl<'a> Iterator for PreOrder<'a> {
    type Item = &'a TreeNode;

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, self.tree.size())
    }

    fn next(&mut self) -> Option<Self::Item> {
        if self.indices.is_empty() {
            return None;
        }

        let node = self.indices.pop().unwrap();
        self.tree
            .get(node)
            .children
            .iter()
            .rev()
            .map(|idx| self.indices.push(*idx))
            .count();

        Some(self.tree.get(node))
    }
}

/// A struct to implement the iterator trait on a pre-order  tree traversal
pub struct PostOrder<'a> {
    tree: &'a Tree,
    traversal: Vec<usize>,
}

impl<'a> PostOrder<'a> {
    fn new(tree: &'a Tree) -> Self {
        Self {
            tree,
            traversal: tree.postorder(0).into_iter().rev().collect(),
        }
    }
}

impl<'a> Iterator for PostOrder<'a> {
    type Item = &'a TreeNode;

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, self.tree.size())
    }

    fn next(&mut self) -> Option<Self::Item> {
        let node = self.traversal.pop();
        match node {
            None => None,
            Some(i) => Some(self.tree.get(i)),
        }
    }
}

/// Genereates a random binary tree of a given size
pub fn generate_tree(n_leaves: usize, brlens: bool) -> Tree {
    let mut tree = Tree::new(None);
    let mut rng = thread_rng();

    let mut next_deq = VecDeque::new();
    next_deq.push_back(0);

    for _ in 0..(n_leaves - 1) {
        let parent_idx = if rng.gen_bool(0.5) {
            next_deq.pop_front()
        } else {
            next_deq.pop_back()
        }
        .unwrap();
        let l1: Option<f32> = if brlens { Some(rng.gen()) } else { None };
        let l2: Option<f32> = if brlens { Some(rng.gen()) } else { None };
        next_deq.push_back(tree.add_child_with_len(None, parent_idx, l1));
        next_deq.push_back(tree.add_child_with_len(None, parent_idx, l2));
    }

    for (i, idx) in next_deq.iter().enumerate() {
        tree.get_mut(*idx).set_name(format!("Tip_{i}"));
    }

    tree
}

/// A node of the Tree
pub struct TreeNode {
    /// Index of the node
    pub idx: usize,
    /// Name of the node
    pub name: Option<String>,
    /// Index of the parent node
    pub parent: Option<usize>,
    /// Indices of child nodes
    pub children: Vec<usize>,
    /// Length of branch between parent and node
    pub length: Option<f32>,
    /// Distance to root
    distance_to_root: usize,
}

impl TreeNode {
    /// Creates a new TreeNode
    pub fn new(idx: usize, name: Option<String>, parent: Option<usize>) -> Self {
        Self {
            idx,
            name,
            parent,
            children: vec![],
            length: None,
            distance_to_root: 0,
        }
    }

    /// Creates a new TreeNode with a branch length
    pub fn new_with_length(
        idx: usize,
        name: Option<String>,
        parent: Option<usize>,
        length: Option<f32>,
    ) -> Self {
        Self {
            idx,
            name,
            parent,
            children: vec![],
            length,
            distance_to_root: 0,
        }
    }

    /// Sets the internal TreeNode name
    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }
}

impl Display for TreeNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.length {
            Some(l) => write!(f, "{:?} ({:.3})", self.name, l),
            None => write!(f, "{:?}", self.name),
        }
    }
}

impl Debug for TreeNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?} <I:{}> (L: {:?})[P: {:?}][Root: {:?}]",
            self.name, self.idx, self.length, self.parent, self.distance_to_root
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generates example tree from the tree traversal wikipedia page
    /// https://en.wikipedia.org/wiki/Tree_traversal#Depth-first_search
    fn build_simple_tree() -> Tree {
        let mut tree = Tree::new(Some("F")); // 0
        tree.add_child(Some("B"), 0); // 1
        tree.add_child(Some("G"), 0); // 2
        tree.add_child(Some("A"), 1); // 3
        tree.add_child(Some("D"), 1); // 4
        tree.add_child(Some("I"), 2); // 5
        tree.add_child(Some("C"), 4); // 6
        tree.add_child(Some("E"), 4); // 7
        tree.add_child(Some("H"), 5); // 8

        tree
    }

    /// Generates example tree from the newick format wikipedia page
    /// https://en.wikipedia.org/wiki/Newick_format#Examples
    fn build_tree_with_lengths() -> Tree {
        let mut tree = Tree::new(Some("F")); // 0
        tree.add_child_with_len(Some("A"), 0, Some(0.1)); // 1
        tree.add_child_with_len(Some("B"), 0, Some(0.2)); // 2
        tree.add_child_with_len(Some("E"), 0, Some(0.5)); // 3
        tree.add_child_with_len(Some("C"), 3, Some(0.3)); // 4
        tree.add_child_with_len(Some("D"), 3, Some(0.4)); // 5

        tree
    }

    fn get_values(indices: &[usize], tree: &Tree) -> Vec<Option<String>> {
        indices
            .iter()
            .map(|idx| tree.get(*idx).name.clone())
            .collect()
    }

    #[test]
    fn test_tips() {
        let mut tree = Tree::new(Some("root"));
        assert_eq!(tree.tips, HashSet::from_iter(vec![0]));

        tree.add_child_with_len(Some("A"), 0, Some(0.1)); // 1
        tree.add_child_with_len(Some("B"), 0, Some(0.2)); // 2
        tree.add_child_with_len(Some("E"), 0, Some(0.5)); // 3

        assert_eq!(tree.tips, HashSet::from_iter(vec![1, 2, 3]));

        tree.add_child_with_len(Some("C"), 3, Some(0.3)); // 4
        tree.add_child_with_len(Some("D"), 3, Some(0.4)); // 5

        assert_eq!(tree.tips, HashSet::from_iter(vec![1, 2, 4, 5]));
    }

    #[test]
    fn test_binary() {
        let mut tree = Tree::new(Some("root"));
        tree.add_child(Some("0L"), 0); //1
        tree.add_child(Some("0R"), 0); //2

        assert!(tree.is_binary());

        tree.add_child(Some("1L"), 1); //3
        tree.add_child(Some("1R"), 1); //4

        assert!(tree.is_binary());

        tree.add_child(Some("3L"), 3); //5
        tree.add_child(Some("3R"), 3); //6
        assert!(tree.is_binary());

        tree.add_child(Some("3?"), 3); //7
        assert!(!tree.is_binary());
    }

    #[test]
    fn binary_from_newick() {
        let test_cases = vec![("(A,B,(C,D)E)F;", false), ("((D,E)B,(F,G)C)A;", true)];

        for (newick, is_binary) in test_cases {
            assert_eq!(Tree::from_newick(newick).unwrap().is_binary(), is_binary)
        }
    }

    #[test]
    fn traverse_preorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.preorder(0)), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["F", "B", "A", "D", "C", "E", "G", "I", "H"])
    }

    #[test]
    fn iter_preorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = tree
            .iter_preorder()
            .flat_map(|node| node.name.clone())
            .collect();
        assert_eq!(values, vec!["F", "B", "A", "D", "C", "E", "G", "I", "H"])
    }

    #[test]
    fn traverse_postorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.postorder(0)), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"])
    }

    #[test]
    fn iter_postorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = tree
            .iter_postorder()
            .flat_map(|node| node.name.clone())
            .collect();
        assert_eq!(values, vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"])
    }

    #[test]
    fn traverse_levelorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.levelorder(0)), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["F", "B", "G", "A", "D", "I", "C", "E", "H"])
    }

    #[test]
    fn prune_tree() {
        let mut tree = build_simple_tree();
        tree.prune(4); // prune D subtree
        let values: Vec<_> = get_values(&(tree.preorder(0)), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["F", "B", "A", "G", "I", "H"]);
    }

    #[test]
    fn path_from_root() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.get_path_from_root(7)), &tree)
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
        let tree = build_simple_tree();
        for ((source, target), ancestor) in test_cases {
            println!(
                "Testing: ({:?}, {:?}) -> {:?}",
                tree.get(source).name,
                tree.get(target).name,
                tree.get(ancestor).name
            );
            assert_eq!(ancestor, tree.get_common_ancestor(source, target));
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
        let tree = build_tree_with_lengths();

        for ((idx_s, idx_t), (dist, branches)) in test_cases {
            let (d_pred, b_pred) = tree.get_distance(idx_s, idx_t);
            assert_eq!(branches, b_pred);
            match dist {
                None => assert!(d_pred.is_none()),
                Some(d) => {
                    assert!(d_pred.is_some_and(|x| (x - d).abs() < f32::EPSILON))
                }
            }
        }
    }

    #[test]
    fn get_correct_leaves() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.get_leaves()), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["A", "C", "E", "H"])
    }

    #[test]
    fn generate_random_correct_size() {
        use rand::prelude::*;
        let mut rng = thread_rng();

        for size in (0..20).map(|_| rng.gen_range(10..=100)) {
            let tree = generate_tree(size, false);
            assert_eq!(tree.get_leaves().len(), size);
        }
    }

    #[test]
    fn to_newick() {
        let tree = build_tree_with_lengths();
        assert_eq!("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;", tree.to_newick());
    }

    // test cases from https://github.com/ila/Newick-validator
    // For now it does not parse the root node branch length...
    #[test]
    fn read_newick() {
        let newick_strings = vec![
            "((D,E)B,(F,G)C)A;",
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
            // "(dog:20,(elephant:30,horse:60):20):50;",
            "(A,B,(C,D));",
            "(A,B,(C,D)E)F;",
            // "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;",
            "(:0.1,:0.2,(:0.3,:0.4):0.5);",
            // "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;",
            "(,,(,));",
        ];
        for newick in newick_strings {
            let tree = Tree::from_newick(newick).unwrap();
            assert_eq!(newick, tree.to_newick());
        }
    }

    #[test]
    fn test_height() {
        // heights computed with ete3
        let test_cases = vec![
            ("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;", 0.9),
            ("((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;", 1.0),
            ("(A,B,(C,D)E)F;", 2.0),
            (
                "((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip0,Tip1);",
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
            assert!(sackin.is_some());
            assert_eq!(tree.sackin().unwrap(), true_sackin);
        }
    }

    #[test]
    fn test_sackin_unrooted() {
        let test_cases = vec![
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;",
            "(A,B,(C,D)E)F;",
            "((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip0,Tip1);",
        ];

        for newick in test_cases {
            let tree = Tree::from_newick(newick).unwrap();
            assert!(tree.sackin().is_none());
        }
    }
}
