#![feature(is_some_and)]

use std::{
    collections::VecDeque,
    fmt::{Debug, Display},
    fs,
    iter::zip,
    path::Path,
};

use ptree::{print_tree, TreeBuilder};
use rand::prelude::*;

/// A Vector backed Tree structure
#[derive(Debug)]
pub struct Tree<T> {
    nodes: Vec<TreeNode<T>>,
}

impl<T> Tree<T> {
    /// Creates a Tree with a single root node
    pub fn new(val: T) -> Self {
        Self {
            nodes: vec![TreeNode::new(0, val, None)],
        }
    }

    /// Creates a node and appends it as a child of the specified parent
    pub fn add_child(&mut self, val: T, parent: usize) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(TreeNode::new(idx, val, Some(parent)));
        self.nodes[parent].children.push(idx);
        idx
    }

    /// Creates a node and appends it as a child of the specified parent
    pub fn add_child_with_len(&mut self, val: T, parent: usize, len: Option<f32>) -> usize {
        let idx = self.nodes.len();
        self.nodes
            .push(TreeNode::new_with_length(idx, val, Some(parent), len));
        self.nodes[parent].children.push(idx);
        idx
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
    pub fn get(&self, node: usize) -> &TreeNode<T> {
        &self.nodes[node]
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

    /// Gets mutable reference to a specified node in the tree
    pub fn get_mut(&mut self, node: usize) -> &mut TreeNode<T> {
        &mut self.nodes[node]
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
}

impl<T> Tree<T>
where
    T: Debug,
{
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
}

impl<T> Tree<T>
where
    T: Display,
{
    /// Generate newick representation of tree
    fn to_newick_impl(&self, root: usize) -> String {
        if self.get(root).children.is_empty() {
            self.get(root).val.to_string()
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
                + &(self.get(root).val.to_string())
        }
    }

    /// Generate newick representation of tree
    pub fn to_newick(&self) -> String {
        self.to_newick_impl(0) + ";"
    }

    /// Saves the tree to a newick file
    pub fn to_file(&self, path: &Path) {
        fs::write(path, self.to_newick()).unwrap()
    }
}

/// Genereates a random binary tree of a given size
pub fn generate_tree(n_leaves: usize, brlens: bool) -> Tree<String> {
    let mut tree = Tree::new(String::from("root"));
    let mut rng = thread_rng();

    let mut next_deq = VecDeque::new();
    next_deq.push_back(0);

    let mut counter = 1;
    for _ in 0..(n_leaves - 1) {
        let parent_idx = if rng.gen_bool(0.5) {
            next_deq.pop_front()
        } else {
            next_deq.pop_back()
        }
        .unwrap();
        let l1: Option<f32> = if brlens { Some(rng.gen()) } else { None };
        let l2: Option<f32> = if brlens { Some(rng.gen()) } else { None };
        next_deq.push_back(tree.add_child_with_len(format!("Node_{counter}"), parent_idx, l1));
        next_deq.push_back(tree.add_child_with_len(
            format!("Node_{}", counter + 1),
            parent_idx,
            l2,
        ));
        counter += 2;
    }

    for (i, idx) in next_deq.iter().enumerate() {
        tree.get_mut(*idx).set_val(format!("Tip_{i}"));
    }

    tree
}

/// A node of the Tree
pub struct TreeNode<T> {
    /// Index of the node
    pub idx: usize,
    /// Value stored in the node (a name)
    pub val: T,
    /// Index of the parent node
    pub parent: Option<usize>,
    /// Indices of child nodes
    children: Vec<usize>,
    /// Length of branch between parent and node
    pub length: Option<f32>,
}

impl<T> TreeNode<T> {
    /// Creates a new TreeNode
    pub fn new(idx: usize, val: T, parent: Option<usize>) -> Self {
        Self {
            idx,
            val,
            parent,
            children: vec![],
            length: None,
        }
    }

    /// Creates a new TreeNode with a branch length
    pub fn new_with_length(idx: usize, val: T, parent: Option<usize>, length: Option<f32>) -> Self {
        Self {
            idx,
            val,
            parent,
            children: vec![],
            length,
        }
    }

    /// Sets the internal TreeNode value
    pub fn set_val(&mut self, val: T) {
        self.val = val;
    }
}

impl<T> Display for TreeNode<T>
where
    T: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.length {
            Some(l) => write!(f, "{:?} ({:.3})", self.val, l),
            None => write!(f, "{:?}", self.val),
        }
    }
}

impl<T> Debug for TreeNode<T>
where
    T: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?} <I:{}> (L: {:?})[P: {:?}]",
            self.val, self.idx, self.length, self.parent
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generates example tree from the tree traversal wikipedia page
    /// https://en.wikipedia.org/wiki/Tree_traversal#Depth-first_search
    fn build_simple_tree() -> Tree<&'static str> {
        let mut tree = Tree::new("F"); // 0
        tree.add_child("B", 0); // 1
        tree.add_child("G", 0); // 2
        tree.add_child("A", 1); // 3
        tree.add_child("D", 1); // 4
        tree.add_child("I", 2); // 5
        tree.add_child("C", 4); // 6
        tree.add_child("E", 4); // 7
        tree.add_child("H", 5); // 8

        tree
    }

    /// Generates example tree from the newick format wikipedia page
    /// https://en.wikipedia.org/wiki/Newick_format#Examples
    fn build_tree_with_lengths() -> Tree<&'static str> {
        let mut tree = Tree::new("F"); // 0
        tree.add_child_with_len("A", 0, Some(0.1)); // 1
        tree.add_child_with_len("B", 0, Some(0.2)); // 2
        tree.add_child_with_len("E", 0, Some(0.5)); // 3
        tree.add_child_with_len("C", 3, Some(0.3)); // 4
        tree.add_child_with_len("D", 3, Some(0.4)); // 5

        tree
    }

    fn get_values<T>(indices: &[usize], tree: &Tree<T>) -> Vec<T>
    where
        T: Clone,
        T: Debug,
    {
        indices
            .iter()
            .map(|idx| tree.get(*idx).val.clone())
            .collect()
    }

    #[test]
    fn traverse_preorder() {
        let tree = build_simple_tree();
        let values = get_values(&(tree.preorder(0)), &tree);
        assert_eq!(values, vec!["F", "B", "A", "D", "C", "E", "G", "I", "H"])
    }

    #[test]
    fn traverse_postorder() {
        let tree = build_simple_tree();
        let values = get_values(&(tree.postorder(0)), &tree);
        assert_eq!(values, vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"])
    }

    #[test]
    fn traverse_levelorder() {
        let tree = build_simple_tree();
        let values = get_values(&(tree.levelorder(0)), &tree);
        assert_eq!(values, vec!["F", "B", "G", "A", "D", "I", "C", "E", "H"])
    }

    #[test]
    fn prune_tree() {
        let mut tree = build_simple_tree();
        tree.prune(4); // prune D subtree
        let values = get_values(&(tree.preorder(0)), &tree);
        assert_eq!(values, vec!["F", "B", "A", "G", "I", "H"]);
    }

    #[test]
    fn path_from_root() {
        let tree = build_simple_tree();
        let values = get_values(&(tree.get_path_from_root(7)), &tree);
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
                "Testing: ({}, {}) -> {}",
                tree.get(source).val,
                tree.get(target).val,
                tree.get(ancestor).val
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
        let values = get_values(&(tree.get_leaves()), &tree);
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
}
