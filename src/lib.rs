use std::{
    collections::VecDeque,
    fmt::{Debug, Display}, path::Path, fs,
};

use ptree::{print_tree, TreeBuilder};
use rand::prelude::*;

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

    /// Gets reference to a specified node in the tree
    pub fn get(&self, node: usize) -> &TreeNode<T> {
        &self.nodes[node]
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

    fn build_simple_tree() -> Tree<i32> {
        let mut tree = Tree::new(1);
        tree.add_child(2, 0); // root.left (1)
        tree.add_child(3, 0); // root.right (2)
        tree.add_child(4, 1); // root.left.left (3)
        tree.add_child(5, 1); // root.left.right (4)

        tree
    }

    /// Generates example tree from the newick format wikipedia page
    /// https://en.wikipedia.org/wiki/Newick_format#Examples
    fn build_tree_with_lengths() -> Tree<&'static str> {
        let mut tree = Tree::new("F");
        tree.add_child_with_len("A", 0, Some(0.1));
        tree.add_child_with_len("B", 0, Some(0.2));
        tree.add_child_with_len("E", 0, Some(0.5));
        tree.add_child_with_len("C", 3, Some(0.3));
        tree.add_child_with_len("D", 3, Some(0.4));

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
        assert_eq!(values, vec![1, 2, 4, 5, 3])
    }

    #[test]
    fn traverse_postorder() {
        let tree = build_simple_tree();
        let values = get_values(&(tree.postorder(0)), &tree);
        assert_eq!(values, vec![4, 5, 2, 3, 1])
    }

    #[test]
    fn get_correct_leaves() {
        let tree = build_simple_tree();
        let values = get_values(&(tree.get_leaves()), &tree);
        assert_eq!(values, vec![3, 4, 5])
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
