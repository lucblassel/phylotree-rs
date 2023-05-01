// #![feature(is_some_and)]

use std::{
    cell::RefCell,
    collections::{HashMap, HashSet, VecDeque},
    fmt::{Debug, Display},
    fs,
    iter::zip,
    path::Path,
};

use distance::DistanceMatrix;
use distr::{Distr, Sampler};
use errors::{ParseError, TreeError};
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use ptree::{print_tree, TreeBuilder};
use rand::prelude::*;

// pub mod data;
pub mod distance;
pub mod distr;
pub mod errors;

type Error = Box<dyn std::error::Error>;
type ResultStd<T, E> = std::result::Result<T, E>;
type Result<T> = std::result::Result<T, Error>;
type Edge = usize;

/// A Vector backed Tree structure
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<TreeNode>,
    tips: HashSet<usize>,
    is_binary: bool,
    is_rooted: bool,
    height: RefCell<Option<f32>>,
    diameter: RefCell<Option<f32>>,
    pub leaf_index: RefCell<Option<Vec<String>>>,
    partitions: RefCell<Option<HashMap<usize, Option<f32>>>>,
}

impl Tree {
    /// Creates a Tree with a single root node
    pub fn new(name: Option<&str>) -> Self {
        Self {
            nodes: vec![TreeNode::new(0, name.map(String::from), None)],
            tips: HashSet::from_iter(vec![0]),
            is_binary: true,
            is_rooted: true,
            height: RefCell::new(None),
            diameter: RefCell::new(None),
            leaf_index: RefCell::new(None),
            partitions: RefCell::new(None),
        }
    }

    /// Creates an empty Tree
    pub fn new_empty() -> Self {
        Self {
            nodes: vec![],
            tips: HashSet::new(),
            is_binary: true,
            is_rooted: true,
            height: RefCell::new(None),
            diameter: RefCell::new(None),
            leaf_index: RefCell::new(None),
            partitions: RefCell::new(None),
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
        self.height = RefCell::new(None);
        self.diameter = RefCell::new(None);
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
        self.nodes[parent].set_internal();
        if self.nodes[parent].children.len() > 2 {
            self.is_binary = false
        }

        if parent == 0 && self.nodes[parent].children.len() == 3 {
            self.is_rooted = false
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
    pub fn height(&self) -> Option<f32> {
        if (*self.height.borrow()).is_some() {
            return *self.height.borrow();
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

            (*self.height.borrow_mut()) = height;
            height
        }
    }

    /// Gets the diameter of the tree (i.e. the longest distance between tips)
    pub fn diameter(&self) -> Option<f32> {
        if (*self.diameter.borrow()).is_some() {
            return *self.diameter.borrow();
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

            (*self.diameter.borrow_mut()) = diameter;

            diameter
        }
    }

    /// Computes the number of cherries in a tree
    pub fn cherries(&self) -> Result<usize> {
        self.check_rooted_binary()?;
        if !self.nodes.is_empty() {
            let mut n = 0;
            for node in self.nodes.iter() {
                if node.children.len() == 2
                    && self.get(node.children[0]).tip
                    && self.get(node.children[1]).tip
                {
                    n += 1;
                }
            }
            Ok(n)
        } else {
            Err(Box::new(TreeError::IsEmpty))
        }
    }

    /// Computes the Colless index for the tree
    pub fn colless(&self) -> Result<usize> {
        self.check_rooted_binary()?;

        let colless = self
            .nodes
            .iter()
            .filter(|node| !node.is_tip())
            .map(|node| {
                if node.children.is_empty() {
                    return 0;
                }
                let left = self.get_subtree_leaves(node.children[0]).len();
                let right = if node.children.len() > 1 {
                    self.get_subtree_leaves(node.children[1]).len()
                } else {
                    0
                };

                left.abs_diff(right)
            })
            .sum();

        Ok(colless)
    }

    /// Computes the normalized colless statistic with a Yule null model
    pub fn colless_yule(&self) -> Result<f64> {
        self.colless().map(|i_c| {
            let n = self.tips.len() as f64;
            let e_i_c = n * n.ln() + (0.57721566 - 1. - f64::ln(2.0)) * n;

            (i_c as f64 - e_i_c) / n
        })
    }

    /// Computes the normalized colless statistic with a PDA null model
    pub fn colless_pda(&self) -> Result<f64> {
        self.colless()
            .map(|i_c| i_c as f64 / f64::powf(self.tips.len() as f64, 3.0 / 2.0))
    }

    /// Computes the sackin statistic for tree
    pub fn sackin(&self) -> Result<usize> {
        self.check_rooted_binary()?;

        Ok(self
            .tips
            .iter()
            .map(|&tip_idx| self.get(tip_idx).distance_to_root)
            .sum())
    }

    /// Computes the normalized sackin statistic with a Yule null model
    pub fn sackin_yule(&self) -> Result<f64> {
        self.sackin().map(|i_n| {
            let n = self.tips.len();
            let sum: f64 = (2..=n).map(|i| 1.0 / (i as f64)).sum();

            (i_n as f64 - 2.0 * (n as f64) * sum) / n as f64
        })
    }

    /// Computes the normalized sackin statistic with a PDA null model
    pub fn sackin_pda(&self) -> Result<f64> {
        self.sackin()
            .map(|i_n| i_n as f64 / f64::powf(self.tips.len() as f64, 3.0 / 2.0))
    }

    /// Checks if the tree is rooted and binary
    pub fn check_rooted_binary(&self) -> Result<()> {
        if !self.is_rooted() {
            Err(Box::new(TreeError::IsNotRooted))
        } else if !self.is_binary() {
            Err(Box::new(TreeError::IsNotBinary))
        } else {
            Ok(())
        }
    }

    /// Checks if the tree is binary or not
    pub fn is_binary(&self) -> bool {
        self.is_binary
    }

    /// Checks if the tree is rooted
    pub fn is_rooted(&self) -> bool {
        !self.nodes.is_empty() && self.nodes[0].children.len() == 2
    }

    /// Checks if the tree tips are uniquely named
    pub fn are_tip_names_unique(&self) -> bool {
        let n_names = self
            .tips
            .iter()
            .filter_map(|&idx| self.get(idx).name.clone())
            .sorted()
            .dedup()
            .count();

        n_names == self.tips.len()
    }

    /// Returns indices in the pre-order traversal
    pub fn preorder(&self, root: usize) -> Result<Vec<usize>> {
        if root >= self.nodes.len() {
            return Err(format!("Leaf number {root} does not exist in this tree").into());
        }

        let mut indices = vec![root];
        for child in self.nodes[root].children.iter() {
            indices.extend(self.preorder(*child)?);
        }

        Ok(indices)
    }

    /// Returns indices in the post-order traversal
    pub fn postorder(&self, root: usize) -> Result<Vec<usize>> {
        if root >= self.nodes.len() {
            return Err(format!("Leaf number {root} does not exist in this tree").into());
        }

        let mut indices = vec![];
        for child in self.nodes[root].children.iter() {
            indices.extend(self.postorder(*child)?);
        }
        indices.push(root);

        Ok(indices)
    }

    pub fn inorder(&self, root: usize) -> Result<Vec<usize>> {
        if root >= self.nodes.len() {
            return Err(format!("Leaf number {root} does not exist in this tree").into());
        }
        if !self.is_binary() {
            return Err("In-order traversal is only possible on binary trees.".into());
        }

        let mut indices = vec![];

        let children = self.nodes[root].children.clone();
        // There is a left child
        if !children.is_empty() {
            indices.extend(self.inorder(children[0])?);
        }
        indices.push(root);
        // There is a right child
        if children.len() > 1 {
            indices.extend(self.inorder(children[1])?);
        }
        Ok(indices)
    }

    /// Returns the indices in the level-order traversal
    pub fn levelorder(&self, root: usize) -> Result<Vec<usize>> {
        if root >= self.nodes.len() {
            return Err(format!("Leaf number {root} does not exist in this tree").into());
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

        Ok(indices)
    }

    /// Rescales the branch lengths of a tree
    pub fn rescale(&mut self, scale: f32) {
        if let Some(diam) = self.diameter() {
            for node in self.nodes.iter_mut().skip(1) {
                node.length = node.length.map(|l| l * scale / diam);
            }
        }
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

    /// Remove nodes with only one child
    pub fn compress(&mut self) -> Result<()> {
        for node in self.nodes.clone().into_iter() {
            if node.children.len() != 1 {
                continue;
            }

            let parent = self.get(node.idx).parent;
            let parent_len = self.get(node.idx).length;
            let child_idx = self.get(node.idx).children[0];
            let child_len = self.get(child_idx).length;

            // Compress tree
            if let Some(parent_idx) = parent {
                // Make parent node point to new child
                let parent_node = self.get_mut(parent_idx);
                let parent_node_idx = parent_node
                    .children
                    .iter()
                    .position(|&idx| idx == node.idx)
                    .unwrap();
                parent_node.children[parent_node_idx] = child_idx;
            }

            let child = self.get_mut(child_idx);
            child.parent = parent;
            child.length = match (parent_len, child_len) {
                (Some(d1), Some(d2)) => Some(d1 + d2),
                (None, None) => None,
                (Some(d), None) | (None, Some(d)) => Some(d),
            };

            let node = self.get_mut(node.idx);
            node.length = None;
            node.parent = None;
            node.children = vec![];
        }
        Ok(())
    }

    /// Gets the index of leaf nodes in the tree
    pub fn get_leaves(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .filter(|node| node.children.is_empty())
            .map(|node| node.idx)
            .collect()
    }

    /// Get indices of descendant nodes.
    pub fn get_descendants(&self, index: usize) -> Vec<usize> {
        let mut indices = vec![];
        if self.nodes[index].is_tip() {
            return indices;
        }

        indices.extend(self.nodes[index].children.clone());

        for &child_idx in self.nodes[index].children.iter() {
            indices.extend(self.get_descendants(child_idx))
        }

        indices
    }

    /// Gets the leaf indices in the subtree rooted at note of a specified index
    pub fn get_subtree_leaves(&self, index: usize) -> Vec<usize> {
        let mut indices = vec![];
        if self.nodes[index].is_tip() {
            return vec![index];
        }

        for &child_idx in self.nodes[index].children.iter() {
            indices.extend(self.get_subtree_leaves(child_idx))
        }

        indices
    }

    /// Get left subtree
    pub fn get_left(&self, index: usize) -> Result<Vec<usize>> {
        if !self.is_binary() {
            return Err(Box::new(TreeError::IsNotBinary));
        }
        if self.nodes[index].is_tip() {
            return Ok(vec![]);
        }

        Ok(self.get_descendants(self.nodes[index].children[0]))
    }

    /// Get right subtree
    pub fn get_right(&self, index: usize) -> Result<Vec<usize>> {
        if !self.is_binary() {
            return Err(Box::new(TreeError::IsNotBinary));
        }
        if self.nodes[index].is_tip() {
            return Ok(vec![]);
        }

        Ok(self.get_descendants(self.nodes[index].children[1]))
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

    /// Gets the distance between 2 nodes, returns the sum of branch lengths (if all
    /// branches in the path have lengths) and the numebr of edges in the path
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

    /// Get Vec containing leaf names
    pub fn get_leaf_names(&self) -> Vec<String> {
        self.tips
            .iter()
            .filter_map(|&tip_index| self.get(tip_index).name.clone())
            .collect()
    }

    /// Initializes the leaf index
    fn init_leaf_index(&self) -> Result<()> {
        if self.nodes.is_empty() {
            return Err(Box::new(TreeError::IsEmpty));
        }
        if self.leaf_index.borrow().is_some() {
            return Ok(());
        }

        let names = self.get_leaf_names();
        if names.len() != self.tips.len() {
            return Err(Box::new(TreeError::UnnamedLeaves));
        }

        if !self.are_tip_names_unique() {
            return Err(Box::new(TreeError::DuplicateLeafNames));
        }

        (*self.leaf_index.borrow_mut()) = Some(names.into_iter().sorted().collect());

        Ok(())
    }

    /// Resets the leaf index
    pub fn reset_leaf_index(&self) {
        (*self.leaf_index.borrow_mut()) = None;
    }

    /// Get the partition corresponding to the branch associated to the node at index
    pub fn get_partition(&self, index: usize) -> Result<Edge> {
        self.init_leaf_index()?;
        let indices = self
            .get_subtree_leaves(index)
            .into_iter()
            .filter_map(|index| self.get(index).name.clone())
            .map(|name| {
                let v = self.leaf_index.borrow().clone();
                v.map(|v| v.iter().position(|n| *n == name).unwrap())
                    .unwrap()
            });

        let mut hash = 0;
        for index in indices {
            hash ^= 1 << index;
        }

        Ok(hash)
    }

    /// Get the partition corresponding to the branch associated to the node at index
    pub fn get_partition_new(&self, index: usize) -> Result<FixedBitSet> {
        self.init_leaf_index()?;

        let indices = self
            .get_subtree_leaves(index)
            .into_iter()
            .filter_map(|index| self.get(index).name.clone())
            .map(|name| {
                let v = self.leaf_index.borrow().clone();
                v.map(|v| v.iter().position(|n| *n == name).unwrap())
                    .unwrap()
            });

        let mut bitset = FixedBitSet::with_capacity(self.tips.len());
        for index in indices {
            bitset.insert(index);
        }

        let mut toggled = bitset.clone();
        toggled.toggle_range(..);

        Ok(toggled.min(bitset))
    }

    /// Caches partitions for distance computation
    fn init_partitions_new(&self) -> Result<HashMap<FixedBitSet, Option<f32>>> {
        self.init_leaf_index()?;

        let mut partitions: HashMap<FixedBitSet, Option<f32>> = HashMap::new();

        for node in self
            .nodes
            .iter()
            .filter(|n| !(n.deleted || n.parent.is_none() || n.is_tip()))
        {
            let part = self.get_partition_new(node.idx)?;

            if part.count_ones(..) == 1 {
                continue;
            }

            let new_len = node.length;
            let old_len = partitions.get(&part);

            let len = match (new_len, old_len) {
                (None, None) => None,
                (Some(new_len), Some(old_len)) => old_len.map(|v| v + new_len),
                (Some(new_len), None) => Some(new_len),
                (None, Some(old_len)) => *old_len,
            };

            partitions.insert(part, len);
        }

        Ok(partitions)
    }

    fn get_partitions_new(&self) -> Result<HashSet<FixedBitSet>> {
        self.init_leaf_index()?;

        let partitions = self.init_partitions_new()?;

        Ok(HashSet::from_iter(partitions.keys().cloned()))
    }

    /// Computes the Robinson Foulds distance between two trees
    pub fn robinson_foulds_new(&self, other: &Self) -> Result<usize> {
        let partitions_s = self.get_partitions_new()?;
        let partitions_o = other.get_partitions_new()?;

        if *(self.leaf_index.borrow()) != *(other.leaf_index.borrow()) {
            return Err(TreeError::DifferentTipIndices.into());
        }

        let mut root_s = HashSet::new();
        for i in self.nodes[0].children.iter() {
            root_s.insert(self.get_partition_new(*i)?);
        }
        let mut root_o = HashSet::new();
        for i in other.nodes[0].children.iter() {
            root_o.insert(other.get_partition_new(*i)?);
        }

        let same_root = root_s == root_o;

        let i = partitions_o.intersection(&partitions_s).count();
        let rf = partitions_o.len() + partitions_s.len() - 2 * i;

        // Hacky...
        if self.is_rooted && rf != 0 && !same_root {
            Ok(rf + 2)
        } else {
            Ok(rf)
        }
    }

    /// Computes the normalized Robinson Foulds distance between two trees
    pub fn robinson_foulds_norm_new(&self, other: &Self) -> Result<f32> {
        let rf = self.robinson_foulds_new(other)?;

        let partitions_s = self.get_partitions_new()?;
        let partitions_o = other.get_partitions_new()?;

        let tot = partitions_o.len() + partitions_s.len();

        Ok((rf as f32) / (tot as f32))
    }

    /// Caches partitions for distance computation
    fn init_partitions(&self) -> Result<()> {
        self.init_leaf_index()?;

        if self.partitions.borrow().is_some() {
            return Ok(());
        }

        let mut partitions: HashMap<Edge, Option<f32>> = HashMap::new();

        let m: usize = 2u64.pow(self.tips.len() as u32) as usize - 1;
        for node in self
            .nodes
            .iter()
            .filter(|n| !(n.deleted || n.parent.is_none() || n.is_tip()))
        {
            let part = self.get_partition(node.idx)?;
            let part = part.min(m ^ part);

            let new_len = node.length;
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

    pub fn get_partitions(&self) -> Result<HashSet<Edge>> {
        self.init_leaf_index()?;

        if self.partitions.borrow().is_none() {
            self.init_partitions()?;
        }

        Ok(HashSet::from_iter(
            self.partitions
                .borrow()
                .as_ref()
                .unwrap()
                .iter()
                .map(|(k, _)| *k),
        ))
    }

    /// Get all partitions of a tree
    pub fn get_partitions_with_lengths(&self) -> Result<HashMap<Edge, f32>> {
        self.init_leaf_index()?;

        if self.partitions.borrow().is_none() {
            self.init_partitions()?;
        }

        let mut partitions = HashMap::new();
        for (hash, len) in self.partitions.borrow().as_ref().unwrap().iter() {
            if len.is_none() {
                return Err(TreeError::MissingBranchLengths.into());
            }
            partitions.insert(*hash, len.unwrap());
        }

        Ok(partitions)
    }

    /// Computes the Robinson Foulds distance between two trees
    pub fn robinson_foulds(&self, other: &Self) -> Result<usize> {
        let partitions_s = self.get_partitions()?;
        let partitions_o = other.get_partitions()?;

        if *(self.leaf_index.borrow()) != *(other.leaf_index.borrow()) {
            return Err(TreeError::DifferentTipIndices.into());
        }

        let i = partitions_o.intersection(&partitions_s).count();

        Ok(partitions_o.len() + partitions_s.len() - 2 * i)
    }

    /// Computes the normalized Robinson Foulds distance between two trees
    pub fn robinson_foulds_norm(&self, other: &Self) -> Result<f32> {
        let partitions_s = self.get_partitions()?;
        let partitions_o = other.get_partitions()?;

        if *(self.leaf_index.borrow()) != *(other.leaf_index.borrow()) {
            return Err(TreeError::DifferentTipIndices.into());
        }

        let i = partitions_o.intersection(&partitions_s).count();

        let tot = partitions_o.len() + partitions_s.len();
        let rf = tot - 2 * i;

        Ok((rf as f32) / (tot as f32))
    }

    /// Computes the weighted Robinson Foulds distance between two trees
    pub fn weighted_robinson_foulds(&self, other: &Self) -> Result<f32> {
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

    /// Computes the khuner felsenstein branch score between two trees
    pub fn khuner_felsenstein(&self, other: &Self) -> Result<f32> {
        let partitions_s = self.get_partitions_with_lengths()?;
        let partitions_o = other.get_partitions_with_lengths()?;

        let mut dist = 0.;

        for (edge, len_s) in partitions_s.iter() {
            if let Some(len_o) = partitions_o.get(edge) {
                dist += f32::powi(len_s - len_o, 2)
            } else {
                dist += f32::powi(*len_s, 2)
            }
        }

        for (edge, len_o) in partitions_o.iter() {
            if !partitions_s.contains_key(edge) {
                dist += f32::powi(*len_o, 2)
            }
        }

        Ok(dist.sqrt())
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
        let root = self.get(root);
        if root.children.is_empty() {
            root.to_newick()
        } else {
            "(".to_string()
                + &(root
                    .children
                    .iter()
                    .map(|child_idx| self.to_newick_impl(*child_idx)))
                .collect::<Vec<String>>()
                .join(",")
                + ")"
                + &(root.to_newick())
        }
    }

    fn distance_matrix_recursive_impl(
        &self,
        current: usize,
        prev: Option<usize>,
        lengths: &mut [f32],
        currlength: f32,
    ) -> Result<()> {
        if prev.is_some() && self.tips.contains(&current) {
            lengths[current] = currlength;
            return Ok(());
        }

        let mut neighbors: Vec<_> = self
            .get(current)
            .children
            .clone()
            .into_iter()
            .map(|idx| (idx, self.get(idx).length))
            .collect();
        if let Some(parent) = self.get(current).parent {
            neighbors.push((parent, self.get(current).length))
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
                    return Err(TreeError::MissingBranchLengths.into());
                }
            }
        }

        Ok(())
    }

    pub fn distance_matrix_recursive(&self) -> Result<DistanceMatrix<f32>> {
        let size = self.nodes.len();
        let mut matrix = DistanceMatrix::new(self.tips.len());
        let mut cache: Vec<Vec<_>> = vec![vec![f32::INFINITY; size]; size];

        for tip in self.tips.iter() {
            self.distance_matrix_recursive_impl(*tip, None, &mut cache[*tip], 0.0)?
        }

        for pair in self.tips.iter().combinations(2) {
            let (i1, i2) = (pair[0], pair[1]);
            let d = cache[*i1][*i2];
            let name1 = self.get(*i1).name.clone().unwrap();
            let name2 = self.get(*i2).name.clone().unwrap();

            matrix.set(&name1, &name2, d, false)?;
        }

        Ok(matrix)
    }

    /// Compute distance matrix from tree without a cache
    pub fn distance_matrix(&self) -> Result<DistanceMatrix<f32>> {
        let mut matrix = DistanceMatrix::new(self.tips.len());

        for pair in self.tips.iter().combinations(2) {
            let (i1, i2) = (pair[0], pair[1]);
            if let (Some(d), _) = self.get_distance(*i1, *i2) {
                let name1 = self.get(*i1).name.clone().unwrap();
                let name2 = self.get(*i2).name.clone().unwrap();

                matrix.set(&name1, &name2, d, false)?;
            } else {
                return Err(TreeError::MissingBranchLengths.into());
            }
        }

        Ok(matrix)
    }

    /// Generate newick representation of tree
    pub fn to_newick(&self) -> String {
        self.to_newick_impl(0) + ";"
    }

    /// Parses a newick string into to Tree
    pub fn from_newick(newick: &str) -> Result<Self> {
        #[derive(Debug)]
        enum Field {
            Name,
            Length,
            Comment,
        }

        let mut tree = Tree::new_empty();

        let mut parsing = Field::Name;
        let mut current_name: Option<String> = None;
        let mut current_length: Option<String> = None;
        let mut current_index = None;
        let mut parent_stack = Vec::new();

        let mut open_delimiters = Vec::new();

        for c in newick.chars() {
            // eprintln!("Parsing: {c}");
            match c {
                '(' => {
                    // Start subtree
                    match parent_stack.last() {
                        None => parent_stack.push(tree.add_root(None)?),
                        Some(parent) => parent_stack.push(tree.add_child(None, *parent)),
                    };
                    open_delimiters.push(0);
                    // eprintln!(
                    //     "\t\tOPEN: Adding empty node to parent stack -> index {:?}",
                    //     parent_stack.last()
                    // )
                }
                ':' => {
                    parsing = Field::Length;
                }
                ',' => {
                    // eprintln!("\t\tSIBL: Adding name and length to current node {current_index:?}");
                    let node = if let Some(index) = current_index {
                        tree.get_mut(index)
                    } else {
                        if let Some(parent) = parent_stack.last() {
                            current_index = Some(tree.add_child(None, *parent));
                            // eprintln!("\t\tSIBL: Adding Child node of parent {parent} -> index is set to {current_index:?}");
                        } else {
                            unreachable!("Sould not be possible to have named child with no parent")
                        };
                        tree.get_mut(current_index.unwrap())
                    };

                    node.name = current_name;
                    if let Some(length) = current_length {
                        node.length = Some(length.parse()?);
                    }

                    // eprintln!("\t\tSIBL: Resetting name and length and index");
                    current_name = None;
                    current_length = None;
                    current_index = None;

                    parsing = Field::Name;
                }
                ')' => {
                    open_delimiters.pop();
                    // eprintln!("\t\tCLOSE: Adding name and length to current node {current_index:?}");
                    let node = if let Some(index) = current_index {
                        tree.get_mut(index)
                    } else {
                        if let Some(parent) = parent_stack.last() {
                            current_index = Some(tree.add_child(None, *parent));
                            // eprintln!("\t\tCLOSE: Adding Child node of parent {parent} -> index is set to {current_index:?}");
                        } else {
                            unreachable!("Sould not be possible to have named child with no parent")
                        };
                        tree.get_mut(current_index.unwrap())
                    };

                    node.name = current_name;
                    if let Some(length) = current_length {
                        node.length = Some(length.parse()?);
                    }

                    // eprintln!("\t\tCLOSE: Resetting name and index");
                    current_name = None;
                    current_length = None;

                    parsing = Field::Name;

                    // eprintln!(
                    //     "\t\tCLOSE: Setting current index to last parent: {:?}",
                    //     parent_stack.last()
                    // );
                    if let Some(parent) = parent_stack.pop() {
                        current_index = Some(parent)
                    } else {
                        return Err("Ending SubTree with no parent node...".into());
                    }
                }
                ';' => {
                    if !open_delimiters.is_empty() {
                        return Err(ParseError::UnclosedBracket.into());
                    }
                    // Finish parsing the Tree
                    // eprintln!("\t\tEND: Adding name and length to current node {current_index:?}");
                    let node = tree.get_mut(current_index.unwrap());
                    node.name = current_name;
                    if let Some(length) = current_length {
                        node.length = Some(length.parse()?);
                    }

                    return Ok(tree);
                }
                _ => {
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
                                return Err(ParseError::WhiteSpaceInNumber.into());
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

            // eprintln!("\tparsing: {parsing:?}");
            // eprintln!("\tcurrent_name: {current_name:?}");
            // eprintln!("\tcurrent_length: {current_length:?}");
            // eprintln!("\tcurrent_index: {current_index:?}");
            // eprintln!("\tparent_stack: {parent_stack:?}");
        }

        Err(ParseError::NoClosingSemicolon.into())
    }

    /// Saves the tree to a newick file
    pub fn to_file(&self, path: &Path) -> Result<()> {
        match fs::write(path, self.to_newick()) {
            Ok(_) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    /// Loads a tree from a newick file
    pub fn from_file(path: &Path) -> Result<Self> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_newick(&newick_string)
    }

    /// returns a preorder traversal iterator
    pub fn iter_preorder(&self) -> PreOrder<'_> {
        PreOrder {
            tree: self,
            indices: vec![0],
        }
    }

    /// returns a postorder traversal iterator
    pub fn iter_postorder(&self) -> Result<PostOrder<'_>> {
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
    fn new(tree: &'a Tree) -> Result<Self> {
        Ok(Self {
            tree,
            traversal: tree.postorder(0)?.into_iter().rev().collect(),
        })
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

/// Genereates a random binary tree of a given size. Branch lengths are uniformly distributed
pub fn generate_tree(n_leaves: usize, brlens: bool, sampler_type: Distr) -> Tree {
    let mut tree = Tree::new(None);
    let mut rng = thread_rng();

    let sampler = Sampler::new(sampler_type);

    let mut next_deq = VecDeque::new();
    next_deq.push_back(0);

    for _ in 0..(n_leaves - 1) {
        let parent_idx = if rng.gen_bool(0.5) {
            next_deq.pop_front()
        } else {
            next_deq.pop_back()
        }
        .unwrap();
        let l1: Option<f32> = if brlens {
            Some(sampler.sample(&mut rng))
        } else {
            None
        };
        let l2: Option<f32> = if brlens {
            Some(sampler.sample(&mut rng))
        } else {
            None
        };
        next_deq.push_back(tree.add_child_with_len(None, parent_idx, l1));
        next_deq.push_back(tree.add_child_with_len(None, parent_idx, l2));
    }

    for (i, idx) in next_deq.iter().enumerate() {
        tree.get_mut(*idx).set_name(format!("Tip_{i}"));
    }

    tree
}

/// Generates a caterpillar tree by adding children to the last node addesd to the tree
/// until we reach the desired numebr of leaves. Branch lengths are uniformly distributed
pub fn generate_caterpillar(n_leaves: usize, brlens: bool) -> Tree {
    let mut tree = Tree::new(None);
    let mut rng = thread_rng();

    let mut parent = 0;
    for i in 1..n_leaves {
        let parent_bkp = parent;
        let l1: Option<f32> = if brlens { Some(rng.gen()) } else { None };
        let l2: Option<f32> = if brlens { Some(rng.gen()) } else { None };
        if i == n_leaves - 1 {
            // Adding tip
            tree.add_child_with_len(Some(&format!("Tip_{i}")), parent, l1);
            tree.add_child_with_len(Some(&format!("Tip_{}", i + 1)), parent, l2);
        } else {
            // Adding parent node
            parent = tree.add_child_with_len(None, parent, l1);
            tree.add_child_with_len(Some(&format!("Tip_{i}")), parent_bkp, l2);
        }
    }

    tree
}

#[derive(Clone)]
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
    /// Is a tip node
    tip: bool,
    /// Distance to root
    distance_to_root: usize,
    ///
    deleted: bool,
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
            tip: true,
            distance_to_root: 0,
            deleted: false,
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
            tip: true,
            distance_to_root: 0,
            deleted: false,
        }
    }

    /// Sets the internal TreeNode name
    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }

    /// Sets this node to internal
    pub fn set_internal(&mut self) {
        self.tip = false;
    }

    /// Sets this node to tip
    pub fn set_tip(&mut self) {
        self.tip = true;
    }

    /// Check if the node is a tip node
    pub fn is_tip(&self) -> bool {
        self.tip
    }

    /// Returns String with node in newick format
    pub fn to_newick(&self) -> String {
        match (&(self.name), &(self.length)) {
            (None, None) => "".into(),
            (None, Some(l)) => format!(":{l}"),
            (Some(n), None) => n.clone(),
            (Some(n), Some(l)) => format!("{n}:{l}"),
        }
    }
}

impl PartialEq for TreeNode {
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

        let lengths_equal = match (self.length, other.length) {
            (None, None) => true,
            (Some(l1), Some(l2)) => (l1 - l2).abs() < f32::EPSILON,
            _ => false,
        };

        self.name == other.name && self.children.len() == other.children.len() && lengths_equal
    }
}

impl Eq for TreeNode {}

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
            "{:?} <I:{}> (L: {:?})[P: {:?}][Root: {:?}] (C: {:?})",
            self.name, self.idx, self.length, self.parent, self.distance_to_root, self.children,
        )
    }
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {

    use super::*;

    /// Generates example tree from the tree traversal wikipedia page
    /// https://en.wikipedia.org/wiki/Tree_traversal#Depth-first_search
    /// The difference is that I is the left child of G since my tree structure
    /// cannot represent a right child only.
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
        let values: Vec<_> = get_values(&(tree.preorder(0).unwrap()), &tree)
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
        let values: Vec<_> = get_values(&(tree.postorder(0).unwrap()), &tree)
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
            .unwrap()
            .flat_map(|node| node.name.clone())
            .collect();
        assert_eq!(values, vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"])
    }

    #[test]
    fn traverse_inorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.inorder(0).unwrap()), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["A", "B", "C", "D", "E", "F", "H", "I", "G"])
    }

    #[test]
    fn traverse_levelorder() {
        let tree = build_simple_tree();
        let values: Vec<_> = get_values(&(tree.levelorder(0).unwrap()), &tree)
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(values, vec!["F", "B", "G", "A", "D", "I", "C", "E", "H"])
    }

    #[test]
    fn prune_tree() {
        let mut tree = build_simple_tree();
        tree.prune(4); // prune D subtree
        let values: Vec<_> = get_values(&(tree.preorder(0).unwrap()), &tree)
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
                    assert!(d_pred.is_some());
                    assert!((d_pred.unwrap() - d).abs() < f32::EPSILON);
                    // assert!(.is_some_and(|x| (x - d).abs() < f32::EPSILON))
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
            let tree = generate_tree(size, false, Distr::Uniform);
            assert_eq!(tree.get_leaves().len(), size);
        }
    }

    #[test]
    fn genera_gamma() {
        use rand::prelude::*;
        let mut rng = thread_rng();
        for size in (0..20).map(|_| rng.gen_range(10..=100)) {
            let tree = generate_tree(size, true, Distr::Gamma);
            let tree2 = generate_tree(size, true, Distr::Uniform);
            assert_eq!(tree.get_leaves().len(), size);
            assert_eq!(tree2.get_leaves().len(), size);
            println!("G: {}", tree.to_newick());
            println!("U: {}", tree2.to_newick())
        }
    }

    #[test]
    fn to_newick() {
        let tree = build_tree_with_lengths();
        assert_eq!("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;", tree.to_newick());
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
        ];
        for newick in newick_strings {
            let tree = Tree::from_newick(newick).unwrap();
            assert_eq!(newick, tree.to_newick());
        }
    }

    #[test]
    fn read_newick_fails() {
        let newick_strings = vec![
            ("((D,E)B,(F,G,C)A;", ParseError::UnclosedBracket),
            ("((D,E)B,(F,G)C)A", ParseError::NoClosingSemicolon),
        ];
        for (newick, _error) in newick_strings {
            let tree = Tree::from_newick(newick);
            assert!(tree.is_err());
            assert!(tree.err().unwrap().is::<ParseError>())
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
                tree.print_debug();
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
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;",
            "(A,B,(C,D)E)F;",
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

            tree.rescale(scale);

            println!("Dealing with tree: {} and scale {}", orig, scale);
            for (n1, n2) in zip(tree.nodes, rescaled.nodes) {
                assert_eq!(n1, n2)
            }
        }
    }

    #[test]
    fn test_mutability() {
        let tree = Tree::from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
        // First computation
        assert_eq!(tree.diameter().unwrap(), 1.1);
        assert_eq!(tree.height().unwrap(), 0.9);
        // Cached value
        eprintln!("{:?}", tree);
        assert_eq!(tree.diameter().unwrap(), 1.1);
        assert_eq!(tree.height().unwrap(), 0.9);
    }

    #[test]
    fn test_unique_tip_names() {
        let test_cases = vec![
            ("(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);",true),
            ("(((i:0.1,j:0.1):0.1,(a:0.1,b:0.1):0.1):0.1,((c:0.1,d:0.1):0.1,((e:0.1,f:0.1):0.1,(g:0.1,h:0.1):0.1):0.1):0.1);", true),
            ("(((((((((,),),),),),),),),);",false),
            ("(((((((((Tip8,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);",false),
        ];

        for (newick, is_unique) in test_cases {
            assert_eq!(
                Tree::from_newick(newick).unwrap().are_tip_names_unique(),
                is_unique,
                "Failed on: {newick}"
            )
        }
    }

    #[test]
    fn test_descendants() {
        let tree = build_simple_tree();
        let descendants_b: Vec<_> = get_values(&tree.get_descendants(1), &tree)
            .into_iter()
            .flatten()
            .sorted()
            .collect();
        let descendants_g: Vec<_> = get_values(&tree.get_descendants(2), &tree)
            .into_iter()
            .flatten()
            .sorted()
            .collect();

        assert_eq!(descendants_b, vec!["A", "C", "D", "E"]);
        assert_eq!(descendants_g, vec!["H", "I"]);
    }

    #[test]
    fn test_compress() {
        let mut tree = Tree::new(Some("root"));
        tree.add_child_with_len(Some("tip_A"), 0, Some(1.0));
        tree.add_child_with_len(Some("in_B"), 0, Some(1.0));
        tree.add_child_with_len(Some("in_C"), 2, Some(1.0));
        tree.add_child_with_len(Some("tip_D"), 3, Some(1.0));

        tree.compress().unwrap();

        assert_eq!(tree.to_newick(), "(tip_A:1,tip_D:3)root;");
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

        for indices in (0..trees.len()).into_iter().combinations(2) {
            let (i0, i1) = (indices[0], indices[1]);
            let t0 = Tree::from_newick(trees[i0]).unwrap();
            let t1 = Tree::from_newick(trees[i1]).unwrap();

            assert!((t0.weighted_robinson_foulds(&t1).unwrap() - rfs[i0][i1]).abs() <= f32::EPSILON)
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
        let rfs: Vec<Vec<f32>> = vec![
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

        for indices in (0..trees.len()).into_iter().combinations(2) {
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
    fn test_rf_unrooted() {
        let ref_s = "(((aaaaaaaaad:0.18749,aaaaaaaaae:0.18749):0.18749,((aaaaaaaaaf:0.18749,(aaaaaaaaag:0.18749,(aaaaaaaaah:0.18749,(aaaaaaaaai:0.18749,aaaaaaaaaj:0.18749):0.18749):0.18749):0.18749):0.18749,(aaaaaaaaak:0.18749,(aaaaaaaaal:0.18749,aaaaaaaaam:0.18749):0.18749):0.18749):0.18749):0.18749,((aaaaaaaaan:0.18749,aaaaaaaaao:0.18749):0.18749,(aaaaaaaaaa:0.18749,(aaaaaaaaab:0.18749,aaaaaaaaac:0.18749):0.18749):0.18749):0.18749);";
        let prd_s = "(aaaaaaaaag:0.24068,(aaaaaaaaah:0.21046,(aaaaaaaaai:0.15487,aaaaaaaaaj:0.17073)1.000:0.22813)0.999:0.26655,(aaaaaaaaaf:0.27459,((((aaaaaaaaan:0.17964,aaaaaaaaao:0.13686)0.994:0.18171,(aaaaaaaaaa:0.19386,(aaaaaaaaab:0.15663,aaaaaaaaac:0.20015)1.000:0.26799)0.981:0.15442)0.999:0.38320,(aaaaaaaaad:0.18133,aaaaaaaaae:0.17164)0.990:0.18734)0.994:0.18560,(aaaaaaaaak:0.24485,(aaaaaaaaal:0.17930,aaaaaaaaam:0.22072)1.000:0.22274)0.307:0.05569)1.000:0.22736)0.945:0.12401);";

        let ref_tree = Tree::from_newick(ref_s).unwrap();
        let prd_tree = Tree::from_newick(prd_s).unwrap();

        let ref_parts: HashSet<_> = HashSet::from_iter(ref_tree.get_partitions().unwrap());
        let prd_parts: HashSet<_> = HashSet::from_iter(prd_tree.get_partitions().unwrap());

        // let common = ref_parts.intersection(&prd_parts).count();

        println!("Leaf indices:");
        println!("Ref: {:#?}", ref_tree.leaf_index);

        println!("\nPartitions: ");
        println!("Ref: ");
        for i in ref_parts.iter().sorted() {
            println!("\t{i:#018b}");
        }
        println!("\nPrd: ");
        for i in prd_parts.iter().sorted() {
            println!("\t{i:#018b}");
        }

        // println!("tree\treference\tcommon\tcompared\trf\trf_comp");
        // println!(
        //     "0\t{}\t{}\t{}\t{}\t{}",
        //     ref_parts.len() - common,
        //     common,
        //     prd_parts.len() - common,
        //     rf,
        //     ref_parts.len() + prd_parts.len() - 2*common,
        // );

        // panic!()
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
                t1.robinson_foulds_new(&t2).unwrap()
            )
        }
    }

    #[test]
    fn medium() {
        fn get_bitset(hash: usize, len: usize) -> FixedBitSet {
            // eprintln!("Converting : {hash}");
            let mut set = FixedBitSet::with_capacity(len);
            let mut hash = hash;
            for i in 0..len {
                // eprintln!("\t{hash:#0len$b}", len = len);
                if hash & 1 == 1 {
                    set.insert(i)
                }
                hash >>= 1
            }
            let mut toggled = set.clone();
            toggled.toggle_range(..);

            // println!("Done\n");

            set.min(toggled)
        }

        let n1 = "((((Tip_13,Tip_14),(Tip_15,(Tip_16,Tip_17))),(((Tip_18,Tip_19),Tip_0),(Tip_1,Tip_2))),((Tip_3,(Tip_4,Tip_5)),(Tip_6,(Tip_7,(Tip_8,(Tip_9,(Tip_10,(Tip_11,Tip_12))))))));";
        let n2 = "(((Tip_7,(Tip_8,Tip_9)),((Tip_10,(Tip_11,Tip_12)),((Tip_13,(Tip_14,Tip_15)),(Tip_16,Tip_17)))),((Tip_18,Tip_19),(Tip_0,(Tip_1,(Tip_2,(Tip_3,(Tip_4,(Tip_5,Tip_6))))))));";
        let rf_true = 26;

        let reftree = Tree::from_newick(n1).unwrap();
        let compare = Tree::from_newick(n2).unwrap();

        let index: Vec<_> = reftree.get_leaf_names().into_iter().sorted().collect();
        let index2: Vec<_> = compare.get_leaf_names().into_iter().sorted().collect();

        let p1 = reftree.get_partitions_new().unwrap();
        println!("REF: [");
        for p in p1 {
            print!("(");
            for b in p.ones() {
                print!("{:?}, ", index[b]);
            }
            println!("),");
        }
        println!("]\n");

        let p2 = compare.get_partitions_new().unwrap();
        println!("COMP: [");
        for p in p2 {
            print!("(");
            for b in p.ones() {
                print!("{:?}, ", index2[b]);
            }
            println!("),");
        }
        println!("]\n");

        for node in compare.nodes.iter() {
            if node.is_tip() {
                println!("COMP TIP: {node:?}")
            }
        }

        // dbg!(&reftree);
        // dbg!(&compare);

        assert_eq!(reftree.robinson_foulds_new(&compare).unwrap(), rf_true)
    }

    #[test]
    // Robinson foulds distances according to
    // https://evolution.genetics.washington.edu/phylip/doc/treedist.html
    fn robinson_foulds_treedist_new() {
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

            assert_eq!(t0.robinson_foulds_new(&t1).unwrap(), rfs[i0][i1])
        }
    }

    // the reference distance matrix was computed with ete3
    #[test]
    fn compute_distance_matrix() {
        let tree = Tree::from_newick("((A:0.1,B:0.2)F:0.6,(C:0.3,D:0.4)E:0.5)G;").unwrap();
        let true_dists: HashMap<(String, String), f32> = HashMap::from_iter(vec![
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
                (dist - matrix.get(&n1, &n2).unwrap()) <= f32::EPSILON,
                "d({n1},{n2}) want:{dist} got:{}",
                matrix.get(&n1, &n2).unwrap()
            )
        }
    }
}
