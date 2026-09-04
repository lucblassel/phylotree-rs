use std::collections::HashMap;
use std::io::{Read, Write};
use std::iter::FromIterator;
use std::ops::{Index, IndexMut};
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

use fixedbitset::FixedBitSet;
use thiserror::Error;

use super::bipartition::BipartitionProfile;
use super::newick::{NewickParser, NewickTokenizer};
use super::taxon_index::TaxonIndex;
use super::{NewickParseError, NewickParseOptions, Tree, TreeError};

/// A canonical nontrivial split of the taxa in a [`TaxonIndex`].
///
/// Equality and hashing use the taxon-to-bit mapping and the underlying bitset.
/// Taxon labels are only accessed when presenting the partition to a user.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bipartition {
    taxon_index: TaxonIndex,
    bits: FixedBitSet,
}

impl Bipartition {
    pub(crate) fn new(taxon_index: TaxonIndex, bits: FixedBitSet) -> Self {
        Self { taxon_index, bits }
    }

    /// Returns the index that defines this partition's bit positions.
    pub fn taxon_index(&self) -> &TaxonIndex {
        &self.taxon_index
    }

    /// Returns whether another bipartition uses the same taxon-to-bit mapping.
    pub fn is_compatible(&self, other: &Self) -> bool {
        self.taxon_index.is_compatible(&other.taxon_index)
    }

    /// Returns the number of taxa on the stored canonical side of the split.
    pub fn len(&self) -> usize {
        self.bits.count_ones(..)
    }

    /// Returns `true` if the stored side of the split contains no taxa.
    pub fn is_empty(&self) -> bool {
        self.bits.is_clear()
    }

    /// Iterates over taxon labels on the stored canonical side of the split.
    pub fn taxa(&self) -> impl Iterator<Item = &str> {
        self.bits
            .ones()
            .map(|index| self.taxon_index.taxon(index).unwrap())
    }

    /// Iterates over taxon labels on the other side of the split.
    pub fn complement_taxa(&self) -> impl Iterator<Item = &str> {
        self.taxon_index
            .taxa()
            .iter()
            .enumerate()
            .filter(|(index, _)| !self.bits.contains(*index))
            .map(|(_, taxon)| taxon.as_str())
    }
}

/// The observed frequency of a distinct tree topology.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopologyFrequency {
    /// Index of the first tree observed with this topology.
    pub representative: usize,
    /// Number of trees with this topology.
    pub count: usize,
    /// Proportion of trees with this topology.
    pub frequency: f64,
}

/// Errors produced by operations involving multiple trees.
#[derive(Debug, Error)]
pub enum TreeListError {
    /// A tree cannot participate in the requested operation.
    #[error("tree {index} is invalid: {source}")]
    InvalidTree {
        /// Index of the invalid tree.
        index: usize,
        /// Underlying tree error.
        #[source]
        source: TreeError,
    },
    /// A tree has a different taxon set from the first tree in the list.
    #[error("tree {index} has a different taxon set")]
    DifferentTaxa {
        /// Index of the tree with different taxa.
        index: usize,
    },
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct TopologyKey {
    rooted: bool,
    partitions: Vec<FixedBitSet>,
    root_partitions: Vec<FixedBitSet>,
}

/// An ordered, in-memory collection of phylogenetic trees.
///
/// A `TreeList` is backed by a [`Vec<Tree>`] and preserves insertion order.
#[derive(Debug, Clone, Default)]
pub struct TreeList {
    trees: Vec<Tree>,
    taxon_index: Option<TaxonIndex>,
}

impl TreeList {
    /// Creates an empty tree list.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates an empty tree list with space for at least `capacity` trees.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            trees: Vec::with_capacity(capacity),
            ..Self::default()
        }
    }

    /// Creates a tree list from a vector of trees.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::{Tree, TreeList};
    ///
    /// let trees = vec![
    ///     Tree::from_newick("(A,B);").unwrap(),
    ///     Tree::from_newick("((A,B),C);").unwrap(),
    /// ];
    /// let tree_list = TreeList::from_trees(trees);
    ///
    /// assert_eq!(tree_list.len(), 2);
    /// assert_eq!(tree_list[1].n_leaves(), 3);
    /// ```
    pub fn from_trees(trees: Vec<Tree>) -> Self {
        Self {
            trees,
            ..Self::default()
        }
    }

    fn invalidate_taxon_index(&mut self) {
        self.taxon_index = None;
    }

    fn invalid_tree(index: usize, source: TreeError) -> TreeListError {
        TreeListError::InvalidTree { index, source }
    }

    fn compute_taxon_index(&self) -> Result<TaxonIndex, TreeListError> {
        let taxon_index = match self.trees.first() {
            Some(tree) => tree
                .taxon_index()
                .map_err(|error| Self::invalid_tree(0, error))?,
            None => TaxonIndex::empty(),
        };
        for (index, tree) in self.trees.iter().enumerate().skip(1) {
            let current = tree
                .taxon_index()
                .map_err(|error| Self::invalid_tree(index, error))?;
            if !taxon_index.is_compatible(&current) {
                return Err(TreeListError::DifferentTaxa { index });
            }
        }
        Ok(taxon_index)
    }

    /// Returns the stored list-wide taxon index, if it is currently available.
    pub fn cached_taxon_index(&self) -> Option<&TaxonIndex> {
        self.taxon_index.as_ref()
    }

    /// Returns a compatible taxon index after validating every tree's taxa.
    ///
    /// The empty list has an empty index. If no index has been stored, this
    /// method computes one without mutating the list or its trees.
    pub fn taxon_index(&self) -> Result<TaxonIndex, TreeListError> {
        self.taxon_index
            .clone()
            .map(Ok)
            .unwrap_or_else(|| self.compute_taxon_index())
    }

    /// Recomputes and stores one shared taxon index on the list and its trees.
    pub fn rebuild_taxon_index(&mut self) -> Result<&TaxonIndex, TreeListError> {
        self.taxon_index = None;
        let taxon_index = self.compute_taxon_index()?;
        for (index, tree) in self.trees.iter_mut().enumerate() {
            tree.bind_taxon_index(taxon_index.clone())
                .map_err(|error| Self::invalid_tree(index, error))?;
        }
        self.taxon_index = Some(taxon_index);
        Ok(self.taxon_index.as_ref().unwrap())
    }

    /// Validates that every tree has the same uniquely named set of taxa.
    pub fn validate_same_taxa(&self) -> Result<(), TreeListError> {
        self.taxon_index().map(|_| ())
    }

    /// Counts how many trees contain each nontrivial bipartition.
    pub fn partition_counts(&self) -> Result<HashMap<Bipartition, usize>, TreeListError> {
        let taxon_index = self.taxon_index()?;
        let mut counts: HashMap<FixedBitSet, usize> = HashMap::new();

        for (index, tree) in self.trees.iter().enumerate() {
            let profile = BipartitionProfile::from_tree(tree)
                .map_err(|error| Self::invalid_tree(index, error))?;
            for bits in profile.partition_bits() {
                *counts.entry(bits.clone()).or_insert(0) += 1;
            }
        }
        Ok(counts
            .into_iter()
            .map(|(bits, count)| (Bipartition::new(taxon_index.clone(), bits), count))
            .collect())
    }

    /// Returns the proportion of trees containing each nontrivial bipartition.
    pub fn partition_frequencies(&self) -> Result<HashMap<Bipartition, f64>, TreeListError> {
        if self.is_empty() {
            self.validate_same_taxa()?;
            return Ok(HashMap::new());
        }

        let denominator = self.len() as f64;
        Ok(self
            .partition_counts()?
            .into_iter()
            .map(|(partition, count)| (partition, count as f64 / denominator))
            .collect())
    }

    fn topology_key(profile: &BipartitionProfile) -> TopologyKey {
        let mut partitions: Vec<_> = profile.partition_bits().cloned().collect();
        partitions.sort_unstable();

        let mut root_partitions = if profile.is_rooted() {
            profile.root_partition_bits().cloned().collect()
        } else {
            Vec::new()
        };
        root_partitions.sort_unstable();

        TopologyKey {
            rooted: profile.is_rooted(),
            partitions,
            root_partitions,
        }
    }

    /// Returns one frequency record per distinct topology, most frequent first.
    ///
    /// Branch lengths, node names, comments, support values, and child order do
    /// not affect topology identity. Root placement is significant for rooted trees.
    pub fn topology_frequencies(&self) -> Result<Vec<TopologyFrequency>, TreeListError> {
        self.validate_same_taxa()?;
        let mut counts: HashMap<TopologyKey, (usize, usize)> = HashMap::new();

        for (index, tree) in self.trees.iter().enumerate() {
            let profile = BipartitionProfile::from_tree(tree)
                .map_err(|error| Self::invalid_tree(index, error))?;
            let key = Self::topology_key(&profile);
            counts
                .entry(key)
                .and_modify(|(_, count)| *count += 1)
                .or_insert((index, 1));
        }

        let denominator = self.len() as f64;
        let mut frequencies: Vec<_> = counts
            .into_values()
            .map(|(representative, count)| TopologyFrequency {
                representative,
                count,
                frequency: count as f64 / denominator,
            })
            .collect();
        frequencies.sort_by(|left, right| {
            right
                .count
                .cmp(&left.count)
                .then_with(|| left.representative.cmp(&right.representative))
        });
        Ok(frequencies)
    }

    /// Returns the most frequent topology, or `None` for an empty list.
    ///
    /// Ties are resolved in favor of the topology encountered first.
    pub fn most_frequent_topology(&self) -> Result<Option<TopologyFrequency>, TreeListError> {
        Ok(self.topology_frequencies()?.into_iter().next())
    }

    /// Returns the number of distinct topologies in the list.
    pub fn n_topologies(&self) -> Result<usize, TreeListError> {
        Ok(self.topology_frequencies()?.len())
    }

    /// Returns the number of trees in the list.
    pub fn len(&self) -> usize {
        self.trees.len()
    }

    /// Returns `true` if the list contains no trees.
    pub fn is_empty(&self) -> bool {
        self.trees.is_empty()
    }

    /// Returns the number of trees the list can hold without reallocating.
    pub fn capacity(&self) -> usize {
        self.trees.capacity()
    }

    /// Returns a reference to the tree at `index`, or `None` if it is out of bounds.
    pub fn get(&self, index: usize) -> Option<&Tree> {
        self.trees.get(index)
    }

    /// Returns a mutable reference to the tree at `index`, or `None` if it is out of bounds.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Tree> {
        if index < self.trees.len() {
            self.invalidate_taxon_index();
        }
        let tree = self.trees.get_mut(index)?;
        tree.reset_bipartition_cache();
        Some(tree)
    }

    /// Returns the first tree, or `None` if the list is empty.
    pub fn first(&self) -> Option<&Tree> {
        self.trees.first()
    }

    /// Returns the last tree, or `None` if the list is empty.
    pub fn last(&self) -> Option<&Tree> {
        self.trees.last()
    }

    /// Appends a tree to the back of the list.
    pub fn push(&mut self, tree: Tree) {
        self.invalidate_taxon_index();
        self.trees.push(tree);
    }

    /// Removes and returns the last tree, or `None` if the list is empty.
    pub fn pop(&mut self) -> Option<Tree> {
        let tree = self.trees.pop();
        if tree.is_some() {
            self.invalidate_taxon_index();
        }
        tree
    }

    /// Removes all trees from the list.
    pub fn clear(&mut self) {
        if !self.trees.is_empty() {
            self.trees.clear();
            self.invalidate_taxon_index();
        }
    }

    /// Returns an iterator over the trees.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::{Tree, TreeList};
    ///
    /// let trees = TreeList::from_trees(vec![
    ///     Tree::from_newick("(A,B);").unwrap(),
    ///     Tree::from_newick("((A,B),C);").unwrap(),
    /// ]);
    /// let leaf_counts: Vec<_> = trees.iter().map(Tree::n_leaves).collect();
    ///
    /// assert_eq!(leaf_counts, vec![2, 3]);
    /// ```
    pub fn iter(&self) -> std::slice::Iter<'_, Tree> {
        self.trees.iter()
    }

    /// Returns a mutable iterator over the trees.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::TreeList;
    ///
    /// let mut trees = TreeList::from_newick("(A:1,B:2); (A:2,B:3);").unwrap();
    /// for tree in trees.iter_mut() {
    ///     tree.rescale(2.0);
    /// }
    ///
    /// assert_eq!(trees[0].length().unwrap(), 6.0);
    /// assert_eq!(trees[1].length().unwrap(), 10.0);
    /// ```
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Tree> {
        self.invalidate_taxon_index();
        for tree in &mut self.trees {
            tree.reset_bipartition_cache();
        }
        self.trees.iter_mut()
    }

    /// Returns the trees as a slice.
    pub fn as_slice(&self) -> &[Tree] {
        self.trees.as_slice()
    }

    /// Returns the trees as a mutable slice.
    pub fn as_mut_slice(&mut self) -> &mut [Tree] {
        self.invalidate_taxon_index();
        for tree in &mut self.trees {
            tree.reset_bipartition_cache();
        }
        self.trees.as_mut_slice()
    }

    /// Consumes the list and returns its backing vector.
    pub fn into_trees(self) -> Vec<Tree> {
        self.trees
    }

    /// Parses zero or more consecutive Newick trees from a string.
    ///
    /// Trees may be separated by arbitrary whitespace. An empty or whitespace-only
    /// string produces an empty list.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::TreeList;
    ///
    /// let trees = TreeList::from_newick("(A,B);\n((A,B),C);").unwrap();
    ///
    /// assert_eq!(trees.len(), 2);
    /// assert_eq!(trees[0].to_newick().unwrap(), "(A,B);");
    /// assert_eq!(trees[1].n_leaves(), 3);
    /// ```
    pub fn from_newick(newick: &str) -> Result<Self, NewickParseError> {
        Self::from_reader(newick.as_bytes())
    }

    /// Parses zero or more consecutive Newick trees from a string using explicit options.
    ///
    /// The same options are applied to every tree in the input.
    ///
    /// # Example
    /// ```
    /// use phylotree::tree::{
    ///     InternalNodeLabelMode, NewickParseOptions, TreeList,
    /// };
    ///
    /// let options = NewickParseOptions {
    ///     internal_node_labels: InternalNodeLabelMode::NumericSupport,
    ///     ..Default::default()
    /// };
    /// let trees = TreeList::from_newick_with_options(
    ///     "(A,B)90; (A,B)95;",
    ///     options,
    /// ).unwrap();
    ///
    /// let supports: Vec<_> = trees.iter().map(|tree| {
    ///     let root = tree.get_root().unwrap();
    ///     tree.get(&root).unwrap().support
    /// }).collect();
    /// assert_eq!(supports, vec![Some(90.0), Some(95.0)]);
    /// ```
    pub fn from_newick_with_options(
        newick: &str,
        options: NewickParseOptions,
    ) -> Result<Self, NewickParseError> {
        Self::from_reader_with_options(newick.as_bytes(), options)
    }

    /// Parses zero or more consecutive Newick trees from a reader.
    ///
    /// Input is consumed incrementally and does not need to be loaded into a string first.
    ///
    /// # Example
    /// ```
    /// use std::io::Cursor;
    /// use phylotree::tree::TreeList;
    ///
    /// let reader = Cursor::new(b"(A,B); ((A,B),C);");
    /// let trees = TreeList::from_reader(reader).unwrap();
    ///
    /// assert_eq!(trees.iter().map(|tree| tree.n_leaves()).collect::<Vec<_>>(), vec![2, 3]);
    /// ```
    pub fn from_reader<R: Read>(reader: R) -> Result<Self, NewickParseError> {
        Self::from_reader_with_options(reader, NewickParseOptions::default())
    }

    /// Parses zero or more consecutive Newick trees from a reader using explicit options.
    ///
    /// # Example
    /// ```
    /// use std::io::Cursor;
    /// use phylotree::tree::{
    ///     InternalNodeLabelMode, NewickParseOptions, TreeList,
    /// };
    ///
    /// let reader = Cursor::new(b"(A,B)90; (A,B)95;");
    /// let options = NewickParseOptions {
    ///     internal_node_labels: InternalNodeLabelMode::NumericSupport,
    ///     ..Default::default()
    /// };
    /// let trees = TreeList::from_reader_with_options(reader, options).unwrap();
    ///
    /// assert_eq!(trees.len(), 2);
    /// let root = trees[0].get_root().unwrap();
    /// assert_eq!(trees[0].get(&root).unwrap().support, Some(90.0));
    /// ```
    pub fn from_reader_with_options<R: Read>(
        reader: R,
        options: NewickParseOptions,
    ) -> Result<Self, NewickParseError> {
        let mut tokenizer = NewickTokenizer::new(reader);
        let mut trees = Vec::new();

        loop {
            let nodes = NewickParser::with_options(&mut tokenizer, options.clone()).parse_next()?;
            match nodes {
                Some(nodes) => trees.push(Tree::from_nodes(nodes)),
                None => return Ok(Self::from_trees(trees)),
            }
        }
    }

    /// Reads zero or more consecutive Newick trees from a file.
    ///
    /// # Example
    /// ```no_run
    /// use std::path::Path;
    /// use phylotree::tree::TreeList;
    ///
    /// let trees = TreeList::from_file(Path::new("posterior.trees"))?;
    /// println!("read {} trees", trees.len());
    /// # Ok::<(), phylotree::tree::NewickParseError>(())
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &Path) -> Result<Self, NewickParseError> {
        let file = std::fs::File::open(path)?;
        Self::from_reader(file)
    }

    /// Reads zero or more consecutive Newick trees from a file using explicit options.
    ///
    /// # Example
    /// ```no_run
    /// use std::path::Path;
    /// use phylotree::tree::{
    ///     InternalNodeLabelMode, NewickParseOptions, TreeList,
    /// };
    ///
    /// let options = NewickParseOptions {
    ///     internal_node_labels: InternalNodeLabelMode::NumericSupport,
    ///     ..Default::default()
    /// };
    /// let trees = TreeList::from_file_with_options(
    ///     Path::new("posterior.trees"),
    ///     options,
    /// )?;
    /// println!("read {} trees", trees.len());
    /// # Ok::<(), phylotree::tree::NewickParseError>(())
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file_with_options(
        path: &Path,
        options: NewickParseOptions,
    ) -> Result<Self, NewickParseError> {
        let file = std::fs::File::open(path)?;
        Self::from_reader_with_options(file, options)
    }

    /// Serializes the trees as one Newick tree per line.
    pub fn to_newick(&self) -> Result<String, TreeError> {
        let mut output = Vec::new();
        self.write_newick(&mut output)?;
        Ok(String::from_utf8(output).expect("Newick serialization always produces valid UTF-8"))
    }

    /// Writes the trees to a stream as one Newick tree per line.
    pub fn write_newick<W: Write>(&self, mut writer: W) -> Result<(), TreeError> {
        for tree in &self.trees {
            writer.write_all(tree.to_newick()?.as_bytes())?;
            writer.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Writes the trees to a file as one Newick tree per line.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn to_file(&self, path: &Path) -> Result<(), TreeError> {
        let file = std::fs::File::create(path)?;
        self.write_newick(std::io::BufWriter::new(file))
    }
}

impl From<Vec<Tree>> for TreeList {
    fn from(trees: Vec<Tree>) -> Self {
        Self::from_trees(trees)
    }
}

impl From<TreeList> for Vec<Tree> {
    fn from(tree_list: TreeList) -> Self {
        tree_list.into_trees()
    }
}

impl FromIterator<Tree> for TreeList {
    fn from_iter<T: IntoIterator<Item = Tree>>(iter: T) -> Self {
        Self::from_trees(iter.into_iter().collect())
    }
}

impl Extend<Tree> for TreeList {
    fn extend<T: IntoIterator<Item = Tree>>(&mut self, iter: T) {
        self.invalidate_taxon_index();
        self.trees.extend(iter);
    }
}

impl IntoIterator for TreeList {
    type Item = Tree;
    type IntoIter = std::vec::IntoIter<Tree>;

    fn into_iter(self) -> Self::IntoIter {
        self.trees.into_iter()
    }
}

impl<'a> IntoIterator for &'a TreeList {
    type Item = &'a Tree;
    type IntoIter = std::slice::Iter<'a, Tree>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut TreeList {
    type Item = &'a mut Tree;
    type IntoIter = std::slice::IterMut<'a, Tree>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

impl Index<usize> for TreeList {
    type Output = Tree;

    fn index(&self, index: usize) -> &Self::Output {
        &self.trees[index]
    }
}

impl IndexMut<usize> for TreeList {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.invalidate_taxon_index();
        self.trees[index].reset_bipartition_cache();
        &mut self.trees[index]
    }
}

impl AsRef<[Tree]> for TreeList {
    fn as_ref(&self) -> &[Tree] {
        self.as_slice()
    }
}

impl AsMut<[Tree]> for TreeList {
    fn as_mut(&mut self) -> &mut [Tree] {
        self.as_mut_slice()
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::tree::InternalNodeLabelMode;

    #[test]
    fn parses_multiple_trees_from_a_reader() {
        let input = "(A,B);\n((A,C),B);\n('A;1',B);";
        let trees = TreeList::from_reader(Cursor::new(input)).unwrap();

        assert_eq!(trees.len(), 3);
        assert_eq!(trees[0].n_leaves(), 2);
        assert_eq!(trees[1].n_leaves(), 3);
        assert_eq!(trees[2].get_leaf_names()[0].as_deref(), Some("A;1"));
    }

    #[test]
    fn parses_an_empty_stream_as_an_empty_list() {
        let trees = TreeList::from_newick("  \n\t").unwrap();
        assert!(trees.is_empty());
    }

    #[test]
    fn applies_parse_options_to_every_tree() {
        let trees = TreeList::from_newick_with_options(
            "(A,B)90; (A,B)95;",
            NewickParseOptions {
                internal_node_labels: InternalNodeLabelMode::NumericSupport,
                ..Default::default()
            },
        )
        .unwrap();

        let supports: Vec<_> = trees
            .iter()
            .map(|tree| tree.get(&tree.get_root().unwrap()).unwrap().support)
            .collect();
        assert_eq!(supports, vec![Some(90.0), Some(95.0)]);
    }

    #[test]
    fn implements_collection_conversions_and_iteration() {
        let input = vec![
            Tree::from_newick("(A,B);").unwrap(),
            Tree::from_newick("(A,B,C);").unwrap(),
        ];
        let mut trees: TreeList = input.into_iter().collect();

        assert_eq!(
            trees.iter().map(Tree::n_leaves).collect::<Vec<_>>(),
            vec![2, 3]
        );
        for tree in &mut trees {
            tree.rescale(2.0);
        }

        let trees: Vec<Tree> = trees.into_iter().collect();
        assert_eq!(trees.len(), 2);
    }

    #[test]
    fn serializes_one_tree_per_line() {
        let trees = TreeList::from_newick("(A,B); (A,(B,C));").unwrap();
        assert_eq!(trees.to_newick().unwrap(), "(A,B);\n(A,(B,C));\n");
    }

    #[test]
    fn reports_an_incomplete_last_tree() {
        let error = TreeList::from_newick("(A,B); (A,B)").unwrap_err();
        assert!(matches!(error, NewickParseError::UnexpectedEOF));
    }

    #[test]
    fn validates_identical_taxon_sets_independent_of_order() {
        let trees = TreeList::from_newick("((A,B),C); (C,(B,A));").unwrap();
        let taxon_index = trees.taxon_index().unwrap();

        assert_eq!(taxon_index.taxa(), &["A", "B", "C"]);
        assert_eq!(taxon_index.index_of("B"), Some(1));
        assert_eq!(taxon_index.taxon(2), Some("C"));
        assert!(trees.validate_same_taxa().is_ok());
    }

    #[test]
    fn rejects_different_or_invalid_taxon_sets() {
        let different = TreeList::from_newick("((A,B),C); ((A,B),D);").unwrap();
        assert!(matches!(
            different.validate_same_taxa(),
            Err(TreeListError::DifferentTaxa { index: 1 })
        ));

        let unnamed = TreeList::from_newick("((A,B),C); ((A,B),);").unwrap();
        assert!(matches!(
            unnamed.validate_same_taxa(),
            Err(TreeListError::InvalidTree {
                index: 1,
                source: TreeError::UnnamedLeaves
            })
        ));

        let duplicate = TreeList::from_newick("((A,B),C); ((A,A),C);").unwrap();
        assert!(matches!(
            duplicate.validate_same_taxa(),
            Err(TreeListError::InvalidTree {
                index: 1,
                source: TreeError::DuplicateLeafNames
            })
        ));
    }

    #[test]
    fn counts_and_decodes_bipartitions() {
        let trees = TreeList::from_newick("((A,B),(C,D)); ((B,A),(D,C)); ((A,C),(B,D));").unwrap();
        let counts = trees.partition_counts().unwrap();
        let frequencies = trees.partition_frequencies().unwrap();

        assert_eq!(counts.len(), 2);
        let ab = counts
            .iter()
            .find(|(partition, _)| partition.taxa().collect::<Vec<_>>() == vec!["A", "B"])
            .unwrap();
        assert_eq!(*ab.1, 2);
        assert_eq!(frequencies[ab.0], 2.0 / 3.0);
        assert_eq!(ab.0.complement_taxa().collect::<Vec<_>>(), vec!["C", "D"]);
    }

    #[test]
    fn mutable_access_invalidates_taxon_and_partition_caches() {
        let mut trees = TreeList::from_newick("((A,B),C); ((A,B),C);").unwrap();
        trees.partition_counts().unwrap();

        trees[1].get_by_name_mut("C").unwrap().name = Some("D".to_owned());

        assert!(matches!(
            trees.validate_same_taxa(),
            Err(TreeListError::DifferentTaxa { index: 1 })
        ));
    }

    #[test]
    fn list_taxon_index_rebuilding_is_explicit() {
        let mut trees = TreeList::from_newick("((A,B),C); (C,(B,A));").unwrap();
        let root = trees[0].get_root().unwrap();
        trees[0].get_mut(&root).unwrap().support = Some(1.0);
        assert!(trees.cached_taxon_index().is_none());
        assert!(trees[0].cached_taxon_index().is_none());

        trees.taxon_index().unwrap();
        assert!(trees.cached_taxon_index().is_none());
        assert!(trees[0].cached_taxon_index().is_none());

        trees.rebuild_taxon_index().unwrap();
        assert!(trees.cached_taxon_index().is_some());
        assert!(trees.iter().all(|tree| tree.cached_taxon_index().is_some()));
    }

    #[test]
    fn independently_created_taxon_indices_are_compatible() {
        let trees = TreeList::from_newick("((A,B),(C,D));").unwrap();
        let clone = trees.clone();
        let independent = TreeList::from_newick("((B,A),(D,C));").unwrap();

        let index = trees.taxon_index().unwrap();
        assert_eq!(index, clone.taxon_index().unwrap());
        assert!(index.is_compatible(&independent.taxon_index().unwrap()));

        let first = trees[0].bipartitions().unwrap();
        let second = independent[0].bipartitions().unwrap();
        assert_eq!(first, second);
        assert!(first
            .iter()
            .next()
            .unwrap()
            .is_compatible(second.iter().next().unwrap()));
        assert_eq!(trees[0].robinson_foulds(&independent[0]).unwrap(), 0);
    }

    #[test]
    fn groups_topologies_independent_of_child_order() {
        let trees = TreeList::from_newick("((A,B),(C,D)); ((D,C),(B,A)); ((A,C),(B,D));").unwrap();
        let frequencies = trees.topology_frequencies().unwrap();

        assert_eq!(frequencies.len(), 2);
        assert_eq!(frequencies[0].representative, 0);
        assert_eq!(frequencies[0].count, 2);
        assert_eq!(frequencies[0].frequency, 2.0 / 3.0);
        assert_eq!(
            trees.most_frequent_topology().unwrap(),
            Some(frequencies[0])
        );
        assert_eq!(trees.n_topologies().unwrap(), 2);
    }

    #[test]
    fn rooted_topology_includes_root_placement() {
        let trees = TreeList::from_newick("((A,B),(C,D)); (A,(B,(C,D)));").unwrap();
        let frequencies = trees.topology_frequencies().unwrap();

        assert_eq!(frequencies.len(), 2);
        assert_eq!(frequencies[0].count, 1);
        assert_eq!(frequencies[1].count, 1);
        assert_eq!(frequencies[0].representative, 0);
    }

    #[test]
    fn empty_list_has_empty_frequency_results() {
        let trees = TreeList::new();

        assert!(trees.validate_same_taxa().is_ok());
        assert!(trees.partition_counts().unwrap().is_empty());
        assert!(trees.partition_frequencies().unwrap().is_empty());
        assert!(trees.topology_frequencies().unwrap().is_empty());
        assert_eq!(trees.most_frequent_topology().unwrap(), None);
        assert_eq!(trees.n_topologies().unwrap(), 0);
    }
}
