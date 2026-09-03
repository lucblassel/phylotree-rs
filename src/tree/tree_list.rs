use std::io::{Read, Write};
use std::iter::FromIterator;
use std::ops::{Index, IndexMut};
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

use super::newick::{NewickParser, NewickTokenizer};
use super::{NewickParseError, NewickParseOptions, Tree, TreeError};

/// An ordered, in-memory collection of phylogenetic trees.
///
/// A `TreeList` is backed by a [`Vec<Tree>`] and preserves insertion order.
#[derive(Debug, Clone, Default)]
pub struct TreeList {
    trees: Vec<Tree>,
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
        Self { trees }
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
        self.trees.get_mut(index)
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
        self.trees.push(tree);
    }

    /// Removes and returns the last tree, or `None` if the list is empty.
    pub fn pop(&mut self) -> Option<Tree> {
        self.trees.pop()
    }

    /// Removes all trees from the list.
    pub fn clear(&mut self) {
        self.trees.clear();
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
        self.trees.iter_mut()
    }

    /// Returns the trees as a slice.
    pub fn as_slice(&self) -> &[Tree] {
        self.trees.as_slice()
    }

    /// Returns the trees as a mutable slice.
    pub fn as_mut_slice(&mut self) -> &mut [Tree] {
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
                None => return Ok(Self { trees }),
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
}
