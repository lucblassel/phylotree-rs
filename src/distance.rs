//! Compute and manipulate phylogenetic distance matrices
//!

use std::{
    collections::{HashMap, HashSet},
    fmt::{Debug, Display},
    fs,
    path::Path,
    str::FromStr,
};

use itertools::Itertools;
use num_traits::{zero, Float, Zero};
use thiserror::Error;
use trait_set::trait_set;

use crate::tree::{Node, Tree};

trait_set! {
    /// Trait describing objects that can be used as branch lengths
    /// in phylogenetic trees.
    pub trait PairwiseDist = Display + Debug + Float + Zero + FromStr + Clone + Copy;
}

/// Errors that can occur when reading, writing and manipulating [`DistanceMatrix`] structs.
#[derive(Error, Debug)]
pub enum MatrixError {
    /// We are trying to insert a distance that already exists in the matrix
    #[error("Distance already exists, cannot overwrite it.")]
    OverwritingNotPermitted,
    /// We are trying to add more taxa than we allocated space in the matrix for
    #[error("You added more taxa than there is room for in the matrix.")]
    SizeExceeded,
    /// We are trying to access a distance between two taxa that does not
    /// exist in the matrix
    #[error("Missing distance between {0} and {1}.")]
    MissingDistance(String, String),
    /// There was an [`std::io::Error`] when writing the matrix to a phylip file
    #[error("Error writing file")]
    IoError(#[from] std::io::Error),
    /// We are trying to access a taxon that does not exist
    #[error("Missing taxon {0}")]
    MissingTaxon(String),
    /// We are trying to get the pair index for the same leaf
    #[error("Pair index only exists for pairs of different leaves")]
    IndexError,
    /// We are trying to set a non zero distance for an identical taxa pair
    #[error("Identical taxa cannot have a non zero distance")]
    NonZeroIdenticalDistance,
    /// We are trying to add a different number of taxa than what we alloted
    #[error("Trying to add {n_taxa} taxa to a matrix of size {size}")]
    SizeError {
        /// Size of the distance matrix
        size: usize,
        /// Number of taxa we are trying to add
        n_taxa: usize,
    },
}

/// Errors that can occur when parsing phylip distance matrix files.
#[derive(Error, Debug)]
pub enum PhylipParseError<T>
where
    T: Debug,
{
    /// The phylip file is empty
    #[error("Matrix file is empty.")]
    EmptyMatrixFile,
    /// There was a [`std::num::ParseIntError`] when reading the number of taxa
    #[error("Could not parse size from file.")]
    SizeParseError(#[from] std::num::ParseIntError),
    /// One of the matrix rows is empty
    #[error("Row {0} is empty.")]
    EmptyRow(usize),
    /// There was an error when reading a distance.
    #[error("Could not parse distance from file.")]
    DistParseError,
    /// There is a missing distance from one of the matrix rows
    #[error("Missing distance from matrix row {0}")]
    MissingDistance(usize),
    /// The size of the matrix and the number of rows do not match
    #[error("Size and number of rows do not match: {0} rows for size {1}")]
    SizeAndRowsMismatch(usize, usize),
    /// A value on the diagonal is non 0
    #[error("Non 0 diagonal distance for taxa {0}")]
    NonZeroDiagonalValue(String),
    /// The square phylip matrix is not symmetric
    #[error("Non symetric matrix: {0} and {1} are different")]
    NonSymmetric(T, T),
    /// The square phylip matrix is not symmetric
    #[error("Matrix is not symetric.")]
    NonSymmetricMat,
    /// There was a [`MatrixError`] when create the distance matrix object
    #[error("Error creating matrix.")]
    MatrixError(#[from] MatrixError),
    /// There was a [`std::io::Error`] when reading the phylip file
    #[error("Error reading file")]
    IoError(#[from] std::io::Error),
}

#[derive(Debug, Clone)]
/// A phylogenetic distance matrix
pub struct DistanceMatrix<T> {
    /// Number of taxa in the matrix
    pub size: usize,
    /// Identifiers of the taxa
    pub taxa: Vec<String>,
    /// Distances between taxa
    matrix: Vec<T>,
    /// Distance value for identical taxa
    zero: T,
}

impl<T> DistanceMatrix<T>
where
    T: PairwiseDist,
{
    /// Create a new distance matrix for a certain number of sequences
    pub fn new(taxa: Vec<String>, matrix: &[T]) -> Self {
        Self {
            size: taxa.len(),
            taxa,
            matrix: matrix.to_vec(),
            zero: zero(),
        }
    }

    /// Create an empty distance matrix with a given size
    pub fn new_with_size(size: usize) -> Self {
        Self {
            size,
            taxa: Vec::with_capacity(size),
            matrix: vec![zero(); size * (size - 1) / 2],
            zero: zero(),
        }
    }

    /// Create a distance matrix from pre-computed values. The `matrix` parameter
    /// represents the upper triangle of the distance matrix as a single vector.
    /// The matrix index to vector index formula can be found in the [`DistanceMatrix.get_index`]
    /// function.
    pub(crate) fn from_precomputed(taxa: Vec<String>, matrix: Vec<T>) -> Result<Self, MatrixError> {
        // Check that taxa and distances are coherent
        let n = taxa.len();
        let n_pairs = (n * (n - 1)) / 2;
        if matrix.len() != n_pairs {
            return Err(MatrixError::SizeError {
                size: {
                    let delta = (8.0 * (n_pairs as f64) + 1.).sqrt() as usize;
                    (delta + 1) / 2
                },
                n_taxa: n,
            });
        }

        Ok(Self {
            size: taxa.len(),
            taxa,
            matrix,
            zero: zero(),
        })
    }

    /// Set the taxa of the matrix
    pub fn set_taxa(&mut self, taxa: Vec<String>) -> Result<(), MatrixError> {
        if taxa.len() != self.size {
            return Err(MatrixError::SizeError {
                size: self.size,
                n_taxa: taxa.len(),
            });
        }
        self.taxa = taxa;
        Ok(())
    }

    /// Get numerical index associated to taxon identifier
    pub fn get_taxa_index(&self, id: &str) -> Result<usize, MatrixError> {
        self.taxa
            .iter()
            .position(|v| *v == id)
            .ok_or_else(|| MatrixError::MissingTaxon(id.to_string()))
    }

    /// Get the index in the distance vector for 2 sequences
    fn get_pair_index(&self, taxon1: &str, taxon2: &str) -> Result<usize, MatrixError> {
        if taxon1 == taxon2 {
            return Err(MatrixError::IndexError);
        }

        let i = self.get_taxa_index(taxon1)?;
        let j = self.get_taxa_index(taxon2)?;

        self.tril_to_vec_index(i, j)
    }

    // Convert (row,col) matrix index to lower triangular linear index
    fn tril_to_vec_index(&self, i: usize, j: usize) -> Result<usize, MatrixError> {
        if i == j || i >= self.size || j >= self.size {
            Err(MatrixError::IndexError)
        } else {
            Ok(tril_to_rowvec_index(self.size, i, j))
        }
    }

    // Convert lower triangular linear index to (row,col) matrix index
    fn vec_to_tril_index(&self, k: usize) -> Result<(usize, usize), MatrixError> {
        if k >= self.matrix.len() {
            Err(MatrixError::IndexError)
        } else {
            Ok(rowvec_to_tril_index(self.size, k))
        }
    }

    /// Get the distance between two sequences
    pub fn get(&self, id_1: &str, id_2: &str) -> Result<&T, MatrixError> {
        if id_1 == id_2 {
            Ok(&self.zero)
        } else {
            let idx = self.get_pair_index(id_1, id_2)?;
            Ok(&self.matrix[idx])
        }
    }

    /// Set the distance between two sequences
    pub fn set(&mut self, id_1: &str, id_2: &str, dist: T) -> Result<(), MatrixError> {
        if id_1 == id_2 {
            return if dist != self.zero {
                Err(MatrixError::NonZeroIdenticalDistance)
            } else {
                Ok(())
            };
        }

        let idx = self.get_pair_index(id_1, id_2)?;
        self.matrix[idx] = dist;
        Ok(())
    }

    /// Get the distance matrix as a HashMap containing taxa pairs as keys
    /// and pairwise distances as values
    pub fn to_map(&self) -> HashMap<(String, String), T> {
        HashMap::from_iter(self.taxa.iter().cartesian_product(self.taxa.iter()).map(
            |(taxon1, taxon2)| {
                let idx = self.get_pair_index(taxon1, taxon2).unwrap();
                ((taxon1.clone(), taxon2.clone()), self.matrix[idx])
            },
        ))
    }

    /// Outputs a matrix as a phylip formatted string
    pub fn to_phylip(&self, square: bool) -> Result<String, MatrixError> {
        let body = self
            .taxa
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let lim = if square { self.size } else { i };
                let row_s = (0..lim)
                    .map(|j| {
                        let d = if i == j {
                            zero()
                        } else {
                            let idx = self.tril_to_vec_index(i, j).unwrap();
                            self.matrix[idx]
                        };
                        format!("{d}")
                    })
                    .join("  ");
                let mut out = name.clone();
                if !row_s.is_empty() {
                    out += &format!("    {row_s}");
                }
                out
            })
            .join("\n");

        Ok(format!("{}\n{body}\n", self.size))
    }

    /// Writes the matrix to a phylip file
    pub fn to_file(&self, path: &Path, square: bool) -> Result<(), MatrixError> {
        match fs::write(path, self.to_phylip(square)?) {
            Ok(_) => Ok(()),
            Err(e) => Err(MatrixError::IoError(e)),
        }
    }

    fn read_phylip_row(
        row: &str,
        row_num: usize,
        size: usize,
        tril: bool,
    ) -> Result<(String, Vec<T>), PhylipParseError<T>> {
        let mut fields = row.split_whitespace();
        let name = fields.next().ok_or(PhylipParseError::EmptyRow(row_num))?;
        let dists = fields
            .map(|d| d.parse().map_err(|_| PhylipParseError::DistParseError))
            .take(if tril { row_num } else { size })
            .collect::<Result<Vec<_>, _>>()?;

        Ok((name.to_string(), dists))
    }

    /// Build a distance matrix from the lower triangle of a phylip file.
    /// Does not enforce symmetry in the phylip file, it should be faster than
    /// [from_phylip_strict]. Also does not check if all the rows are the same size
    pub fn from_phylip_tril(phylip: &str) -> Result<Self, PhylipParseError<T>> {
        let mut lines = phylip.lines();
        let size: usize = lines
            .next()
            .ok_or(PhylipParseError::EmptyMatrixFile)?
            .parse()
            .map_err(PhylipParseError::SizeParseError)?;

        let mut taxa = vec![];
        let mut max_v = Vec::with_capacity(size * (size - 1) / 2);

        for (i, line) in lines.enumerate() {
            let (name, dists) = Self::read_phylip_row(line, i, size, true)?;

            if dists.len() != i {
                return Err(PhylipParseError::MissingDistance(i + 1));
            }

            taxa.push(name.to_string());
            max_v.extend(&dists[..]);
        }

        if taxa.len() != size {
            return Err(PhylipParseError::SizeAndRowsMismatch(taxa.len(), size));
        }

        Ok(Self::from_precomputed(taxa, max_v)?)
    }

    /// Build a distance matrix from a phylip formatted string, checks that the phylip matrix is
    /// symmetric
    pub fn from_phylip_strict(phylip: &str, square: bool) -> Result<Self, PhylipParseError<T>> {
        let mut lines = phylip.lines();
        let size = lines
            .next()
            .ok_or(PhylipParseError::EmptyMatrixFile)?
            .parse()
            .map_err(PhylipParseError::SizeParseError)?;

        let mut names = vec![];
        let mut rows = vec![];

        for (i, line) in lines.enumerate() {
            let (name, dists) = Self::read_phylip_row(line, i, size, false)?;

            if square && dists.len() != size || !square && dists.len() != i {
                return Err(PhylipParseError::MissingDistance(i + 1));
            }

            if square && dists[i] != zero() {
                return Err(PhylipParseError::NonZeroDiagonalValue(name.to_string()));
            }

            names.push(name);
            rows.push(dists);
        }

        if names.len() != size {
            return Err(PhylipParseError::SizeAndRowsMismatch(names.len(), size));
        }

        let mut matrix = Self::new_with_size(size);
        matrix.set_taxa(names.iter().cloned().map(|v| v.to_string()).collect_vec())?;

        let mut seen = HashSet::new();

        for (n1, row) in names.iter().zip(rows) {
            for (n2, dist) in names.iter().zip(row) {
                if seen.contains(&(n2.to_string(), n1.to_string())) {
                    let known = matrix.get(n1, n2)?;
                    if *known != dist {
                        return Err(PhylipParseError::NonSymmetric(*known, dist));
                    }
                } else {
                    seen.insert((n1.to_string(), n2.to_string()));
                    matrix.set(n1, n2, dist)?;
                }
            }
        }

        Ok(matrix)
    }

    /// Find smallest distance in the distance matrix as well as the
    /// corresponding index
    pub fn min(&self) -> Option<((usize, usize), T)> {
        self.indexed_iter().fold(None, |a, ((i, j), v)| match a {
            None => Some(((i, j), *v)),
            Some((_, av)) => {
                if *v < av {
                    Some(((i, j), *v))
                } else {
                    a
                }
            }
        })
    }

    /// Find largest distance in the distance matrix as well as the
    /// corresponding index
    pub fn max(&self) -> Option<((usize, usize), T)> {
        self.indexed_iter().fold(None, |a, ((i, j), v)| match a {
            None => Some(((i, j), *v)),
            Some((_, av)) => {
                if *v > av {
                    Some(((i, j), *v))
                } else {
                    a
                }
            }
        })
    }

    /// Reads the matrix from a phylip file
    pub fn from_file(path: &Path, square: bool) -> Result<Self, PhylipParseError<T>> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_phylip_strict(&newick_string, square)
    }

    /// Iterator over the lower triangle of the matrix
    pub fn iter(&self) -> impl Iterator<Item = &'_ T> {
        self.matrix.iter()
    }

    /// Iterator over the lower triangle of the matrix with the corresponding
    /// (row,col) matrix indices
    pub fn indexed_iter(&self) -> impl Iterator<Item = ((usize, usize), &'_ T)> {
        self.iter()
            .enumerate()
            .map(|(k, v)| (self.vec_to_tril_index(k).unwrap(), v))
    }
}

impl<T> IntoIterator for DistanceMatrix<T>
where
    T: PairwiseDist,
{
    type Item = T;
    type IntoIter = std::vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.matrix.into_iter()
    }
}

/// Tree building methods
impl<T> DistanceMatrix<T>
where
    T: PairwiseDist,
{
    /// Build a [`Tree`] from the distance matrix using the
    /// [UPGMA](https://en.wikipedia.org/wiki/UPGMA) algorithm
    pub fn upgma(&self) -> Result<Tree, MatrixError> {
        // Setup
        let mut dm = (*self).clone();
        let mut card = vec![1; dm.size];
        //let mut taxa = self.taxa.clone();
        let mut merged = vec![false; dm.size];
        let mut heights = vec![0.; dm.size];
        let mut n_clusters = self.size;

        // Build starting star tree
        let mut tree = Tree::new();
        let virt_root = tree.add(Node::new());
        let mut node_ids = self
            .taxa
            .iter()
            .map(|n| tree.add_child(Node::new_named(n), virt_root, None).unwrap())
            .collect_vec();

        let max_val = T::infinity();

        while n_clusters > 2 {
            let ((a, b), d_ab) = dm.min().ok_or(MatrixError::IndexError)?; // TODO: Add error
            let u = a;
            let node_height = d_ab.to_f64().unwrap() / 2.;
            let d_au = node_height - heights[a];
            let d_bu = node_height - heights[b];

            merged[b] = true;
            heights[u] += d_au;

            // Merge nodes a and b in to new parent u
            let u_node = tree
                .merge_children(
                    &node_ids[a],
                    &node_ids[b],
                    Some(d_au),
                    Some(d_bu),
                    None,
                    None,
                )
                .unwrap();

            // get cardinality of new parent node u
            let c_a = T::from(card[a]).unwrap();
            let c_b = T::from(card[b]).unwrap();
            for (x, _) in merged.iter().enumerate().filter(|(_, &p)| !p) {
                if x != a && b != x {
                    // Set distance from new parent to other unmerged nodes
                    let d_ax = dm.matrix[self.tril_to_vec_index(x, a)?];
                    let d_bx = dm.matrix[self.tril_to_vec_index(x, b)?];
                    let d_ux = (c_a * d_ax + c_b * d_bx) / (c_a + c_b);
                    dm.matrix[self.tril_to_vec_index(u, x)?] = d_ux;
                }

                // Set distance to merge node as max if necessary
                if b != x {
                    let b_idx = dm.tril_to_vec_index(b, x)?;
                    dm.matrix[b_idx] = max_val;
                }
            }

            // Replace a with new node u + update cardinality
            node_ids[u] = u_node;
            card[u] = card[a] + card[b];

            // Update number of clusters
            n_clusters -= 1;
        }

        // Merge last nodes
        let Some((a_i, b_i)) = merged
            .into_iter()
            .enumerate()
            .filter_map(|(i, v)| if v { None } else { Some(i) })
            .next_tuple()
        else {
            return Err(MatrixError::IndexError);
        };

        let (a, b) = (node_ids[a_i], node_ids[b_i]);

        let d_ab = dm.matrix[self.tril_to_vec_index(a_i, b_i)?];
        let tree_height = d_ab.to_f64().unwrap() / 2.;
        let d_ar = tree_height - heights[a_i];
        let d_br = tree_height - heights[b_i];

        let root_node = tree.get_mut(&virt_root).unwrap();
        root_node.set_child_edge(&a, Some(d_ar));
        root_node.set_child_edge(&b, Some(d_br));
        tree.get_mut(&a).unwrap().parent_edge = Some(d_ar);
        tree.get_mut(&b).unwrap().parent_edge = Some(d_br);

        Ok(tree)
    }

    /// Build a [`Tree`] from the distance matrix using the
    /// [neighbor joining](https://en.wikipedia.org/wiki/Neighbor_joining) algorithm
    pub fn neighbor_joining(&self) -> Result<Tree, MatrixError> {
        todo!()
    }
}

////////////////////////
// INDEXING UTILITIES //
////////////////////////

// Convert (row,col) matrix index of triangular matrix to the index
// in the row-wise vector representation of the matrix
pub(crate) fn tril_to_rowvec_index(_size: usize, i: usize, j: usize) -> usize {
    let (i, j) = if i > j { (i, j) } else { (j, i) };

    ((i.saturating_sub(1)) * i) / 2 + j
}

// Convert index in row-wise vector representation of the lower
// triangular matrix to the (row,col) matrix index.
// Formula from [hal-02047514](https://hal.science/hal-02047514)
pub(crate) fn rowvec_to_tril_index(_size: usize, k: usize) -> (usize, usize) {
    let p = (((1.0 + 8.0 * (k as f64)).sqrt() - 1.0) / 2.0).floor() as usize;
    let i = p + 1;
    let j = k - p * (p + 1) / 2;

    (i, j)
}

#[allow(dead_code)]
// Convert (row,col) matrix index of triangular matrix to the index
// in the row-wise vector representation of the matrix
pub(crate) fn tril_to_colvec_index(size: usize, i: usize, j: usize) -> usize {
    let (i, j) = if j < i { (j, i) } else { (i, j) };
    ((2 * size - 3 - i) * i) / 2 + j - 1
}

#[allow(dead_code)]
// Convert index in col-wise vector representation of the lower
// triangular matrix to the (row,col) matrix index.
// Formula from [hal-02047514](https://hal.science/hal-02047514)
pub(crate) fn colvec_to_tril_index(size: usize, k: usize) -> (usize, usize) {
    // k + 1 because formula is for 1-indexed
    let k = k + 1;
    let n = size;

    let kp = n * (n - 1) / 2 - k;
    let p = (((1. + 8. * (kp as f64)).sqrt() - 1.) / 2.).floor() as usize;

    let i = n - (kp - p * (p + 1) / 2) - 1;
    let j = n - 1 - p - 1;

    (i, j)
}

#[cfg(test)]
mod tests {

    use core::panic;

    use ndarray::Array2;

    use super::*;

    const SQUARE: &str = "4
s1    0  2  3  5
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";

    const TRIANGLE: &str = "4
s1
s2    2
s3    3  6
s5    5  10  15
";

    fn build_matrix() -> DistanceMatrix<f32> {
        let names = [("s1", 1.0), ("s2", 2.0), ("s3", 3.0), ("s5", 5.0)];
        let mut matrix = DistanceMatrix::new_with_size(names.len());
        matrix
            .set_taxa(
                names
                    .iter()
                    .cloned()
                    .map(|(n, _)| n.to_string())
                    .collect_vec(),
            )
            .unwrap();

        for pair in names.iter().combinations(2) {
            let (n1, d1) = pair[0];
            let (n2, d2) = pair[1];
            matrix.set(n1, n2, d1 * d2).unwrap();
        }

        matrix
    }

    #[test]
    fn test_to_phylip() {
        let matrix = build_matrix();
        // let matrix = build_matrix_new();

        eprintln!("{:?} ({:?})", matrix.matrix, matrix.taxa);

        assert_eq!(
            SQUARE,
            matrix.to_phylip(true).unwrap(),
            "True:\n{SQUARE}\nPred:\n{}",
            matrix.to_phylip(true).unwrap(),
        );
        assert_eq!(TRIANGLE, matrix.to_phylip(false).unwrap());
    }

    #[test]
    fn from_phylip() -> Result<(), PhylipParseError<f32>> {
        let build: DistanceMatrix<f32> = DistanceMatrix::from_phylip_strict(SQUARE, true)?;
        assert_eq!(
            SQUARE,
            build.to_phylip(true).unwrap(),
            "{SQUARE}\n{}",
            build.to_phylip(true).unwrap()
        );

        let build = DistanceMatrix::from_phylip_strict(TRIANGLE, false)?;
        assert_eq!(
            TRIANGLE,
            build.to_phylip(false).unwrap(),
            "{TRIANGLE}\n{}",
            build.to_phylip(false).unwrap()
        );

        Ok(())
    }

    #[test]
    fn from_phylip_strict_errors() {
        let square_nonsym = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        let matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip_strict(square_nonsym, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::NonSymmetric(_, _) => {}
            PhylipParseError::NonSymmetricMat => {}
            _ => panic!("Error should be 'ParseError::NonSymmetric' not: {err}"),
        }

        let square_missing_dist = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0
s5    5  10  15  0
";
        let mut matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip_strict(square_missing_dist, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::MissingDistance(_) => {}
            _ => panic!("Error should be 'ParseError::MissingDistance' not: {err}"),
        }

        let square_missing_row = "4
s1    0  2  3  7
s2    2  0  6  10
s5    4  6  0  5
";
        matrix = DistanceMatrix::from_phylip_strict(square_missing_row, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::SizeAndRowsMismatch(_, _) => {}
            _ => panic!("Error should be 'ParseError::SizeAndRowsMismatch' not: {err}"),
        }

        let square_missing_size = "s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        matrix = DistanceMatrix::from_phylip_strict(square_missing_size, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::SizeParseError(_) => {}
            _ => panic!("Error should be 'ParseError::SizeParseError' not: {err}"),
        }

        let square_nonzero = "4
s1    1  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    7  10  15  0
";
        let matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip_strict(square_nonzero, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::NonZeroDiagonalValue(_) => {}
            _ => panic!("Error should be 'ParseError::NonZeroDiagonalValue' not: {err}"),
        }
    }

    #[test]
    fn from_phylip_tril_errors() {
        // These errors should be ignored

        let square_nonsym = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        let mut matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip_tril(square_nonsym);
        assert!(matrix.is_ok());

        let square_missing_dist = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0
s5    5  10  15  0
";
        matrix = DistanceMatrix::from_phylip_tril(square_missing_dist);
        assert!(matrix.is_ok());

        let square_nonzero = "4
s1    1  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    7  10  15  0
";
        matrix = DistanceMatrix::from_phylip_tril(square_nonzero);
        assert!(matrix.is_ok());

        // These should still raise errors

        let square_missing_row = "4
s1    0  2  3  7
s2    2  0  6  10
s5    4  6  0  5
";
        matrix = DistanceMatrix::from_phylip_tril(square_missing_row);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::SizeAndRowsMismatch(_, _) => {}
            _ => panic!("Error should be 'ParseError::SizeAndRowsMismatch' not: {err}"),
        }

        let square_missing_size = "s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        matrix = DistanceMatrix::from_phylip_tril(square_missing_size);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            PhylipParseError::SizeParseError(_) => {}
            _ => panic!("Error should be 'ParseError::SizeParseError' not: {err}"),
        }
    }

    #[test]
    fn index_conversions() {
        for size in [5, 10, 20, 50] {
            let mut c = 0;
            for i in 1..size {
                for j in 0..i {
                    let k = tril_to_rowvec_index(size, i, j);
                    assert_eq!(k, c, "Expected {c} got {k}");
                    assert_eq!((i, j), rowvec_to_tril_index(size, k));
                    c += 1;
                }
            }

            c = 0;
            let mut m = Array2::zeros((size, size));
            for j in 0..size {
                for i in (j + 1)..size {
                    *(m.get_mut((i, j)).unwrap()) = c + 1;
                    let k = tril_to_colvec_index(size, i, j);
                    assert_eq!(k, c, "Expected {c} got {k}");
                    assert_eq!((i, j), colvec_to_tril_index(size, k));
                    c += 1;
                }
            }
        }
    }

    #[test]
    fn iterators() {
        let dm = build_matrix();
        let tril = vec![2., 3., 6., 5., 10., 15.];
        let indices = vec![(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)];

        assert!(dm.iter().zip(tril.clone()).all(|(v1, v2)| { *v1 == v2 }));

        assert!(dm
            .indexed_iter()
            .zip(indices.into_iter().zip(tril.clone()))
            .all(|((i1, v1), (i2, v2))| i1 == i2 && *v1 == v2));

        assert!(dm.into_iter().zip(tril).all(|(v1, v2)| { v1 == v2 }));
    }

    #[test]
    fn min_max() {
        let dm = build_matrix();
        let min = Some(((1, 0), 2.));
        let max = Some(((3, 2), 15.));

        assert_eq!(dm.min(), min);
        assert_eq!(dm.max(), max);
    }

    #[test]
    fn build_upgma() {
        // Expected tree
        let expected =
            Tree::from_newick("((c:14.0,d:14.0):2.5,(e:11.0,(a:8.5,b:8.5):2.5):5.5);").unwrap();

        // Corresponding distance matrix
        let p_str = r"5
a    0  17 21 31 23
b    17 0  30 34 21
c    21 30 0  28 39
d    31 34 28 0  43
e    23 21 39 43 0 
";
        let dm = DistanceMatrix::<f64>::from_phylip_strict(p_str, true).unwrap();
        let built = dm.upgma().unwrap();

        // Check topologies and branch lengths are the same
        let rf = expected.robinson_foulds(&built).unwrap();
        let wrf = expected.weighted_robinson_foulds(&built).unwrap();
        assert_eq!(wrf, 0.0);
        assert_eq!(rf, 0,);
    }
}
