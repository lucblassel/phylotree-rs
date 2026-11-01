//! Compute and manipulate phylogenetic distance matrices
//!

use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    str::FromStr,
};
#[cfg(not(target_arch = "wasm32"))]
use std::{fs, path::Path};

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
    /// The compact lower triangle has the wrong number of distances.
    #[error(
        "Distance matrix for {taxa} taxa requires {expected} distances, but received {actual}"
    )]
    DistanceCountMismatch {
        /// Number of taxa in the matrix.
        taxa: usize,
        /// Required number of pairwise distances.
        expected: usize,
        /// Supplied number of pairwise distances.
        actual: usize,
    },
    /// The number of taxa is too large to represent its pair count as a [`usize`].
    #[error("Distance matrix for {taxa} taxa is too large")]
    MatrixTooLarge {
        /// Number of taxa requested.
        taxa: usize,
    },
    /// Taxon labels must uniquely identify matrix rows and columns.
    #[error("Duplicate taxon label '{0}'")]
    DuplicateTaxon(String),
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
    /// A matrix index is outside the valid range.
    #[error("Matrix index {index} is out of bounds for a matrix with {size} taxa")]
    IndexOutOfBounds {
        /// Requested row or column index.
        index: usize,
        /// Number of taxa in the matrix.
        size: usize,
    },
    /// We are trying to get the compact pair index for the same leaf.
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
    #[error("Could not parse distance '{value}' at row {row}, column {column}")]
    DistParseError {
        /// One-based matrix row.
        row: usize,
        /// One-based distance column within the row.
        column: usize,
        /// Input token that could not be parsed.
        value: String,
    },
    /// A parsed distance is not finite.
    #[error("Distance at row {row}, column {column} is not finite")]
    NonFiniteDistance {
        /// One-based matrix row.
        row: usize,
        /// One-based distance column within the row.
        column: usize,
    },
    /// There are too few distances in one of the matrix rows.
    #[error("Missing distance from matrix row {0}")]
    MissingDistance(usize),
    /// There are too many distances in one of the matrix rows.
    #[error("Unexpected extra distance in matrix row {0}")]
    UnexpectedDistance(usize),
    /// The size of the matrix and the number of rows do not match
    #[error("Size and number of rows do not match: {0} rows for size {1}")]
    SizeAndRowsMismatch(usize, usize),
    /// Taxon labels must be unique.
    #[error("Duplicate taxon label '{0}'")]
    DuplicateTaxon(String),
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
    /// Identifiers of the taxa. Its length defines the matrix dimensions.
    taxa: Vec<String>,
    /// Distances in row-major order over the strict lower triangle.
    matrix: Vec<T>,
    /// Distance value for identical taxa.
    zero: T,
}

impl<T> DistanceMatrix<T>
where
    T: PairwiseDist,
{
    /// Creates a distance matrix from taxon labels and a compact lower triangle.
    ///
    /// The distances must be in row-major order over the strict lower triangle,
    /// containing exactly `n * (n - 1) / 2` values for `n` taxa.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let matrix = DistanceMatrix::new(vec!["A".into(), "B".into(), "C".into()], &[1.0, 2.0, 3.0]).unwrap();
    /// assert_eq!(*matrix.get("B", "C").unwrap(), 3.0);
    /// ```
    pub fn new(taxa: Vec<String>, matrix: &[T]) -> Result<Self, MatrixError> {
        Self::from_precomputed(taxa, matrix.to_vec())
    }

    /// Creates a zero-filled distance matrix for the supplied taxa.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let matrix = DistanceMatrix::<f64>::zeros(vec!["A".into(), "B".into()]).unwrap();
    /// assert_eq!(*matrix.get_by_index(0, 1).unwrap(), 0.0);
    /// ```
    pub fn zeros(taxa: Vec<String>) -> Result<Self, MatrixError> {
        let n_pairs = Self::pair_count(taxa.len())?;
        Self::from_precomputed(taxa, vec![zero(); n_pairs])
    }

    pub(crate) fn pair_count(size: usize) -> Result<usize, MatrixError> {
        let (left, right) = if size & 1 == 0 {
            (size / 2, size.saturating_sub(1))
        } else {
            (size.saturating_sub(1) / 2, size)
        };
        left.checked_mul(right)
            .ok_or(MatrixError::MatrixTooLarge { taxa: size })
    }

    fn validate_taxa(taxa: &[String]) -> Result<(), MatrixError> {
        let mut seen = HashMap::with_capacity(taxa.len());
        for (index, taxon) in taxa.iter().enumerate() {
            if seen.insert(taxon.as_str(), index).is_some() {
                return Err(MatrixError::DuplicateTaxon(taxon.clone()));
            }
        }
        Ok(())
    }

    /// Creates a checked distance matrix without cloning its compact storage.
    ///
    /// Values are stored in row-major order over the strict lower triangle.
    pub(crate) fn from_precomputed(taxa: Vec<String>, matrix: Vec<T>) -> Result<Self, MatrixError> {
        Self::validate_taxa(&taxa)?;
        let expected = Self::pair_count(taxa.len())?;
        if matrix.len() != expected {
            return Err(MatrixError::DistanceCountMismatch {
                taxa: taxa.len(),
                expected,
                actual: matrix.len(),
            });
        }

        Ok(Self {
            taxa,
            matrix,
            zero: zero(),
        })
    }

    /// Replaces the taxon labels without changing matrix dimensions or values.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// matrix.set_taxa(vec!["X".into(), "Y".into()]).unwrap();
    /// assert_eq!(*matrix.get("X", "Y").unwrap(), 2.0);
    /// ```
    pub fn set_taxa(&mut self, taxa: Vec<String>) -> Result<(), MatrixError> {
        if taxa.len() != self.len() {
            return Err(MatrixError::SizeError {
                size: self.len(),
                n_taxa: taxa.len(),
            });
        }
        Self::validate_taxa(&taxa)?;
        self.taxa = taxa;
        Ok(())
    }

    /// Returns the number of taxa in the matrix.
    pub fn len(&self) -> usize {
        self.taxa.len()
    }

    /// Returns `true` if the matrix contains no taxa.
    pub fn is_empty(&self) -> bool {
        self.taxa.is_empty()
    }

    /// Returns the taxon labels in matrix order.
    pub fn taxa(&self) -> &[String] {
        &self.taxa
    }

    /// Returns the compact strict lower-triangle storage.
    pub fn as_slice(&self) -> &[T] {
        &self.matrix
    }

    /// Get numerical index associated to taxon identifier
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(matrix.get_taxa_index("B").unwrap(), 1);
    /// ```
    pub fn get_taxa_index(&self, id: &str) -> Result<usize, MatrixError> {
        self.taxa
            .iter()
            .position(|v| *v == id)
            .ok_or_else(|| MatrixError::MissingTaxon(id.to_string()))
    }

    // Convert (row,col) matrix index to lower triangular linear index
    fn tril_to_vec_index(&self, i: usize, j: usize) -> Result<usize, MatrixError> {
        self.validate_index(i)?;
        self.validate_index(j)?;
        if i == j {
            Err(MatrixError::IndexError)
        } else {
            Ok(tril_to_rowvec_index(self.len(), i, j))
        }
    }

    fn validate_index(&self, index: usize) -> Result<(), MatrixError> {
        if index < self.len() {
            Ok(())
        } else {
            Err(MatrixError::IndexOutOfBounds {
                index,
                size: self.len(),
            })
        }
    }

    fn distance_by_index_unchecked(&self, i: usize, j: usize) -> &T {
        if i == j {
            &self.zero
        } else {
            &self.matrix[tril_to_rowvec_index(self.len(), i, j)]
        }
    }

    // Convert lower triangular linear index to (row,col) matrix index
    fn vec_to_tril_index(&self, k: usize) -> Result<(usize, usize), MatrixError> {
        if k >= self.matrix.len() {
            Err(MatrixError::IndexError)
        } else {
            Ok(rowvec_to_tril_index(self.len(), k))
        }
    }

    /// Gets the distance between two taxa by matrix index.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(*matrix.get_by_index(0, 1).unwrap(), 2.0);
    /// ```
    pub fn get_by_index(&self, i: usize, j: usize) -> Result<&T, MatrixError> {
        self.validate_index(i)?;
        self.validate_index(j)?;
        Ok(self.distance_by_index_unchecked(i, j))
    }

    /// Gets the distance between two taxa by label.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(*matrix.get("A", "B").unwrap(), 2.0);
    /// ```
    pub fn get(&self, id_1: &str, id_2: &str) -> Result<&T, MatrixError> {
        let i = self.get_taxa_index(id_1)?;
        let j = self.get_taxa_index(id_2)?;
        self.get_by_index(i, j)
    }

    /// Sets the distance between two taxa by matrix index.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// matrix.set_by_index(0, 1, 3.0).unwrap();
    /// assert_eq!(*matrix.get_by_index(1, 0).unwrap(), 3.0);
    /// ```
    pub fn set_by_index(&mut self, i: usize, j: usize, dist: T) -> Result<(), MatrixError> {
        self.validate_index(i)?;
        self.validate_index(j)?;
        if i == j {
            return if dist == self.zero {
                Ok(())
            } else {
                Err(MatrixError::NonZeroIdenticalDistance)
            };
        }

        let idx = self.tril_to_vec_index(i, j)?;
        self.matrix[idx] = dist;
        Ok(())
    }

    /// Sets the distance between two taxa by label.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// matrix.set("A", "B", 3.0).unwrap();
    /// assert_eq!(*matrix.get("B", "A").unwrap(), 3.0);
    /// ```
    pub fn set(&mut self, id_1: &str, id_2: &str, dist: T) -> Result<(), MatrixError> {
        let i = self.get_taxa_index(id_1)?;
        let j = self.get_taxa_index(id_2)?;
        self.set_by_index(i, j, dist)
    }

    /// Get the distance matrix as a HashMap containing taxa pairs as keys
    /// and pairwise distances as values. This contains 2 copies of each distance
    /// (one for each permutation of the pair of names), and diagonal zero entries.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(matrix.to_map()[&("A".into(), "B".into())], 2.0);
    /// ```
    pub fn to_map(&self) -> HashMap<(String, String), T> {
        HashMap::from_iter(self.taxa.iter().enumerate().flat_map(|(i, taxon1)| {
            self.taxa.iter().enumerate().map(move |(j, taxon2)| {
                (
                    (taxon1.clone(), taxon2.clone()),
                    *self.distance_by_index_unchecked(i, j),
                )
            })
        }))
    }

    /// Outputs a matrix as a phylip formatted string
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// let text = matrix.to_phylip(true).unwrap();
    /// assert_eq!(DistanceMatrix::<f64>::from_phylip_strict(&text, true).unwrap().len(), 2);
    /// ```
    pub fn to_phylip(&self, square: bool) -> Result<String, MatrixError> {
        if self.is_empty() {
            return Ok("0\n".to_owned());
        }

        let body = self
            .taxa
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let lim = if square { self.len() } else { i };
                let row_s = (0..lim)
                    .map(|j| {
                        let d = if i == j {
                            zero()
                        } else {
                            *self.distance_by_index_unchecked(i, j)
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

        Ok(format!("{}\n{body}\n", self.len()))
    }

    /// Writes the matrix to a phylip file
    ///
    /// ```no_run
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// matrix.to_file(std::path::Path::new("distances.phy"), true).unwrap();
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
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
            .enumerate()
            .map(|(column, d)| {
                d.parse().map_err(|_| PhylipParseError::DistParseError {
                    row: row_num + 1,
                    column: column + 1,
                    value: d.to_owned(),
                })
            })
            .take(if tril { row_num } else { size })
            .collect::<Result<Vec<_>, _>>()?;

        Ok((name.to_string(), dists))
    }

    /// Build a distance matrix from the lower triangle of a phylip file.
    /// Does not enforce symmetry in the phylip file, it should be faster than
    /// [`DistanceMatrix::from_phylip_strict`]. Also does not check if all the rows are the same size
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let matrix = DistanceMatrix::<f64>::from_phylip_tril("2\nA\nB 2\n").unwrap();
    /// assert_eq!(*matrix.get("A", "B").unwrap(), 2.0);
    /// ```
    pub fn from_phylip_tril(phylip: &str) -> Result<Self, PhylipParseError<T>> {
        let mut lines = phylip.lines();
        let size: usize = lines
            .next()
            .ok_or(PhylipParseError::EmptyMatrixFile)?
            .parse()
            .map_err(PhylipParseError::SizeParseError)?;

        let mut taxa = vec![];
        let mut max_v = Vec::with_capacity(Self::pair_count(size)?);

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

    /// Builds a distance matrix from a PHYLIP-formatted string.
    ///
    /// The parser requires exactly the declared number of rows and exactly the
    /// expected number of distances in every row. Square matrices must have a
    /// zero diagonal and exact symmetry. All parsed distances must be finite.
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let matrix = DistanceMatrix::<f64>::from_phylip_strict("2\nA 0 2\nB 2 0\n", true).unwrap();
    /// assert_eq!(*matrix.get("A", "B").unwrap(), 2.0);
    /// ```
    pub fn from_phylip_strict(phylip: &str, square: bool) -> Result<Self, PhylipParseError<T>> {
        let mut lines = phylip.lines();
        let size = lines
            .next()
            .ok_or(PhylipParseError::EmptyMatrixFile)?
            .trim()
            .parse()
            .map_err(PhylipParseError::SizeParseError)?;
        let input_rows: Vec<_> = lines.collect();

        if input_rows.len() != size {
            return Err(PhylipParseError::SizeAndRowsMismatch(
                input_rows.len(),
                size,
            ));
        }

        let mut names = Vec::with_capacity(size);
        let mut rows = Vec::with_capacity(size);
        let mut seen_names = HashMap::with_capacity(size);

        for (row_index, line) in input_rows.into_iter().enumerate() {
            let mut fields = line.split_whitespace();
            let name = fields
                .next()
                .ok_or(PhylipParseError::EmptyRow(row_index + 1))?;
            if seen_names.insert(name, row_index).is_some() {
                return Err(PhylipParseError::DuplicateTaxon(name.to_owned()));
            }

            let mut distances = Vec::new();
            for (column_index, value) in fields.enumerate() {
                let distance =
                    value
                        .parse::<T>()
                        .map_err(|_| PhylipParseError::DistParseError {
                            row: row_index + 1,
                            column: column_index + 1,
                            value: value.to_owned(),
                        })?;
                if !distance.is_finite() {
                    return Err(PhylipParseError::NonFiniteDistance {
                        row: row_index + 1,
                        column: column_index + 1,
                    });
                }
                distances.push(distance);
            }

            let expected = if square { size } else { row_index };
            if distances.len() < expected {
                return Err(PhylipParseError::MissingDistance(row_index + 1));
            }
            if distances.len() > expected {
                return Err(PhylipParseError::UnexpectedDistance(row_index + 1));
            }

            names.push(name.to_owned());
            rows.push(distances);
        }

        if square {
            for i in 0..size {
                if rows[i][i] != zero() {
                    return Err(PhylipParseError::NonZeroDiagonalValue(names[i].clone()));
                }
                for (j, previous_row) in rows.iter().enumerate().take(i) {
                    if rows[i][j] != previous_row[i] {
                        return Err(PhylipParseError::NonSymmetric(rows[i][j], previous_row[i]));
                    }
                }
            }
        }

        let mut compact = Vec::with_capacity(Self::pair_count(size)?);
        for (i, row) in rows.into_iter().enumerate() {
            compact.extend(row.into_iter().take(i));
        }
        Ok(Self::from_precomputed(names, compact)?)
    }

    /// Find smallest distance in the distance matrix as well as the
    /// corresponding index
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(matrix.min().unwrap().1, 2.0);
    /// ```
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
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(matrix.max().unwrap().1, 2.0);
    /// ```
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
    ///
    /// ```no_run
    /// use phylotree::distance::DistanceMatrix;
    /// let matrix = DistanceMatrix::<f64>::from_file(std::path::Path::new("distances.phy"), true).unwrap();
    /// println!("{} taxa", matrix.len());
    /// ```
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &Path, square: bool) -> Result<Self, PhylipParseError<T>> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_phylip_strict(&newick_string, square)
    }

    /// Iterator over the lower triangle of the matrix
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(matrix.iter().copied().collect::<Vec<_>>(), vec![2.0]);
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = &'_ T> {
        self.matrix.iter()
    }

    /// Iterator over the lower triangle of the matrix with the corresponding
    /// (row,col) matrix indices
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// assert_eq!(matrix.indexed_iter().next().unwrap(), ((1, 0), &2.0));
    /// ```
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
    ///
    /// ```
    /// use phylotree::distance::DistanceMatrix;
    /// let mut matrix = DistanceMatrix::new(vec!["A".into(), "B".into()], &[2.0]).unwrap();
    /// let tree = matrix.upgma().unwrap();
    /// assert_eq!(tree.n_leaves(), 2);
    /// ```
    pub fn upgma(&self) -> Result<Tree, MatrixError> {
        // Setup
        let mut dm = (*self).clone();
        let size = self.len();
        let mut card = vec![1; size];
        let mut merged = vec![false; size];
        let mut heights = vec![0.; size];
        let mut n_clusters = size;

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
                .merge_children(node_ids[a], node_ids[b], Some(d_au), Some(d_bu), None, None)
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

        let root_node = tree.get_mut(virt_root).unwrap();
        root_node.set_child_edge(a, Some(d_ar));
        root_node.set_child_edge(b, Some(d_br));
        tree.get_mut(a).unwrap().parent_edge = Some(d_ar);
        tree.get_mut(b).unwrap().parent_edge = Some(d_br);

        Ok(tree)
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
        let mut matrix = DistanceMatrix::zeros(
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
    fn checked_construction_preserves_matrix_invariants() {
        let empty = DistanceMatrix::<f32>::new(vec![], &[]).unwrap();
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
        assert!(empty.taxa().is_empty());
        assert!(empty.as_slice().is_empty());
        assert!(empty.to_map().is_empty());
        assert_eq!(empty.to_phylip(true).unwrap(), "0\n");

        let singleton = DistanceMatrix::<f32>::new(vec!["A".into()], &[]).unwrap();
        assert_eq!(singleton.len(), 1);
        assert_eq!(*singleton.get("A", "A").unwrap(), 0.0);
        assert_eq!(
            singleton.to_map().get(&("A".into(), "A".into())),
            Some(&0.0)
        );

        assert!(matches!(
            DistanceMatrix::<f32>::new(vec!["A".into(), "B".into()], &[]),
            Err(MatrixError::DistanceCountMismatch {
                taxa: 2,
                expected: 1,
                actual: 0,
            })
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::new(vec!["A".into()], &[1.0]),
            Err(MatrixError::DistanceCountMismatch {
                taxa: 1,
                expected: 0,
                actual: 1,
            })
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::new(vec!["A".into(), "A".into()], &[1.0]),
            Err(MatrixError::DuplicateTaxon(taxon)) if taxon == "A"
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::pair_count(usize::MAX),
            Err(MatrixError::MatrixTooLarge { .. })
        ));
    }

    #[test]
    fn indexed_and_named_access_validate_inputs() {
        let mut matrix =
            DistanceMatrix::<f32>::new(vec!["A".into(), "B".into(), "C".into()], &[1.0, 2.0, 3.0])
                .unwrap();

        assert_eq!(*matrix.get_by_index(0, 0).unwrap(), 0.0);
        assert_eq!(*matrix.get_by_index(0, 2).unwrap(), 2.0);
        assert_eq!(*matrix.get_by_index(2, 0).unwrap(), 2.0);
        assert!(matches!(
            matrix.get_by_index(3, 0),
            Err(MatrixError::IndexOutOfBounds { index: 3, size: 3 })
        ));
        assert!(matches!(
            matrix.get("missing", "missing"),
            Err(MatrixError::MissingTaxon(taxon)) if taxon == "missing"
        ));

        matrix.set_by_index(0, 2, 4.0).unwrap();
        assert_eq!(*matrix.get("A", "C").unwrap(), 4.0);
        assert_eq!(*matrix.get("C", "A").unwrap(), 4.0);
        assert!(matches!(
            matrix.set_by_index(1, 1, 1.0),
            Err(MatrixError::NonZeroIdenticalDistance)
        ));
        assert!(matrix.set_by_index(1, 1, 0.0).is_ok());

        let original_taxa = matrix.taxa().to_vec();
        assert!(matches!(
            matrix.set_taxa(vec!["A".into(), "A".into(), "C".into()]),
            Err(MatrixError::DuplicateTaxon(taxon)) if taxon == "A"
        ));
        assert_eq!(matrix.taxa(), original_taxa);
    }

    #[test]
    fn to_map_contains_both_directions_and_diagonal_zeroes() {
        let matrix = DistanceMatrix::<f32>::new(vec!["A".into(), "B".into()], &[2.5]).unwrap();
        let map = matrix.to_map();

        assert_eq!(map.len(), 4);
        assert_eq!(map.get(&("A".into(), "A".into())), Some(&0.0));
        assert_eq!(map.get(&("B".into(), "B".into())), Some(&0.0));
        assert_eq!(map.get(&("A".into(), "B".into())), Some(&2.5));
        assert_eq!(map.get(&("B".into(), "A".into())), Some(&2.5));
    }

    #[test]
    fn strict_phylip_accepts_empty_and_singleton_matrices() {
        let empty = DistanceMatrix::<f32>::from_phylip_strict("0\n", true).unwrap();
        assert!(empty.is_empty());
        assert_eq!(empty.to_phylip(true).unwrap(), "0\n");

        let square = DistanceMatrix::<f32>::from_phylip_strict("1\nA 0\n", true).unwrap();
        assert_eq!(square.taxa(), &["A"]);
        assert_eq!(*square.get("A", "A").unwrap(), 0.0);

        let triangle = DistanceMatrix::<f32>::from_phylip_strict("1\nA\n", false).unwrap();
        assert_eq!(triangle.taxa(), &["A"]);
        assert_eq!(triangle.to_phylip(false).unwrap(), "1\nA\n");
    }

    #[test]
    fn strict_phylip_rejects_wrong_row_and_column_counts() {
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("1\nA 0\nB 0\n", true),
            Err(PhylipParseError::SizeAndRowsMismatch(2, 1))
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0 1\n", true),
            Err(PhylipParseError::SizeAndRowsMismatch(1, 2))
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0\nB 1 0\n", true),
            Err(PhylipParseError::MissingDistance(1))
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0 1 2\nB 1 0\n", true),
            Err(PhylipParseError::UnexpectedDistance(1))
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0 1 bad\nB 1 0\n", true),
            Err(PhylipParseError::DistParseError {
                row: 1,
                column: 3,
                ..
            })
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA\n\n", false),
            Err(PhylipParseError::EmptyRow(2))
        ));
    }

    #[test]
    fn strict_phylip_rejects_duplicate_nonfinite_and_invalid_square_values() {
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0 1\nA 1 0\n", true),
            Err(PhylipParseError::DuplicateTaxon(taxon)) if taxon == "A"
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0 NaN\nB NaN 0\n", true),
            Err(PhylipParseError::NonFiniteDistance { row: 1, column: 2 })
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 1 2\nB 2 0\n", true),
            Err(PhylipParseError::NonZeroDiagonalValue(taxon)) if taxon == "A"
        ));
        assert!(matches!(
            DistanceMatrix::<f32>::from_phylip_strict("2\nA 0 1\nB 2 0\n", true),
            Err(PhylipParseError::NonSymmetric(2.0, 1.0))
        ));
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
        let rf = expected.unrooted_robinson_foulds(&built).unwrap();
        let wrf = expected.weighted_robinson_foulds(&built).unwrap();
        assert_eq!(wrf, 0.0);
        assert_eq!(rf, 0,);
    }
}
