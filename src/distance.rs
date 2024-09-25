//! Compute and manipulate phylogenetic distance matrices
//!

use std::{
    fmt::{Debug, Display},
    fs,
    path::Path,
    str::FromStr,
    vec,
};

use itertools::{iproduct, Itertools};
use ndarray::{s, Array1, Array2, Axis};
use num_traits::{zero, Float, Zero};
use thiserror::Error;

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
pub enum ParseError<T>
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

#[derive(Debug)]
/// A phylogenetic distance matrix
pub struct DistanceMatrix<T> {
    /// Number of taxa in the matrix
    pub size: usize,
    /// Identifiers of the taxa
    pub taxa: Vec<String>,
    /// Distances between taxa
    matrix: Array2<T>,
}

impl<T> DistanceMatrix<T>
where
    T: Display + Debug + Float + Zero + FromStr,
{
    /// Build a distance matrix from a phylip formatted string
    pub fn from_phylip(phylip: &str, square: bool) -> Result<Self, ParseError<T>> {
        let mut lines = phylip.lines();
        let size: usize = lines
            .next()
            .ok_or(ParseError::EmptyMatrixFile)?
            .parse()
            .map_err(ParseError::SizeParseError)?;

        let mut matrix = Array2::<T>::zeros((size, size));

        let mut taxa = Vec::with_capacity(size);

        for (i, line) in lines.enumerate() {
            let mut fields = line.split_whitespace();
            let name = fields.next().ok_or(ParseError::EmptyRow(i))?;
            let dists = fields
                .map(|d| d.parse().map_err(|_| ParseError::DistParseError))
                .collect::<Result<Vec<_>, _>>()?;
            let dists = Array1::from_vec(dists);

            if square && dists.len() != size || !square && dists.len() != i {
                return Err(ParseError::MissingDistance(i + 1));
            }

            if square {
                matrix.slice_mut(s![i, ..]).assign(&dists);
            } else {
                matrix.slice_mut(s![i, ..i]).assign(&dists);
                matrix.slice_mut(s![..i, i]).assign(&dists);
            }

            if (*matrix.get((i, i)).unwrap()) != zero() {
                return Err(ParseError::NonZeroDiagonalValue(name.to_string()));
            }

            taxa.push(name.to_string());
        }

        if taxa.len() != size {
            return Err(ParseError::SizeAndRowsMismatch(taxa.len(), size));
        }

        // Check that matrix is symmetric
        if !(matrix.t() == matrix) {
            return Err(ParseError::NonSymmetricMat);
        }

        Ok(DistanceMatrix { size, taxa, matrix })
    }

    /// Build phylip formatted representation of the distance matrix
    pub fn to_phylip(&self, square: bool) -> Result<String, MatrixError> {
        let body = self
            .matrix
            .axis_iter(Axis(0))
            .enumerate()
            .map(|(i, row)| {
                let lim = if square { self.size } else { i };
                let row_s = row.iter().take(lim).map(|e| e.to_string()).join("  ");

                let mut out = self.taxa[i].to_string();
                if !row_s.is_empty() {
                    out.push_str(&format!("    {row_s}"));
                }
                out

                //format!("{}    {row_s}", self.taxa[i])
            })
            .join("\n");

        Ok(format!("{}\n{body}\n", self.size))
    }

    /// Write distance matrix to file in Phylip format
    pub fn to_file(&self, path: &Path, square: bool) -> Result<(), MatrixError> {
        match fs::write(path, self.to_phylip(square)?) {
            Ok(_) => Ok(()),
            Err(e) => Err(MatrixError::IoError(e)),
        }
    }

    /// Find minimum non-zero distance and associated indices
    pub fn min(self) -> Option<(T, (usize, usize))> {
        self.matrix.indexed_iter().fold(None, |a, ((i, j), v)| {
            if i == j {
                return a; // Skip diag
            }

            match a {
                None => Some((*v, (i, j))),
                Some((mut a_v, (mut a_i, mut a_j))) => {
                    if *v < a_v {
                        a_i = i;
                        a_j = j;
                        a_v = *v;
                    }

                    Some((a_v, (a_i, a_j)))
                }
            }
        })
    }

    /// Get numerical index associated to taxon identifier
    pub fn get_taxa_index(&self, id: &str) -> Result<usize, MatrixError> {
        self.taxa
            .iter()
            .find_position(|v| *v == id)
            .ok_or(MatrixError::MissingTaxon(id.to_string()))
            .map(|(i, _)| i)
    }

    /// get distance between 2 taxa
    pub fn get(&self, id_1: &str, id_2: &str) -> Result<&T, MatrixError> {
        let i1 = self.get_taxa_index(id_1)?;
        let i2 = self.get_taxa_index(id_2)?;

        self.matrix.get((i1, i2)).ok_or(MatrixError::IndexError)
    }

    /// Set the taxa identifers of the distance matrix
    pub fn set_taxa(&mut self, taxa: Vec<String>) -> Result<(), MatrixError> {
        if self.size != taxa.len() {
            return Err(MatrixError::SizeError {
                size: self.size,
                n_taxa: taxa.len(),
            });
        }

        self.taxa = taxa;

        Ok(())
    }

    /// Set the distance between 2 taxa
    pub fn set(&mut self, id_1: &str, id_2: &str, dist: T) -> Result<(), MatrixError> {
        let i1 = self.get_taxa_index(id_1)?;
        let i2 = self.get_taxa_index(id_2)?;

        // Set value in both symmetric entries
        for index in vec![(i1, i2), (i2, i1)] {
            let v_ptr = self.matrix.get_mut(index).ok_or(MatrixError::IndexError)?;
            *v_ptr = dist;
        }

        Ok(())
    }

    /// Allocate space for a distance matrix of given number of taxa
    pub fn new_with_size(size: usize) -> Self {
        Self {
            size,
            taxa: Vec::with_capacity(size),
            matrix: Array2::<T>::zeros((size, size)),
        }
    }

    pub(crate) fn from_precomputed(taxa: Vec<String>, matrix: Vec<T>) -> Result<Self, MatrixError> {
        // Check that sizes are OK
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

        let mut m = Self::new_with_size(taxa.len());
        m.set_taxa(taxa)?;

        for ((i, j), v) in iproduct!(0..n, 0..n)
            .filter(|(i, j)| i < j)
            .zip(matrix.iter())
        {
            *(m.matrix.get_mut((i, j)).ok_or(MatrixError::IndexError)?) = *v;
            *(m.matrix.get_mut((j, i)).ok_or(MatrixError::IndexError)?) = *v;
        }

        Ok(m)
    }
}

#[cfg(test)]
mod tests {

    use core::panic;

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
    fn from_phylip() -> Result<(), ParseError<f32>> {
        let build: DistanceMatrix<f32> = DistanceMatrix::from_phylip(SQUARE, true)?;
        assert_eq!(
            SQUARE,
            build.to_phylip(true).unwrap(),
            "{SQUARE}\n{}",
            build.to_phylip(true).unwrap()
        );

        let build = DistanceMatrix::from_phylip(TRIANGLE, false)?;
        assert_eq!(
            TRIANGLE,
            build.to_phylip(false).unwrap(),
            "{TRIANGLE}\n{}",
            build.to_phylip(false).unwrap()
        );

        Ok(())
    }

    #[test]
    fn from_phylip_errors() {
        let square_nonsym = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        let matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip(square_nonsym, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            ParseError::NonSymmetric(_, _) => {}
            ParseError::NonSymmetricMat => {}
            _ => panic!("Error should be 'ParseError::NonSymmetric' not: {err}"),
        }

        let square_missing_dist = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0
s5    5  10  15  0
";
        let mut matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip(square_missing_dist, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            ParseError::MissingDistance(_) => {}
            _ => panic!("Error should be 'ParseError::MissingDistance' not: {err}"),
        }

        let square_missing_row = "4
s1    0  2  3  7
s2    2  0  6  10
s5    4  6  0  5
";
        matrix = DistanceMatrix::from_phylip(square_missing_row, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            ParseError::SizeAndRowsMismatch(_, _) => {}
            _ => panic!("Error should be 'ParseError::SizeAndRowsMismatch' not: {err}"),
        }

        let square_missing_size = "s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        matrix = DistanceMatrix::from_phylip(square_missing_size, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            ParseError::SizeParseError(_) => {}
            _ => panic!("Error should be 'ParseError::SizeParseError' not: {err}"),
        }

        let square_nonzero = "4
s1    1  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    7  10  15  0
";
        let matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip(square_nonzero, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            ParseError::NonZeroDiagonalValue(_) => {}
            _ => panic!("Error should be 'ParseError::NonZeroDiagonalValue' not: {err}"),
        }
    }

    #[test]
    fn test_dm2() {
        let mat = DistanceMatrix::<f32>::from_phylip(TRIANGLE, false).unwrap();
        let (v, idx) = mat.min().unwrap();
        assert_eq!(v, 2.0);
        assert!(idx == (0, 1) || idx == (1, 0));
    }

    #[test]
    fn from_phylip2() -> Result<(), ParseError<f64>> {
        let build: DistanceMatrix<f64> = DistanceMatrix::from_phylip(SQUARE, true)?;

        assert_eq!(
            SQUARE,
            build.to_phylip(true).unwrap(),
            "{SQUARE}\n{}",
            build.to_phylip(true).unwrap()
        );

        let build = DistanceMatrix::from_phylip(TRIANGLE, false)?;
        assert_eq!(
            TRIANGLE,
            build.to_phylip(false).unwrap(),
            "{TRIANGLE}\n{}",
            build.to_phylip(false).unwrap()
        );

        Ok(())
    }
}
