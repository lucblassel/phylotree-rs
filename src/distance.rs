//! Compute and manipulate phylogenetic distance matrices
//!

use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    fs,
    path::Path,
    str::FromStr,
    vec,
};

use itertools::Itertools;
use num_traits::{Float, Zero};
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
    /// The square phylip matrix is not symmetric
    #[error("Non symetric matrix: {0} and {1} are different")]
    NonSymmetric(T, T),
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
    matrix: Vec<T>,
    /// Distance value for identical taxa
    zero: T,
}

impl<T> DistanceMatrix<T>
where
    T: Display + Debug + Float + Zero + FromStr,
{
    /// Create a new distance matrix for a certain number of sequences
    pub fn new(taxa: Vec<String>, matrix: Vec<T>) -> Self {
        Self {
            size: taxa.len(),
            taxa,
            matrix,
            zero: Zero::zero(),
        }
    }

    /// Create an empty distance matrix with a given size
    pub fn new_with_size(size: usize) -> Self {
        Self {
            size,
            taxa: Vec::with_capacity(size),
            matrix: vec![Zero::zero(); size * (size - 1) / 2],
            zero: Zero::zero(),
        }
    }

    /// Create a distance matrix from pre-computed values. The `matrix` parameter
    /// represents the upper triangle of the distance matrix as a single vector.
    /// The matrix index to vector index formula can be found in the [`DistanceMatrix.get_index`]
    /// function.
    pub(crate) fn from_precomputed(taxa: Vec<String>, matrix: Vec<T>) -> Self {
        Self {
            size: taxa.len(),
            taxa,
            matrix,
            zero: Zero::zero(),
        }
    }

    /// Set the taxa of the matrix
    pub fn set_taxa(&mut self, taxa: Vec<String>) -> Result<(), MatrixError> {
        if taxa.len() != self.size {
            Err(MatrixError::SizeError {
                size: self.size,
                n_taxa: taxa.len(),
            })
        } else {
            self.taxa = taxa;
            Ok(())
        }
    }

    /// Get the index in the distance vector for 2 sequences
    fn get_index(&self, taxon1: &str, taxon2: &str) -> Result<usize, MatrixError> {
        if taxon1 == taxon2 {
            return Err(MatrixError::IndexError);
        }

        let mut i = self
            .taxa
            .iter()
            .position(|v| v == taxon1)
            .ok_or(MatrixError::MissingTaxon(taxon1.to_string()))?;

        let mut j = self
            .taxa
            .iter()
            .position(|v| v == taxon2)
            .ok_or(MatrixError::MissingTaxon(taxon2.to_string()))?;

        if j < i {
            std::mem::swap(&mut i, &mut j);
        }

        Ok((2 * self.size - 3 - i) * i / 2 + j - 1)
    }

    /// Get the distance between two sequences
    pub fn get(&self, id_1: &str, id_2: &str) -> Result<&T, MatrixError> {
        if id_1 == id_2 {
            Ok(&self.zero)
        } else {
            let idx = self.get_index(id_1, id_2)?;
            Ok(&self.matrix[idx])
        }
    }

    /// Set the distance between two sequences
    pub fn set(&mut self, id_1: &str, id_2: &str, dist: T) -> Result<(), MatrixError> {
        if id_1 == id_2 {
            if dist != self.zero {
                Err(MatrixError::NonZeroIdenticalDistance)
            } else {
                Ok(())
            }
        } else {
            let idx = self.get_index(id_1, id_2)?;
            self.matrix[idx] = dist;
            Ok(())
        }
    }

    /// Get the distance matrix as a HashMap containing taxa pairs as keys
    /// and pairwise distances as values
    pub fn to_map(&self) -> HashMap<(String, String), T> {
        HashMap::from_iter(self.taxa.iter().cartesian_product(self.taxa.iter()).map(
            |(taxon1, taxon2)| {
                let idx = self.get_index(taxon1, taxon2).unwrap();
                ((taxon1.clone(), taxon2.clone()), self.matrix[idx])
            },
        ))
    }

    /// Returns a string representing the distance matrix in square format
    fn to_phylip_square(&self) -> Result<String, MatrixError> {
        let names: Vec<_> = self.taxa.iter().cloned().sorted().collect();
        let mut output = format!("{}\n", self.size);

        for name1 in names.iter() {
            output += &format!("{name1}  ");
            for name2 in names.iter() {
                let d = if name1 != name2 {
                    *self
                        .get(name1, name2)
                        .map_err(|_| MatrixError::MissingDistance(name1.clone(), name2.clone()))?
                } else {
                    Zero::zero()
                };

                output += &format!("  {}", d);
            }
            output += "\n"
        }

        Ok(output)
    }

    /// Returns a string representing the distance matrix in triangle format
    fn to_phylip_triangle(&self) -> Result<String, MatrixError> {
        let names: Vec<_> = self.taxa.iter().cloned().sorted().collect();
        let mut output = format!("{}\n", self.size);

        for name1 in names.iter() {
            output += &format!("{name1}  ");
            for name2 in names.iter() {
                if name1 == name2 {
                    break;
                }
                let d = self
                    .get(name1, name2)
                    .map_err(|_| MatrixError::MissingDistance(name1.clone(), name2.clone()))?;
                output += &format!("  {}", d);
            }
            output += "\n"
        }

        Ok(output)
    }

    /// Outputs the matrix as a phylip formatted string
    pub fn to_phylip(&self, square: bool) -> Result<String, MatrixError> {
        if square {
            self.to_phylip_square()
        } else {
            self.to_phylip_triangle()
        }
    }

    /// Writes the matrix to a phylip file
    pub fn to_file(&self, path: &Path, square: bool) -> Result<(), MatrixError> {
        match fs::write(path, self.to_phylip(square)?) {
            Ok(_) => Ok(()),
            Err(e) => Err(MatrixError::IoError(e)),
        }
    }

    /// Build a distance matrix from a phylip formatted string
    pub fn from_phylip(phylip: &str, square: bool) -> Result<Self, ParseError<T>> {
        let mut lines = phylip.lines();
        let size = lines
            .next()
            .ok_or(ParseError::EmptyMatrixFile)?
            .parse()
            .map_err(ParseError::SizeParseError)?;

        let mut names = vec![];
        let mut rows = vec![];

        for (i, line) in lines.enumerate() {
            let mut fields = line.split_whitespace();
            let name = fields.next().ok_or(ParseError::EmptyRow(i))?;
            let dists: Result<Vec<_>, _> = fields
                .map(|d| d.parse::<T>().map_err(|_| ParseError::DistParseError))
                .collect();

            let dists = dists?;

            if square && dists.len() != size || !square && dists.len() != i {
                return Err(ParseError::MissingDistance(i + 1));
            }

            names.push(name);
            rows.push(dists);
        }

        if names.len() != size {
            return Err(ParseError::SizeAndRowsMismatch(names.len(), size));
        }

        let mut matrix = Self::new_with_size(size);
        matrix.set_taxa(names.iter().cloned().map(|v| v.to_string()).collect_vec())?;

        for (&n1, row) in names.iter().zip(rows) {
            for (&n2, dist) in names.iter().zip(row) {
                matrix.set(n1, n2, dist)?;
            }
        }

        Ok(matrix)
    }

    /// Reads the matrix from a phylip file
    pub fn from_file(path: &Path, square: bool) -> Result<Self, ParseError<T>> {
        let newick_string = fs::read_to_string(path)?;
        Self::from_phylip(&newick_string, square)
    }
}

#[cfg(test)]
mod tests {

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
        let names = vec![("s1", 1.0), ("s2", 2.0), ("s3", 3.0), ("s5", 5.0)];
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

        assert_eq!(SQUARE, matrix.to_phylip(true).unwrap());
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
        //         let square_nonsym = "4
        // s1    0  2  3  7
        // s2    2  0  6  10
        // s3    3  6  0  15
        // s5    5  10  15  0
        // ";
        //         let mut matrix: Result<DistanceMatrix<f32>, _> =
        //             DistanceMatrix::from_phylip(square_nonsym, true);
        //         assert!(matrix.is_err());
        //
        //         let err = matrix.err().unwrap();
        //         match err {
        //             ParseError::NonSymmetric(_, _) => {}
        //             _ => panic!("Error should be 'ParseError::NonSymmetric' not: {err}"),
        //         }

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
s5    5  10  15  0
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
    }
}
