use std::{
    collections::{hash_map::Iter, HashMap, HashSet},
    fmt::{Debug, Display},
    fs,
    path::Path,
    str::FromStr,
    vec,
};

use itertools::Itertools;
use num_traits::{Float, Zero};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Distance already exists, cannot overwrite it.")]
    OverwritingNotPermitted,
    #[error("You added more sequences than there is room for in the matrix.")]
    SizeExceeded,
    #[error("Missing distance between {0} and {1}.")]
    MissingDistance(String, String),
    #[error("Error writing file")]
    IOError(#[from] std::io::Error),
}

#[derive(Error, Debug)]
pub enum ParseError<T>
where
    T: Debug,
{
    #[error("Matrix file is empty.")]
    EmptyMatrixFile,
    #[error("Could not parse size from file.")]
    SizeParseError(#[from] std::num::ParseIntError),
    #[error("Row {0} is empty.")]
    EmptyRow(usize),
    #[error("Could not parse distance from file.")]
    DistParseError,
    #[error("Missing distance from matrix line {0}")]
    MissingDistance(usize),
    #[error("Size and number of rows do not match: {0} rows for size {1}")]
    SizeAndRowsMismatch(usize, usize),
    #[error("Non symetric matrix: {0} and {1} are different")]
    NonSymmetric(T, T),
    #[error("Error creating matrix.")]
    MatrixError(#[from] crate::distance::Error),
    #[error("Error reading file")]
    IOError(#[from] std::io::Error),
}

#[derive(Debug)]
/// Struct representing a phylogenetic distance matrix
pub struct DistanceMatrix<T> {
    pub size: usize,
    pub ids: HashSet<String>,
    matrix: HashMap<(String, String), T>,
}

impl<T> DistanceMatrix<T>
where
    T: Display + Debug + Float + Zero + FromStr,
{
    /// Create a new distance matrix for a certain number of sequences
    pub fn new(size: usize) -> Self {
        Self {
            size,
            ids: HashSet::with_capacity(size),
            matrix: HashMap::with_capacity(size * size - 1),
        }
    }

    /// Gets the correct numbering of sequence ids for the matrix key
    fn get_key(&self, id_1: &str, id_2: &str) -> (String, String) {
        if id_1 < id_2 {
            (id_1.to_owned(), id_2.to_owned())
        } else {
            (id_2.to_owned(), id_1.to_owned())
        }
    }

    /// Get the distance between two sequences
    pub fn get(&self, id_1: &str, id_2: &str) -> Option<&T> {
        self.matrix.get(&self.get_key(id_1, id_2))
    }

    /// Get the mutable distance between two sequences
    pub fn get_mut(&mut self, id_1: &str, id_2: &str) -> Option<&mut T> {
        self.matrix.get_mut(&self.get_key(id_1, id_2))
    }

    /// Return iterator over the distance matrix
    pub fn iter(&self) -> Iter<(String, String), T> {
        self.matrix.iter()
    }

    /// Set an entry in the distance matrix, if overwriting is permitted and the key
    /// exists it returns the old value as `Some(old)`
    pub fn set(
        &mut self,
        id_1: &str,
        id_2: &str,
        dist: T,
        overwrite: bool,
    ) -> Result<Option<T>, Error> {
        if let Some(old_dist) = self.get_mut(id_1, id_2) {
            if overwrite {
                let copy = *old_dist;
                *old_dist = dist;
                return Ok(Some(copy));
            } else {
                return Err(Error::OverwritingNotPermitted);
            }
        }

        self.matrix.insert(self.get_key(id_1, id_2), dist);
        self.ids.insert(id_1.to_owned());
        self.ids.insert(id_2.to_owned());

        if self.ids.len() > self.size {
            Err(Error::SizeExceeded)
        } else {
            Ok(None)
        }
    }

    /// Returns a string representing the distance matrix in square format
    fn to_phylip_square(&self) -> Result<String, Error> {
        let names: Vec<_> = self.ids.iter().cloned().sorted().collect();
        let mut output = format!("{}\n", self.size);

        for name1 in names.iter() {
            output += &format!("{name1}  ");
            for name2 in names.iter() {
                let d = if name1 != name2 {
                    *self
                        .get(name1, name2)
                        .ok_or::<Error>(Error::MissingDistance(name1.clone(), name2.clone()))?
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
    fn to_phylip_triangle(&self) -> Result<String, Error> {
        let names: Vec<_> = self.ids.iter().cloned().sorted().collect();
        let mut output = format!("{}\n", self.size);

        for name1 in names.iter() {
            output += &format!("{name1}  ");
            for name2 in names.iter() {
                if name1 == name2 {
                    break;
                }
                let d = self
                    .get(name1, name2)
                    .ok_or::<Error>(Error::MissingDistance(name1.clone(), name2.clone()))?;
                output += &format!("  {}", d);
            }
            output += "\n"
        }

        Ok(output)
    }

    /// Outputs the matrix in phylip format to a String
    pub fn to_phylip(&self, square: bool) -> Result<String, Error> {
        if square {
            self.to_phylip_square()
        } else {
            self.to_phylip_triangle()
        }
    }

    /// Saves the tree to a newick file
    pub fn to_file(&self, path: &Path, square: bool) -> Result<(), Error> {
        match fs::write(path, self.to_phylip(square)?) {
            Ok(_) => Ok(()),
            Err(e) => Err(Error::IOError(e)),
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

        let mut matrix = Self::new(size);

        for (&n1, row) in names.iter().zip(rows) {
            for (&n2, dist) in names.iter().zip(row) {
                if let Some(old) = matrix.set(n1, n2, dist, true)? {
                    if old != dist {
                        return Err(ParseError::NonSymmetric(old, dist));
                    }
                }
            }
        }

        Ok(matrix)
    }

    /// Loads a tree from a newick file
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
        let mut matrix = DistanceMatrix::new(names.len());

        for pair in names.iter().combinations(2) {
            let (n1, d1) = pair[0];
            let (n2, d2) = pair[1];
            matrix.set(n1, n2, d1 * d2, false).unwrap();
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
        let square_nonsym = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";
        let mut matrix: Result<DistanceMatrix<f32>, _> =
            DistanceMatrix::from_phylip(square_nonsym, true);
        assert!(matrix.is_err());

        let err = matrix.err().unwrap();
        match err {
            ParseError::NonSymmetric(_, _) => {}
            _ => panic!("Error should be 'ParseError::NonSymmetric' not: {err}"),
        }

        let square_missing_dist = "4
s1    0  2  3  7
s2    2  0  6  10
s3    3  6  0
s5    5  10  15  0
";
        matrix = DistanceMatrix::from_phylip(square_missing_dist, true);
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
