use std::{
    collections::{HashMap, HashSet},
    fmt::{Display, Debug},
};

use itertools::Itertools;
use num_traits::{Float, Zero};

use crate::errors::DistanceMatrixError;

type Error = Box<dyn std::error::Error>;
type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
/// Struct representing a phylogenetic distance matrix
pub struct DistanceMatrix<T> {
    pub size: usize,
    pub ids: HashSet<String>,
    matrix: HashMap<(String, String), T>,
}

impl<T> DistanceMatrix<T>
where
    T: Display + Debug + Float + Zero,
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

    /// Set an entry in the distance matrix
    pub fn set(&mut self, id_1: &str, id_2: &str, dist: T, overwrite: bool) -> Result<()> {
        if let Some(old_dist) = self.get_mut(id_1, id_2) {
            if overwrite {
                *old_dist = dist;
                return Ok(());
            } else {
                return Err(DistanceMatrixError::OverwritingNotPermitted.into());
            }
        }

        self.matrix.insert(self.get_key(id_1, id_2), dist);
        self.ids.insert(id_1.to_owned());
        self.ids.insert(id_2.to_owned());

        if self.ids.len() > self.size {
            Err(DistanceMatrixError::SizeExceeded.into())
        } else {
            Ok(())
        }
    }

    /// Returns a string representing the distance matrix in square format
    fn to_phylip_square(&self) -> Result<String> {
        let names: Vec<_> = self.ids.iter().cloned().sorted().collect();
        let mut output = format!("{}\n", self.size);

        for name1 in names.iter() {
            output += &format!("{name1}  ");
            for name2 in names.iter() {
                let d = if name1 != name2 {
                    *self.get(name1, name2).ok_or::<Error>(
                        DistanceMatrixError::MissinDinstance(name1.clone(), name2.clone()).into(),
                    )?
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
    fn to_phylip_triangle(&self) -> Result<String> {
        let names: Vec<_> = self.ids.iter().cloned().sorted().collect();
        let mut output = format!("{}\n", self.size);

        for name1 in names.iter() {
            output += &format!("{name1}  ");
            for name2 in names.iter() {
                if name1 == name2 {
                    break;
                }
                let d = self.get(name1, name2).ok_or::<Error>(
                    DistanceMatrixError::MissinDinstance(name1.clone(), name2.clone()).into(),
                )?;
                output += &format!("  {}", d);
            }
            output += "\n"
        }

        Ok(output)
    }

    /// Outputs the matrix in phylip format to a String
    pub fn to_phylip(&self, square: bool) -> Result<String> {
        if square {
            self.to_phylip_square()
        } else {
            self.to_phylip_triangle()
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_to_phylip() {
        let names = vec![("s1", 1.0), ("s2", 2.0), ("s3", 3.0), ("s5", 5.0)];
        let square = "4
s1    0  2  3  5
s2    2  0  6  10
s3    3  6  0  15
s5    5  10  15  0
";

        let triangle = "4
s1  
s2    2
s3    3  6
s5    5  10  15
";
        let mut matrix = DistanceMatrix::new(names.len());

        for pair in names.iter().combinations(2) {
            let (n1, d1) = pair[0];
            let (n2, d2) = pair[1];
            matrix.set(n1, n2, d1 * d2, false).unwrap()
        }

        assert_eq!(square, matrix.to_phylip(true).unwrap());
        assert_eq!(triangle, matrix.to_phylip(false).unwrap());
    }
}
