use itertools::Itertools;

use super::{Tree, TreeError};

#[derive(Debug, Clone)]
pub struct TreeSession<'a> {
    tree: &'a Tree,
    taxa: Vec<String>,
}

// Create from a tree
impl<'a> TryFrom<&'a Tree> for TreeSession<'a> {
    type Error = TreeError;

    fn try_from(tree: &'a Tree) -> Result<Self, Self::Error> {
        let taxa = tree
            .get_leaf_names()
            .into_iter()
            .flatten()
            .sorted()
            .collect_vec();

        if taxa.len() != tree.n_leaves() {
            Err(TreeError::UnnamedLeaves)
        } else if taxa.iter().duplicates().count() > 0 {
            Err(TreeError::DuplicateLeafNames)
        } else {
            Ok(Self { tree, taxa })
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const NW1: &str = "(((((((((Tip9,Tip8),Tip7),Tip6),Tip5),Tip4),Tip3),Tip2),Tip1),Tip0);";
    const NW2: &str = "(Tip0,((Tip2,(Tip3,(Tip4,(Tip5,(Tip6,((Tip8,Tip9),Tip7)))))),Tip1));";

    #[test]
    fn session_create() {
        let tree = Tree::from_newick(NW1).unwrap();
        let session = TreeSession::try_from(&tree).unwrap();
        eprintln!("{session:?}");

        panic!()
    }
}
