use thiserror::Error;

#[derive(Error, Debug)]
pub enum TreeError {
    #[error("This tree is not Binary.")]
    IsNotBinary,
    #[error("This tree is not rooted.")]
    IsNotRooted,
    #[error("This tree is empty.")]
    IsEmpty,
    #[error("All your leaf nodes must ne named.")]
    UnnamedLeaves,
    #[error("Your leaf names must be unique.")]
    DuplicateLeafNames,
    #[error("The leaf index of the tree is not initialized.")]
    LeafIndexNotInitialized,
    #[error("The tree must have all branch lengths")]
    MissingBranchLengths,
    #[error("The trees have different tips indices")]
    DifferentTipIndices
}

#[derive(Error, Debug, PartialEq)]
pub enum ParseError {
    #[error("Cannot have whitespace in number field")]
    WhiteSpaceInNumber,
    #[error("Missing a closing bracket")]
    UnclosedBracket,
    #[error("The tree is missin a semi colon at the end")]
    NoClosingSemicolon
}