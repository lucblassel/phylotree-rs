use clap::{Parser, Subcommand};
use clap_complete::Shell;
use std::path::PathBuf;

use phylotree::distr::Distr;
use phylotree::TreeShape;

const ABOUT: &str = "A Simple command line tool to manipulate phylogenetic trees";

/// A simple command line tool to manipulate phylogenetic trees
#[derive(Parser, Debug)]
#[command(author, version, long_about=ABOUT)]
pub struct Args {
    #[command(subcommand)]
    /// The command to execute
    pub command: Commands,
}

/// The available commands in the `phylotree` tool
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Generate random tree(s)
    Generate {
        /// Number of tips in the generated tree
        #[arg(short, long, default_value_t = 20)]
        tips: usize,

        /// Generate uniform branch lengths
        #[arg(short, long)]
        branch_lengths: bool,

        /// Output file (directory if generating multiple trees)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Number of trees to generate
        #[arg(short = 'n', long)]
        trees: Option<usize>,

        /// Shape of the tree
        #[arg(value_enum, short, long, default_value_t=TreeShape::Yule)]
        shape: TreeShape,

        /// Distribution of branch lengths
        #[arg(value_enum, short, long, default_value_t=Distr::Uniform)]
        distribution: Distr,
    },

    /// Get statistics about a tree
    Stats {
        /// Input newick file of the tree
        trees: Vec<PathBuf>,
    },

    /// Compare phylogenetic trees to a reference
    ///
    /// This will return:
    ///  - the path of the compared tree
    ///  - the number of bipartitions unique to the reference tree
    ///  - the number of common bipartitions
    ///  - the number of bipartitions unique to the compared tree
    ///  - the Robinson-Foulds distance
    ///  - the normalize Robinson-Foulds distance
    ///  - the weighted Robinson-Foulds distance
    ///  - the Khuner-Felsenstein branch-score
    #[clap(verbatim_doc_comment)]
    Compare {
        /// Reference tree
        reftree: PathBuf,
        /// Tree(s) to compare to reference
        tocompare: Vec<PathBuf>,
    },
    /// Output the phylogenetic distance matrix of the tree
    Matrix {
        /// The phylogenetic tree
        tree: PathBuf,
        /// Output a square matrix instead of a triangular one
        #[arg(short, long)]
        square: bool,
        /// File to save the matrix to
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Outputs a subset of phylogenetic distances
    Distance {
        /// The phylogenetic tree
        tree: PathBuf,
        /// The tips to consider
        tips: Vec<String>,
        /// Tab separated file to save distances to
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Collapse branches that are under a certain branch length threshold
    Collapse {
        /// The phylogenetic tree
        tree: PathBuf,
        /// The threshold under which to collapse the branch
        threshold: phylotree::tree::EdgeLength,
        /// Print the number of collapsed branches at the end
        #[arg(short, long)]
        verbose: bool,
        /// Don't collapse external branches
        #[arg(short, long)]
        exclude_tips: bool,
        /// File to save the tree to
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Randomly resolve a multifurcated tree to a binary one
    Resolve {
        /// The phylogenetic tree
        tree: PathBuf,
        /// The output tree
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Remove tips from the trees
    Remove {
        /// The phylogenetic tree
        tree: PathBuf,
        /// Names of tips to remove
        tips: Vec<String>,
        /// File to save the tree to
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Remove or collapse branches corresponding to identical sequences in a reference alignment
    Deduplicate {
        /// The phylogenetic tree
        tree: PathBuf,
        /// The FASTA alignment
        alignment: PathBuf,
        /// Whether to collapse duplicate branches to 0 instead of removing one of them
        #[arg(short, long)]
        collapse: bool,
        /// File to save the tree to
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Print the number of collapse/removed branches at the end
        #[arg(short, long, action=clap::ArgAction::Count)]
        verbose: u8,
    },
    /// Draw an image of the phylogenetic tree
    Draw {
        /// The phylogenetic tree
        tree: PathBuf,
        /// Width of the output svg
        #[arg(short = 'W', long, default_value_t = 1000.)]
        width: f64,
        /// height of the output svg
        #[arg(short = 'H', long, default_value_t = 1000.)]
        height: f64,
        /// Remove the white background from the svg
        #[arg(short, long)]
        transparent: bool,
        /// Padding as percentage of min(width, height)
        #[arg(short, long, default_value_t = 5.)]
        padding: f64,
        /// Output svg file
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Generate shell-completion scripts
    Completion {
        /// Shell
        #[arg(value_enum)]
        shell: Shell,
    },
    /// Rescale branch lengths of tree(s)
    Rescale {
        /// Rescaling factor by which branch lengtsh will be multiplied
        /// If there are multiple input trees, then it will be treated as a directory
        factor: f64,
        /// Path to tree(s) to rescale
        trees: Vec<PathBuf>,
        /// Path to output.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

impl Commands {}
