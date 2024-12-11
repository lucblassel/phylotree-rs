# Phylotree changelog

## Unreleased

## v0.1.3 - 2024-12-11
### Added
#### Lib
- Adding branch depth information to bipartitions
- Basic nexus file output
- Adding basic CI (`cargo fmt`, `cargo test` and `cargo clippy`)
- Iterators for DistanceMatrix
- min and max methods for DistanceMatrix
- Node merging into parent
- Tree building with UPGMA
#### Bin
- Shell completion 
- Rescale command
- Many to single ref tree comparison
### Changed 
- Better output for `phylotree stats` command.
- More informative type names
- Better output for criterion benchmarks
### Fixed
- Make node searching possible with capturing closures (#9 and #10 🙏 @4less)
- Fixing bad tests
- Triangular Matrix / Vector index conversions
- Removed un-needed allocations and clones to speed up execution of [`phylocompare`](https://github.com/lucblassel/phylocompare) (#12 🙏 @krtab)

## v0.1.2 - 2023-09-20
### Added
- Basic Python bindings to this crate usign PyO3. Behind `python` feature. Used for the [`phylotree` python package](https://pypi.org/project/phylotree/).
- Adding Yule process to tree simulation
- Adding drawing method to visualize trees as svgs. Only radial layout for now. 
- Added part of the ETE3 test suite
- Added possibility to search for nodes with a user specified closure. 
- Possibility to output newick with different formats *(with/without comments, no lengths, only leaf names, ...)* with `tree::Tree::to_formatted_newick`
- Random resolution of polytomies
- Added ladderisation function
- Added commands to CLI:
    - Collapse: set branches under a given threshold to 0
    - Resolve: randomly resolve multifurcations
    - Remove: remove tips from the tree
    - Deduplicate: given a sequence alignment, remove tips corresponding to duplicated sequences
    - Draw: draw a tree to a svg file

### Changed
- More efficient method to compute pairwise distances. Inspired by the `PhyloDM` crate.
- Added new simulation methods to `generate` command of CLI.

## V0.1.1 - 2023-05-06
### Changed
Better documentation and README

## v0.1.0 - 2023-05-02 [Initial Release]
### Added
- Basic `Tree`, `Node` and `DistanceMatrix` structs to deal with phylogenetic trees
- Simple manipulation methods on a tree: 
    - Add nodes to trees and build them programatically
    - Retreive nodes from trees and search for nodes with certain properties (e.g. get root, or all leaf nodes, ...)
    - Get a subtree rooted at a given node
- Tree traversals
- Tree characteristics *(is binary?, is rooted?, ...)* and metrics *(colless index, ...)*
- Tree comparison methods: 
    - Robinson-Foulds (+ normalized and weighted variants)
    - Khuner-Felsenstein
- Methods to find paths within the tree, i.e. last common ancestors, distance betwee 2 nodes, ...
- Methods to compute the whole phylogenetic distance matrix efficiently.
- Methods to modify the tree: prune, rescale, compress. 
- Methods to read/write trees from/to newick files or strings.
- Methods to simulate trees of different shapes with simple processes: Caterpillar or Uniform *(like ete3's populate)*
- Simple CLI that uses these methods with following commands:
    - generate: generate random trees
    - stats: get statistics of a tree
    - compare: compare 2 trees
    - matrix: get distance matrix from tree

