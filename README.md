# PhyloTree

This crate aims to be a general purpose crate to deal with phylogenetic trees and simulate them. 

## Installing
To install this tool you must build it from source with Cargo:
```shell
git clone https://github.com/lucblassel/phylotree.git
cd phylotree
cargo build --release
mv target/release/phylotree <somewhere/in/your/PATH>
```

## CLI
There is a simple CLI that comes with this package:

```
A simple command line tool to generate a random phylogenetic tree

Usage: phylotree <COMMAND>

Commands:
  generate  Generate random tree(s)
  stats     Get statistics about a tree
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

### phylotree generate
```
Generate random tree(s)

Usage: phylotree generate [OPTIONS] --output <OUTPUT>

Options:
  -t, --tips <TIPS>      Number of tips in the generated tree [default: 20]
  -b, --branch-lengths   Generate uniform branch lengths
  -o, --output <OUTPUT>  Output file (directory if generating multiple trees)
  -n, --trees <TREES>    Number of trees to generate
  -h, --help             Print help
```

### phylotree stats
```
Get statistics about a tree

Usage: phylotree stats [TREES]...

Arguments:
  [TREES]...  Input newick file of the tree

Options:
  -h, --help  Print help
```

### phylotree compare
```
Compare two trees

Usage: phylotree compare <REFTREE> <TOCOMPARE>

Arguments:
  <REFTREE>    Reference tree
  <TOCOMPARE>  Tree to compare to reference

Options:
  -h, --help  Print help
```

## API
Rough API ref goes here

## RoadMap
Various goals and aims for this project:
 - [ ] Tree object
    - [x] Create base tree and node objects
    - [x] Enforce tip name uniqueness for distance calculations
    - [x] Rescale tree branch lengths
    - [x] Traversals:
        - [x] Pre-order
        - [x] Post-order
        - [x] Level-order
        - [x] In-order (for binary trees only)
    - [ ] Distances:
        - [x] LCA
        - [x] Distance between nodes
        - [x] RF distance between trees
        - [x] Weighted RF 
        - [x] Branch score index
        - [ ] Transfer distance
        - [ ] Other distances ? 
    - [x] Metrics:
        - [x] Height of the tree
        - [x] Diameter of the tree
        - [x] Sackin index of the tree (with normalized variants)
        - [x] Colless index of the tree
        - [x] Number of cherries
    - [ ] Convert to distance matrix
        - [x] Write to phylip
        - [ ] Read form phylip 
    - [ ] Build from Distance matrix
        - [ ] UPGMA
        - [ ] NJ
        - [ ] BIONJ ?
 - [ ] Tree simulation
    - [x] Random binary topology
    - [x] Caterpillar trees
    - [x] Add different branch length distributions
    - [ ] Add birth death simulation
    - [ ] Add coalescent simulation
    - [ ] Simulate ultrametric trees
 - [x] I/O stuff:
    - [x] Write tree to newick
    - [x] Read tree from newick *(Could be better though...)*