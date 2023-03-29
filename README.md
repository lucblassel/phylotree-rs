# PhyloTree

This crate aims to be a general purpose crate to deal with phylogenetic trees and simulate them. 

## RoadMap
Various goals and aims for this project:
 - [ ] Tree object
    - [x] Create base tree and node objects
    - [ ] Traversals:
        - [x] Pre-order
        - [x] Post-order
        - [x] Level-order
        - [ ] How to do inorder with more than 2 child nodes
    - [ ] Distances:
        - [x] LCA
        - [x] Distance between nodes
        - [ ] Height of the tree
        - [ ] RF distance between trees
        - [ ] Transfer distance
        - [ ] Other distances ? 
 - [ ] Tree simulation
    - [x] Random binary topology
    - [ ] Add different branch length distributions
    - [ ] Add birth death simulation
    - [ ] Add coalescent simulation
    - [ ] Simulate ultrametric trees
 - [ ] I/O stuff:
    - [x] Write tree to newick
    - [ ] Read tree from newick


## API
Rough API ref goes here