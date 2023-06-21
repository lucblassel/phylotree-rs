# PhyloTree

[![license_badge](https://img.shields.io/crates/l/phylotree)](https://choosealicense.com/licenses/gpl-3.0/) [![crate version](https://img.shields.io/crates/v/phylotree)](https://crates.io/crates/phylotree) [![rust doc](https://img.shields.io/docsrs/phylotree)](https://docs.rs/phylotree)

This crate aims to be a general purpose package to deal with phylogenetic trees and simulate them. Is also comes with a simple CLI tool to manipulate and simulate phylogenetic trees directly from the command line.

- [Crate usage](#using-this-crate)
- [CLI usage](#using-the-cli)
- [Python package](#python-package)

## Installing `phylotree`

### Crate

To use this crate just run `cargo add phylotree` in your cargo project or add `phylotree = "0.1.0"` to your `Cargo.toml` file.

### Cli

To install the CLI you can use
`cargo install phylotree`

Or you can build the project from source:

```shell
git clone https://github.com/lucblassel/phylotree.git
cd phylotree
cargo build --release
mv target/release/phylotree <somewhere/in/your/PATH>
```

## Using this crate

Below is some sample usage of the `phylotree` crate, please see [docs.rs/phylotree](https://docs.rs/phylotree) for the full documentation.

```rust
use phylotree::tree::Tree;

//////////////////////////////////
// Building a tree from scratch //
//////////////////////////////////
let mut tree = Tree::new();

// Add the root node
let root = tree.add(Node::new());

// Add a child to the root
let child1 = tree.add_child(Node::new_named("Child_1"), root, None).unwrap();
// Add a child to the root with a branch length
let child2 = tree.add_child(Node::new_named("Child_2"), root, Some(0.5)).unwrap();

// Add more children
let child3 = tree.add_child(Node::new_named("Child_3"), child1, None).unwrap();

// Get depth of child
assert_eq!(tree.get(&child3).unwrap().get_depth(), 2)

///////////////////////////////
// Reading and writing trees //
///////////////////////////////
let newick_str = "((A:0.1,B:0.2)F:0.6,(C:0.3,D:0.4)E:0.5)G;";
let tree = Tree::from_newick(newick_str).unwrap();

assert_eq!(tree.to_newick().unwrap(), newick_string)

//////////////////////
// Traversing trees //
//////////////////////
let newick_str = "((A,B)C,(D,E)F)G;";
let mut tree = Tree::from_newick(newick_str).unwrap();
let root = tree.get_root().unwrap();

let preorder: Vec<_> = tree.preorder(&root).unwrap()
    .iter()
    .map(|node_id| tree.get(node_id).unwrap().name.clone().unwrap())
    .collect();

assert_eq!(preorder, vec!["G", "C", "A", "B", "F", "D", "E"]);


/////////////////////
// Comparing trees //
/////////////////////

// The second tree is just a random rotation of the first,
// they represent the same phylogeney
let newick_orig = "((A:0.1,B:0.2)F:0.6,(C:0.3,D:0.4)E:0.5)G;";
let newick_rota = "((D:0.3,C:0.4)E:0.5,(B:0.2,A:0.1)F:0.6)G;";

let tree_orig = Tree::from_newick(newick_orig).unwrap();
let tree_rota = Tree::from_newick(newick_rota).unwrap();

let rf = tree_orig.robinson_foulds(&tree_rota).unwrap();

assert_eq!(rf, 0)


/////////////////////////////////
// Computing a distance matrix //
/////////////////////////////////
let newick = "((T3:0.2,T1:0.2):0.3,(T2:0.4,T0:0.5):0.6);";
let tree = Tree::from_newick(newick).unwrap();
// Compute the whole distance matrix
let matrix = tree.distance_matrix_recursive().unwrap();
let phylip="\
4
T0    0  1.6  0.9  1.6
T1    1.6  0  1.5  0.4
T2    0.9  1.5  0  1.5
T3    1.6  0.4  1.5  0
";

assert_eq!(matrix.to_phylip(true).unwrap(), phylip)
```

## Using the CLI

There is a simple CLI that comes with this package:

```
A simple command line tool to manipulate phylogenetic trees

Usage: phylotree <COMMAND>

Commands:
  generate     Generate random tree(s)
  stats        Get statistics about a tree
  compare      Compare two phylogenetic trees
  matrix       Output the phylogenetic distance matrix of the tree
  distance     Outputs a subset of phylogenetic distances
  collapse     Collapse branches that are under a certain branch length threshold
  remove       Remove tips from the trees
  deduplicate  Remove or collapse branches corresponding to identical sequences in a reference alignment
  help         Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

## Python package

A python package has also been implemented, using PyO3 to create the python bindings. It is available on [PyPi](https://pypi.org/project/phylotree/0.1.1/), you can see the API of the package [here](./phylotree/__init__.py).

Example usage:

```python
from phylotree import Tree

# Read a newick file
tree = Tree.from_newick("path/to/tree.nwk")

# Read a newick string
tree = Tree.from_string("((A,(C,E)D)B,((H)I)G)F;")

# Print newick formatted string of the tree
tree.to_newick()

# Get information about the tree
tree.n_nodes()
tree.n_tips()
tree.height()
tree.diameter()
tree.n_cherries()
tree.colless()
tree.sackin()
tree.is_binary
tree.is_rooted

# Get information about nodes
names = tree.get_leaf_names()
info = tree.get_node_attributes(name=names[0])
print(f"Parent edge has length: {info['parent_edge']}")
dist_branches, dist_topo = tree.get_distance(names=(names[0], names[1]))

distance_matrix = tree.to_matrix()
dist_mat = distance_matrix[(names[0], names[1])]
assert dist_branches == dist_mat

# Modify the tree
tree.compress() # Remove nodes with 1 parent and 1 child
tree.rescale() # Change branch lengths
tree.prune(name="D") # Remove sbutree rooted a specific node

# Traverse the tree in a given order
for node_id in tree.traversal(order="levelorder"):
    node = tree.get_node_attributes(id=node_id)
    print(node.name)
```
