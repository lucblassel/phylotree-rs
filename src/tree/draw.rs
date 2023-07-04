///! A module to draw phylogenetic trees
use serde::Serialize;

use super::{NodeError, Tree, TreeError};
use std::f64::consts::PI;

#[derive(Serialize)]
/// Represents a branch in a visualization as a line between 2 points
pub struct Branch {
    /// X coordinate of the start of the branch
    pub xstart: f64,
    /// Y coordinate of the start of the branch
    pub ystart: f64,
    /// X coordinate of the end of the branch
    pub xend: f64,
    /// Y coordinate of the end of the branch
    pub yend: f64,
}

#[derive(Serialize)]
/// Represents a node in a visualization as a point and a label
pub struct Node {
    /// X coordinate of the node
    pub x: f64,
    /// Y coordinate of the node
    pub y: f64,
    /// Label of the node
    pub label: Option<String>,
}

/// Represents the whole layout of a visualized tree
#[derive(Default, Serialize)]
pub struct Layout {
    /// Branches of the tree
    pub branches: Vec<Branch>,
    /// Nodes of the tree
    pub nodes: Vec<Node>,
}

impl Branch {
    /// Rescale a branch by applying a multiplicative
    /// factor to it's start and end coordinates
    pub fn rescale(&mut self, factor: f64) {
        self.xstart *= factor;
        self.ystart *= factor;
        self.xend *= factor;
        self.yend *= factor;
    }
}

impl Node {
    /// Rescale a node position by applying a multiplicative
    /// factor to it's start and end coordinates
    pub fn rescale(&mut self, factor: f64) {
        self.x *= factor;
        self.y *= factor;
    }
}

impl Layout {
    /// Rescale a Layout by applying a multiplicative factor
    pub fn rescale(&mut self, factor: f64) {
        self.branches.iter_mut().map(|b| b.rescale(factor)).count();
        self.nodes.iter_mut().map(|n| n.rescale(factor)).count();
    }
}

/// Returns the radial layout of a tree according to the algorithm
/// detailed in:  
/// <https://link.springer.com/chapter/10.1007/11602613_110>
pub fn radial_layout(tree: &Tree) -> Result<Layout, TreeError> {
    let mut l = vec![0.; tree.size()];
    let mut w = vec![2. * PI; tree.size()];
    let mut t = vec![0.; tree.size()];
    let mut x = vec![(0., 0.); tree.size()];

    let root = tree.get_root()?;

    // Get subtree sizes
    for v in tree.postorder(&root)?.iter() {
        let node = tree.get(v)?;
        if node.is_tip() {
            l[*v] = 1 as f64
        } else {
            for w_ in node.children.iter() {
                l[*v] += l[*w_];
            }
        }
    }

    // Compute angles
    for v in tree.preorder(&root)?.iter() {
        let node = tree.get(v)?;
        if *v != root {
            let d = node.parent_edge.ok_or(TreeError::MissingBranchLengths)?;
            let u = node
                .parent
                .ok_or(TreeError::from(NodeError::HasNoParent(*v)))?;

            x[*v] = (
                x[u].0 + d * f64::cos(t[*v] + w[*v] / 2.),
                x[u].1 + d * f64::sin(t[*v] + w[*v] / 2.),
            );
        }

        let mut n = t[*v];

        for w_ in node.children.iter() {
            w[*w_] = (l[*w_] * 2. * PI) / l[root];
            t[*w_] = n;
            n += w[*w_];
        }
    }

    let mut layout = Layout::default();

    for v in tree.preorder(&root)?.iter().skip(1) {
        let u = tree
            .get(v)?
            .parent
            .ok_or(TreeError::from(NodeError::HasNoParent(*v)))?;
        let name = tree.get(v)?.name.clone();

        let line_start = x[u];
        let line_end = x[*v];

        layout.branches.push(Branch {
            xstart: line_start.0,
            ystart: line_start.1,
            xend: line_end.0,
            yend: line_end.1,
        });

        layout.nodes.push(Node {
            x: line_end.0,
            y: line_end.1,
            label: name,
        });
    }

    Ok(layout)
}
