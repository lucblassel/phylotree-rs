use crate::tree::Node as RustNode;
use crate::tree::Tree as RustTree;
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{IntoPyDict, PyDict},
    Python,
};
use std::{collections::HashMap, path::Path};

#[pyclass]
struct Tree {
    tree: RustTree,
}

#[pyclass]
struct Node {
    node: RustNode,
}

struct TreeError(crate::tree::TreeError);
struct ParseError(crate::tree::NewickParseError);
struct NodeError(crate::tree::NodeError);

impl From<TreeError> for PyErr {
    fn from(err: TreeError) -> Self {
        PyValueError::new_err(err.0.to_string())
    }
}

impl From<crate::tree::TreeError> for TreeError {
    fn from(err: crate::tree::TreeError) -> Self {
        Self(err)
    }
}

impl From<ParseError> for PyErr {
    fn from(err: ParseError) -> Self {
        PyValueError::new_err(err.0.to_string())
    }
}

impl From<crate::tree::NewickParseError> for ParseError {
    fn from(err: crate::tree::NewickParseError) -> Self {
        Self(err)
    }
}

impl From<NodeError> for PyErr {
    fn from(err: NodeError) -> Self {
        PyValueError::new_err(err.0.to_string())
    }
}

impl From<crate::tree::NodeError> for NodeError {
    fn from(err: crate::tree::NodeError) -> Self {
        Self(err)
    }
}

#[pymethods]
impl Tree {
    #[new]
    fn new() -> Self {
        Self {
            tree: RustTree::new(),
        }
    }

    #[staticmethod]
    pub fn from_newick(path: &str) -> Result<Self, ParseError> {
        let path = Path::new(path);
        let tree = RustTree::from_file(path)?;

        Ok(Tree { tree })
    }

    #[staticmethod]
    pub fn from_string(string: &str) -> Result<Self, ParseError> {
        let tree = RustTree::from_newick(string)?;
        Ok(Tree { tree })
    }

    fn to_string(&self) -> Result<String, TreeError> {
        let s = self.tree.to_newick()?;
        Ok(s)
    }

    fn to_file(&self, path: &str) -> Result<(), TreeError> {
        let path = Path::new(path);
        self.tree.to_file(&path)?;
        Ok(())
    }

    fn is_binary(&self) -> bool {
        self.tree.is_binary()
    }

    fn is_rooted(&self) -> bool {
        self.tree.is_rooted().unwrap_or(false)
    }

    fn height(&self) -> Option<f64> {
        self.tree.height()
    }

    fn n_tips(&self) -> usize {
        self.tree.n_leaves()
    }

    fn n_nodes(&self) -> usize {
        self.tree.size()
    }

    fn diameter(&self) -> Option<f64> {
        self.tree.diameter()
    }

    fn n_cherries(&self) -> Result<usize, TreeError> {
        let cherries = self.tree.cherries()?;
        Ok(cherries)
    }

    fn colless(&self, normalisation: Option<&str>) -> Result<f64, TreeError> {
        let res = match normalisation {
            Some("yule") => self.tree.colless_yule()?,
            Some("pda") => self.tree.colless_pda()?,
            None => self.tree.colless()? as f64,
            _ => unreachable!(),
        };
        Ok(res)
    }

    fn sackin(&self, normalisation: Option<&str>) -> Result<f64, TreeError> {
        let res = match normalisation {
            Some("yule") => self.tree.sackin_yule()?,
            Some("pda") => self.tree.sackin_pda()?,
            None => self.tree.sackin()? as f64,
            _ => unreachable!(),
        };
        Ok(res)
    }

    fn compare(&self, other: &Self) -> Result<HashMap<&'static str, f64>, TreeError> {
        let cmp = self.tree.compare_topologies(&other.tree)?;
        Ok(HashMap::from_iter([
            ("rf", cmp.rf),
            ("norm_rf", cmp.norm_rf),
            ("weighted_rf", cmp.weighted_rf),
            ("branch_score", cmp.branch_score),
        ]))
    }

    fn compress(&mut self) -> Result<(), TreeError> {
        self.tree.compress()?;
        Ok(())
    }

    fn rescale(&mut self, factor: f64) {
        self.tree.rescale(factor);
    }

    fn get_name_index(&self, name: &str) -> Option<usize> {
        let node = self.tree.get_by_name(name);
        node.map(|n| n.id)
    }

    fn prune(&mut self, root: usize) -> Result<(), TreeError> {
        self.tree.prune(&root)?;
        Ok(())
    }

    fn copy(&self) -> Self {
        let tree = self.tree.clone();
        Self { tree }
    }

    fn get_distance(
        &self,
        source: usize,
        target: usize,
    ) -> Result<(Option<f64>, usize), TreeError> {
        let dist = self.tree.get_distance(&source, &target)?;
        Ok(dist)
    }

    fn get_root_id(&self) -> Result<usize, TreeError> {
        let id = self.tree.get_root()?;
        Ok(id)
    }

    fn inorder(&self, root: usize) -> Result<Vec<usize>, TreeError> {
        let traversal = self.tree.inorder(&root)?;
        Ok(traversal)
    }

    fn postorder(&self, root: usize) -> Result<Vec<usize>, TreeError> {
        let traversal = self.tree.postorder(&root)?;
        Ok(traversal)
    }

    fn preorder(&self, root: usize) -> Result<Vec<usize>, TreeError> {
        let traversal = self.tree.preorder(&root)?;
        Ok(traversal)
    }

    fn levelorder(&self, root: usize) -> Result<Vec<usize>, TreeError> {
        let traversal = self.tree.levelorder(&root)?;
        Ok(traversal)
    }

    fn get_leaf_names(&self) -> Vec<Option<String>> {
        self.tree.get_leaf_names()
    }

    fn get_node_attributes(&self, id: usize, py: Python) -> Result<PyObject, TreeError> {
        let node = self.tree.get(&id)?;
        let name = node.name.clone();
        let parent: Option<usize> = node.parent.clone();
        let children: Vec<usize> = node.children.clone();
        let parent_edge: Option<f64> = node.parent_edge.clone();
        let comment: Option<String> = node.comment.clone();

        let mut key_vals: Vec<(&str, PyObject)> = vec![("id", id.to_object(py))];

        key_vals.push(("name", name.map(|n| n.to_object(py)).unwrap_or(py.None())));
        key_vals.push((
            "parent",
            parent.map(|p| p.to_object(py)).unwrap_or(py.None()),
        ));
        key_vals.push((
            "parent_edge",
            parent_edge.map(|e| e.to_object(py)).unwrap_or(py.None()),
        ));
        key_vals.push((
            "comment",
            comment.map(|c| c.to_object(py)).unwrap_or(py.None()),
        ));
        key_vals.push(("children", children.to_object(py)));

        let dict = key_vals.into_py_dict(py);

        Ok(dict.into())
    }
}

#[pymethods]
impl Node {
    #[new]
    fn new(name: Option<&str>) -> Self {
        let node = match name {
            Some(name) => RustNode::new_named(name),
            None => RustNode::new(),
        };

        Node { node }
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn pytree(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Tree>()?;
    Ok(())
}
