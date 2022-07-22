use numpy::{PyArray2, ToPyArray};
use pyo3::{pymodule, PyResult, Python, types::PyModule};
use pyo3::prelude::*;

use crate::tree::ArenaTree as Tree;
use crate::types::{Edge, NodeId, Taxon};

#[pyclass]
struct PhyloDM {
    tree: Tree,
}

#[pymethods]
impl PhyloDM {
    #[new]
    fn new() -> Self {
        Self {
            tree: Tree::default(),
        }
    }

    pub fn load_from_newick_path(&mut self, path: String) {
        self.tree.load_from_newick_path(&path);
    }

    pub fn add_node(&mut self, taxon: Option<String>) -> usize {
        match taxon {
            Some(taxon) => self.tree.add_node(Some(Taxon(taxon))),
            None => self.tree.add_node(None),
        }.0
    }

    pub fn add_edge(&mut self, parent_id: usize, child_id: usize, length: f64) {
        self.tree.add_edge(NodeId(parent_id), NodeId(child_id), Edge(length));
    }

    pub fn get_nodes(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for node in &self.tree.nodes {
            out.push(node.id.0);
        }
        out
    }

    pub fn dm(&mut self, norm: bool) -> Py<PyArray2<f64>> {
        let array = self.tree.dm(norm);
        Python::with_gil(|py| {
            return Py::from(array.to_pyarray(py));
        })
    }

    pub fn taxa(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for node_id in &self.tree.row_idx_to_leaf_idx {
            let node = self.tree.get_node(*node_id);
            out.push(node.taxon.clone().unwrap().0);
        }
        out
    }

    pub fn length(&self) -> f64 {
        self.tree.length().0
    }

}

#[pymodule]
fn pdm(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PhyloDM>()?;
    Ok(())
}
