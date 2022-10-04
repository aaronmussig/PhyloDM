use numpy::{PyArray2, ToPyArray};
use pyo3::{Py, pyclass, pymethods, pymodule, PyResult, Python, types::PyModule};

use crate::pdm::PDM as RustPhyloDM;
use crate::tree::{Edge, NodeId, Taxon};

#[pyclass]
struct PhyloDM {
    tree: RustPhyloDM,
}

#[pymethods]
impl PhyloDM {
    #[new]
    fn new() -> Self {
        Self {
            tree: RustPhyloDM::default(),
        }
    }

    pub fn load_from_newick_path(&mut self, path: &str) {
        self.tree.load_from_newick_path(path);
    }

    pub fn add_node(&mut self, taxon: Option<&str>) -> usize {
        match taxon {
            Some(taxon) => self.tree.add_node(Some(&Taxon(taxon.to_string()))),
            None => self.tree.add_node(None),
        }
            .0
    }

    pub fn add_edge(&mut self, parent_id: usize, child_id: usize, length: f64) {
        self.tree.add_edge(
            NodeId(parent_id),
            NodeId(child_id),
            Edge(length),
        );
    }

    pub fn get_nodes(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for node in &self.tree.nodes {
            out.push(node.id.0);
        }
        out
    }

    pub fn dm(&mut self, norm: bool) -> Py<PyArray2<f64>> {
        let (_, array) = self.tree.matrix(norm);
        Python::with_gil(|py| {
            return Py::from(array.to_pyarray(py));
        })
    }

    pub fn taxa(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for taxon in self.tree.leaf_nodes() {
            out.push(taxon.0);
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
