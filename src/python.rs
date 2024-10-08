use numpy::{PyArray1, PyArray2, PyArrayMethods, ToPyArray};
use pyo3::{Py, pyclass, pymethods, pymodule, PyResult, Python, types::PyModule, Bound};
use pyo3::exceptions::PyValueError;

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

    pub fn load_from_newick_path(&mut self, path: &str) -> PyResult<()> {
        let result = self.tree.load_from_newick_path(path);
        if result.is_err() {
            return Err(PyValueError::new_err("Unable to load newick."));
        }
        Ok(())
    }

    pub fn add_node(&mut self, taxon: Option<&str>) -> PyResult<usize> {
        let out = match taxon {
            Some(taxon) => self.tree.add_node(Some(&Taxon(taxon.to_string()))),
            None => self.tree.add_node(None),
        };
        if out.is_err() {
            return Err(PyValueError::new_err("Unable to add node."));
        }
        Ok(out.unwrap().0)
    }

    pub fn add_edge(&mut self, parent_id: usize, child_id: usize, length: f64) {
        self.tree.add_edge(
            NodeId(parent_id),
            NodeId(child_id),
            Edge(length),
        );
    }

    pub fn update_edge_lengths(&mut self, child_nodes: &Bound<'_, PyArray1<usize>>, lengths: &Bound<'_, PyArray1<f64>>) -> PyResult<()> {
        
        let binding = lengths.to_vec().unwrap();
        let new_lengths_vec: Vec<Edge> = binding.iter().map(|x| Edge(*x)).collect();
        
        let child_nodes_binding = child_nodes.to_vec().unwrap();
        let child_nodes_vec: Vec<NodeId> = child_nodes_binding.iter().map(|x| NodeId(*x)).collect();
        
        let result = self.tree.update_edge_lengths(&child_nodes_vec, &new_lengths_vec);
        
        if result.is_err() {
            return Err(PyValueError::new_err("Unable to update edge lengths."));
        }
        
        Ok(())
    }
    
    pub fn update_all_edge_lengths(&mut self, length: f64) -> PyResult<()> {
        let result = self.tree.update_all_edge_lengths(Edge(length));
        if result.is_err() {
            return Err(PyValueError::new_err("Unable to update all edge lengths."));
        }
        Ok(())
    }

    pub fn get_nodes(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for node in &self.tree.nodes {
            out.push(node.id.0);
        }
        out
    }

    pub fn dm(&mut self, norm: bool) -> PyResult<Py<PyArray2<f64>>> {
        let matrix = self.tree.matrix(norm);
        if matrix.is_err() {
            return Err(PyValueError::new_err("Unable to compute distance matrix."));
        }
        let (_, array) = matrix.unwrap();
        Ok(Python::with_gil(|py| {
            return Py::from(array.to_pyarray_bound(py));
        }))
    }

    pub fn taxa(&self) -> PyResult<Vec<String>> {
        let mut out: Vec<String> = Vec::new();
        let taxa = self.tree.leaf_nodes();
        if taxa.is_err() {
            return Err(PyValueError::new_err("Unable to get taxa."));
        }
        for taxon in taxa.unwrap() {
            out.push(taxon.0);
        }
        Ok(out)
    }

    pub fn length(&self) -> f64 {
        self.tree.length().0
    }

    pub fn compute_row_vec(&mut self) -> PyResult<()> {
        if self.tree.compute_row_vec().is_err() {
            return Err(PyValueError::new_err("Unable to compute row vector."));
        }
        Ok(())
    }

    pub fn distance(&self, a: &str, b: &str, norm: bool) -> f64 {
        let taxon_a = Taxon(a.to_string());
        let taxon_b = Taxon(b.to_string());
        self.tree.distance(&taxon_a, &taxon_b, norm)
    }
    
    pub fn get_nearest_taxa(&self, taxon: &str) -> Vec<String> {
        let taxon = Taxon(taxon.to_string());
        let result = self.tree.get_nearest_taxa(&taxon);
        let mut out: Vec<String> = Vec::with_capacity(result.len());
        for taxon in result {
            out.push(taxon.0.clone());
        }
        out
    }
    
}

#[pymodule]
fn pdm(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PhyloDM>()?;
    Ok(())
}
