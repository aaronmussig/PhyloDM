#![crate_name = "phylodm"]

// #[cfg(any(feature = "python"))]
// extern crate pyo3;

#[cfg(any(feature = "python"))]
pub mod python;

pub mod tree;
pub mod node;
pub mod types;
pub mod util;
