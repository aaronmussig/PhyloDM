extern crate derive_more;

use derive_more::Add;
use derive_more::Sum;

#[derive(Debug, Default, Clone, Copy, Eq, PartialOrd, PartialEq, Hash)]
pub struct NodeId(pub usize);

#[derive(Debug, Default, Add, Copy, Clone, Sum)]
pub struct Edge(pub f64);

#[derive(Debug, Default, Clone, Copy, Eq, PartialEq, Hash, Add, Ord, PartialOrd)]
pub struct NodeDepth(pub usize);

#[derive(Debug, Eq, PartialEq, Hash, Clone, Ord, PartialOrd)]
pub struct Taxon(pub String);
