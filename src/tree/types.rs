extern crate derive_more;

use derive_more::Add;
use derive_more::Sum;

/// The unique identifier corresponding to a node in the tree.
///
/// # Examples
///
/// ```
/// use phylodm::tree::NodeId;
/// NodeId(0);
/// ```
#[derive(Debug, Default, Clone, Copy, Eq, PartialOrd, PartialEq, Hash, Ord)]
pub struct NodeId(pub usize);

/// A distance between two nodes in the tree.
///
/// # Examples
///
/// ```
/// use phylodm::tree::Edge;
/// Edge(1.0);
/// ```
#[derive(Debug, Default, Add, Copy, Clone, Sum)]
pub struct Edge(pub f64);

/// The depth of a node within the tree, root node is 0.
///
/// # Examples
///
/// ```
/// use phylodm::tree::NodeDepth;
/// NodeDepth(0);
/// ```
#[derive(Debug, Default, Clone, Copy, Eq, PartialEq, Hash, Add, Ord, PartialOrd)]
pub struct NodeDepth(pub usize);

/// A taxon/leaf node within the tree.
///
/// # Examples
///
/// ```
/// use phylodm::tree::Taxon;
/// Taxon("A".to_string());
/// ```
#[derive(Debug, Eq, PartialEq, Hash, Clone, Ord, PartialOrd)]
pub struct Taxon(pub String);
