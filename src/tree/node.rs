use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::tree::{Edge, NodeDepth, NodeId, Taxon};

/// A node in the tree.
#[derive(Default)]
pub struct Node {
    pub id: NodeId,
    pub taxon: Option<Taxon>,
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
    pub parent_distance: Option<Edge>,
    pub depth: Option<NodeDepth>,
    pub desc_distances: Option<HashMap<NodeId, Edge>>,
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for Node {}

impl Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl Node {
    #[must_use]
    pub fn new(id: NodeId, taxon: Option<Taxon>) -> Self {
        Self {
            id,
            taxon,
            parent: None,
            children: vec![],
            parent_distance: None,
            depth: None,
            desc_distances: None,
        }
    }

    /// Add a child to this node.
    pub fn add_child(&mut self, child: NodeId) {
        self.children.push(child);
    }

    /// Set the parent node.
    pub fn set_parent(&mut self, parent: NodeId, length: Edge) {
        self.parent = Some(parent);
        self.parent_distance = Some(length);
    }

    /// Set the depth of this node.
    pub fn set_depth(&mut self, depth: NodeDepth) {
        self.depth = Some(depth);
    }

    /// Set the distances to children to the input value. Implies this is an internal node.
    pub fn set_desc_distances(&mut self, distances: &Option<HashMap<NodeId, Edge>>) {
        self.desc_distances = distances.clone();
    }

    /// Set the distances to children to be zero. Implies that this Node is a leaf node.
    pub fn set_desc_distances_as_leaf(&mut self) {
        let mut map = HashMap::new();
        map.insert(self.id, Edge(0.0));
        self.set_desc_distances(&Some(map));
    }

    /// Check if this is a leaf node (i.e. no children).
    #[must_use]
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// Check if this is a root node (i.e. no parents).
    #[must_use]
    pub fn is_root(&self) -> bool {
        self.parent.is_none()
    }
}
