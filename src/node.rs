use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::types::{BranchLength,  NodeDepth, NodeId};

#[derive(Debug)]
pub struct Node
{
    pub id: NodeId,
    pub taxon: Option<String>,
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
    pub parent_distance: Option<BranchLength>,
    pub depth: Option<NodeDepth>,
    pub desc_distances: Option<HashMap<NodeId, BranchLength>>,
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


impl Node
{
    pub fn new(id: NodeId, taxon: Option<String>) -> Self {
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

    pub fn add_child(&mut self, child: &NodeId) {
        self.children.push(*child);
    }

    pub fn set_parent(&mut self, parent: &NodeId, length: &BranchLength) {
        self.parent = Some(*parent);
        self.parent_distance = Some(*length);
    }

    pub fn set_depth(&mut self, depth: &NodeDepth) {
        self.depth = Some(*depth);
    }

    pub fn set_desc_distances(&mut self, distances: &Option<HashMap<NodeId, BranchLength>>) {
        self.desc_distances = distances.clone();
    }

    pub fn set_desc_distances_as_leaf(&mut self) {
        let mut map = HashMap::new();
        map.insert(self.id, BranchLength(0.0));
        self.set_desc_distances(&Some(map));
    }
}

