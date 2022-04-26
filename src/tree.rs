use std::borrow::BorrowMut;
use std::collections::HashMap;
use itertools::Itertools;
use light_phylogeny::ArenaTree as LpTree;
use light_phylogeny::read_newick;
use numpy::ndarray::Array2;
use crate::node::Node;


use crate::types::{BranchLength, NodeDepth, NodeId};
use crate::util::{row_idx_from_mat_coords, row_vec_to_symmat};

#[derive(Default)]
pub struct ArenaTree
{
    pub nodes: Vec<Node>,
    pub taxon_to_node_id: HashMap<String, NodeId>,
    pub leaf_idx_to_row_idx: HashMap<NodeId, usize>,
    pub row_idx_to_leaf_idx: Vec<NodeId>,
    pub nodes_at_depth: HashMap<NodeDepth, Vec<NodeId>>,
}

impl ArenaTree
{
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn length(&self) -> BranchLength {
        self.nodes.iter().map(|n| n.parent_distance.unwrap_or(BranchLength(0.0))).sum()
    }

    pub fn n_leaf_nodes(&self) -> usize {
        return self.leaf_idx_to_row_idx.len();
    }

    pub fn add_leaf_node(&mut self, taxon: &String) -> &NodeId {
        // Check if the taxon already exists
        if self.taxon_to_node_id.contains_key(taxon) {
            panic!("Taxon {} already exists", taxon);
        }

        // Add the taxon to the tree
        let node_id = NodeId(self.n_nodes());
        self.taxon_to_node_id.insert(taxon.clone(), node_id);
        self.leaf_idx_to_row_idx.insert(node_id, self.leaf_idx_to_row_idx.len());
        self.row_idx_to_leaf_idx.push(node_id);

        // Create the new node
        self.nodes.push(Node::new(node_id, Some(taxon.clone())));
        return &self.nodes.last().unwrap().id;

    }

    pub fn add_internal_node(&mut self) -> &NodeId {
        let node_id = NodeId(self.n_nodes());
        self.nodes.push(Node::new(node_id, None));
        return &self.nodes.last().unwrap().id;
    }

    pub fn add_node(&mut self, taxon: &Option<String>) -> &NodeId {
        return match taxon {
            Some(t) => self.add_leaf_node(t),
            None => self.add_internal_node(),
        }
    }

    pub fn get_node(&self, node_id: &NodeId) -> &Node {
        return &self.nodes[node_id.0];
    }

    pub fn get_node_mut(&mut self, node_id: &NodeId) -> &mut Node {
        return &mut self.nodes[node_id.0];
    }

    pub fn add_edge(&mut self, parent: &NodeId, child: &NodeId, length: &BranchLength) {
        self.get_node_mut(parent).add_child(child);
        self.get_node_mut(child).set_parent(parent, length);
    }

    pub fn assign_node_depth(&mut self) {
        /*
        Set the depth of each node in the tree.
         */
        let mut root = None;
        for node in self.nodes.iter() {
            if node.parent.is_none() {
                if root.is_some() {
                    panic!("Multiple root nodes detected");
                }
                root = Some(node.id);
            }
        }
        if root.is_none() {
            panic!("No root node found");
        }

        // Set the depth of all nodes
        self.set_node_depth_rec(&root.unwrap(), &NodeDepth(0));
    }

    fn set_node_depth_rec(&mut self, node_id: &NodeId, depth: &NodeDepth) {
        let cur_node = self.get_node_mut(node_id);
        cur_node.set_depth(depth);

        for child in cur_node.children.clone().iter() {
            self.set_node_depth_rec(child, &(NodeDepth(1) + *depth));
        }

        // save in arena
        if !self.nodes_at_depth.contains_key(depth) {
            self.nodes_at_depth.insert(*depth, vec![]);
        }
        self.nodes_at_depth.get_mut(depth).unwrap().push(*node_id);
    }


    pub fn process_nodes_at_depth(&mut self, depth: &NodeDepth, row_vec: &mut Vec<f64>) {
        let nodes_at_depth = self.nodes_at_depth.get(depth).unwrap().clone();
        self.mutate_row_vec(row_vec, &nodes_at_depth);
    }


    pub fn mutate_row_vec(&mut self,
                          row_vec: &mut Vec<f64>,
                          nodes_at_depth: &Vec<NodeId>) {
        for node_id in nodes_at_depth {
            let node = self.get_node(node_id);

            if node.children.len() == 0 {
                // Set leaf node descendant distances to be 0 (to itself).
                self.get_node_mut(&node_id).set_desc_distances_as_leaf();
                continue;
            } else if node.children.len() == 1 {
                // This only happens with malformed trees, bring forward the distances.
                let child_node = self.get_node(&node.children[0]);
                if child_node.desc_distances.is_some() {
                    let mut new_desc_distances = child_node.desc_distances.as_ref().unwrap().clone();
                    for (_k, v) in new_desc_distances.iter_mut() {
                        *v = *v + node.parent_distance.unwrap_or(BranchLength(0.0));
                    }
                    self.get_node_mut(&node_id).set_desc_distances(&Some(new_desc_distances));
                    continue;
                } else {
                    panic!("??")
                }
            } else {
                // Calculate pairwise distances to each taxon.
                let child_idxs: Vec<&NodeId> = node.children.iter().map(|child_id| child_id).collect();
                self.get_node_mut(&node_id).desc_distances = Some(self.calc_pairwise_node_data(&child_idxs, row_vec));

                // free memory
                for child_idx in self.get_node(node_id).children.clone().iter() {
                    self.get_node_mut(child_idx).desc_distances = None;
                }
            }
        }
    }

    pub fn pairwise_node_data_worker(&self, i: usize, child_idxs: &Vec<&NodeId>) -> (Vec<(NodeId, BranchLength)>, Vec<(usize, BranchLength)>) {
        let mut node_desc_distances: Vec<(NodeId, BranchLength)> = Vec::new();
        let mut row_vec_data: Vec<(usize, BranchLength)> = Vec::new();

        let n_taxa = self.n_leaf_nodes();

        let child_i = self.get_node(child_idxs[i]);
        for (k, v) in child_i.desc_distances.as_ref().unwrap().iter() {
            node_desc_distances.push((*k, *v + child_i.parent_distance.unwrap_or(BranchLength(0.0))));
        }

        for j in 0..i {
            let child_j = self.get_node(child_idxs[j]);
            for (k, v) in child_j.desc_distances.as_ref().unwrap().iter() {
                node_desc_distances.push((*k, *v + child_j.parent_distance.unwrap_or(BranchLength(0.0))));
            }

            // iterate over all of the taxa nodes and add the branch lengths
            for (child_i_taxon, child_i_dist) in child_i.desc_distances.as_ref().unwrap().iter() {
                let child_i_taxon_row_idx = self.leaf_idx_to_row_idx[child_i_taxon];
                for (child_j_taxon, child_j_dist) in child_j.desc_distances.as_ref().unwrap().iter() {
                    let child_j_taxon_row_idx = self.leaf_idx_to_row_idx[child_j_taxon];
                    let row_idx = row_idx_from_mat_coords(n_taxa as i64, child_i_taxon_row_idx as i64, child_j_taxon_row_idx as i64);
                    let new_dist = *child_i_dist + *child_j_dist + child_i.parent_distance.unwrap_or(BranchLength(0.0))
                        + child_j.parent_distance.unwrap_or(BranchLength(0.0));
                    row_vec_data.push((row_idx, new_dist));
                }
            }
        }
        return (node_desc_distances, row_vec_data);
    }

    pub fn calc_pairwise_node_data(&self, child_idxs: &Vec<&NodeId>, row_vec: &mut Vec<f64>) -> HashMap<NodeId, BranchLength> {
        let mut node_desc_distances: HashMap<NodeId, BranchLength> = HashMap::new();

        let res = (0..child_idxs.len()).into_iter()
            .map(|child_idx| {
                return self.pairwise_node_data_worker(child_idx, child_idxs);
            })
            .collect::<Vec<_>>();

        for cur_res in &res {
            for (k, v) in cur_res.0.iter() {
                node_desc_distances.insert(*k, *v);
            }
            for (k, v) in cur_res.1.iter() {
                row_vec[*k] = v.0;
            }
        }
        return node_desc_distances;
    }

    pub fn order_leaf_node_idx(&mut self) {
        let mut new_leaf_idx_to_row_idx: HashMap<NodeId, usize> = HashMap::new();
        let mut new_row_idx_to_leaf_idx: Vec<NodeId> = Vec::new();
        for (_taxon, node_id) in self.taxon_to_node_id.iter().sorted_by_key(|x| x.0) {
            new_leaf_idx_to_row_idx.insert(*node_id, new_row_idx_to_leaf_idx.len());
            new_row_idx_to_leaf_idx.push(*node_id);
        }
        self.leaf_idx_to_row_idx = new_leaf_idx_to_row_idx;
        self.row_idx_to_leaf_idx = new_row_idx_to_leaf_idx;
    }

    pub fn dm(&mut self, norm: &bool) -> Array2<f64> {
        self.order_leaf_node_idx();

        // Create the output row matrix
        let num_leaf = self.n_leaf_nodes();
        let mut row_vec = vec![0.0; (&num_leaf * (&num_leaf + 1)) / 2];

        // Get the depth of each node
        self.assign_node_depth();

        let depths = self.nodes_at_depth.keys().sorted_by_key(|x| x.0).rev().cloned().collect::<Vec<_>>();

        for cur_depth in depths.iter() {
            self.process_nodes_at_depth(cur_depth, &mut row_vec);
        }

        let mut array = Array2::<f64>::zeros([num_leaf, num_leaf]);

        row_vec_to_symmat(&row_vec, array.borrow_mut());

        if *norm {
            let tree_length = self.length();
            array.mapv_inplace(|x| x / tree_length.0);
        }
        return array;
    }

    pub fn load_from_newick_path(&mut self, path: &str) {
        let mut lp_tree: LpTree<String> = LpTree::default();
        read_newick(path.to_string(), &mut lp_tree);
        for cur_node in &lp_tree.arena {
            if cur_node.name.len() > 0 {
                self.add_node(&Some(cur_node.name.clone()));
            } else {
                self.add_node(&None);
            }
        }
        for cur_node in &lp_tree.arena {
            if cur_node.parent.is_some() {
                let parent_node = cur_node.parent.unwrap();
                self.add_edge(&NodeId(parent_node), &NodeId(cur_node.idx), &BranchLength(cur_node.l as f64));
            }
        }
    }

}

