use std::collections::HashMap;

use itertools::Itertools;
use light_phylogeny::ArenaTree as LpTree;
use light_phylogeny::read_newick;
use numpy::ndarray::Array2;

use crate::node::Node;
use crate::types::{Edge, NodeDepth, NodeId, Taxon};
use crate::util::{create_row_vec_from_mat_dims, row_idx_from_mat_coords, row_vec_to_symmat};

#[derive(Default)]
pub struct ArenaTree
{
    pub nodes: Vec<Node>,
    pub taxon_to_node_id: HashMap<Taxon, NodeId>,
    pub leaf_idx_to_row_idx: HashMap<NodeId, usize>,
    pub row_idx_to_leaf_idx: Vec<NodeId>,
    pub nodes_at_depth: HashMap<NodeDepth, Vec<NodeId>>,
}

impl ArenaTree
{
    /// Return the number of nodes (leaf + internal) in the tree.
    #[must_use]
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn leaf_nodes(&self) -> Vec<Taxon> {
        let mut out = Vec::new();
        for leaf_idx in self.row_idx_to_leaf_idx.iter() {
            let leaf = self.get_node(*leaf_idx);
            out.push(leaf.taxon.clone().unwrap());
        }
        out
    }

    /// Return the sum of all branches in the tree.
    #[must_use]
    pub fn length(&self) -> Edge {
        self.nodes.iter().map(|n| n.parent_distance.unwrap_or(Edge(0.0))).sum()
    }

    /// Returns a vector of pointers to nodes at a specific depth. Panics if depth doesn't exist.
    fn get_node_idxs_at_depth(&self, depth: NodeDepth) -> Vec<NodeId> {
        let nodes_at_depth = match self.nodes_at_depth.get(&depth) {
            Some(node_ids) => node_ids,
            None => panic!("No nodes at depth: {:?}", depth),
        };
        let mut nodes: Vec<NodeId> = vec![];
        for node_id in nodes_at_depth.iter() {
            nodes.push(*node_id);
        }
        nodes
    }

    /// Return the number of leaf nodes in the tree.
    #[must_use]
    pub fn n_leaf_nodes(&self) -> usize {
        self.leaf_idx_to_row_idx.len()
    }

    /// Returns the row vector index for a given leaf node id.
    fn get_row_vec_idx_from_leaf_idx(&self, leaf_id: NodeId) -> usize {
        *self.leaf_idx_to_row_idx.get(&leaf_id).unwrap()
    }

    /// Get the row vector index for two taxa in the tree.
    fn get_row_vec_idx_dist_between_leaf_idx(&self, i: NodeId, j: NodeId) -> usize {
        let n_taxa = self.n_leaf_nodes();
        let i_idx = self.get_row_vec_idx_from_leaf_idx(i);
        let j_idx = self.get_row_vec_idx_from_leaf_idx(j);
        row_idx_from_mat_coords(n_taxa, i_idx, j_idx)
    }

    /// Add a new leaf node to the tree.
    /// Returns the ID of the new node.
    pub fn add_leaf_node(&mut self, taxon: Taxon) -> NodeId {

        // Panic if the taxon is already in the tree.
        if self.taxon_to_node_id.contains_key(&taxon) {
            panic!("Taxon already exists in the tree: '{:?}'", taxon);
        }

        // Create the new node, and place it in the tree.
        let node_id = NodeId(self.n_nodes());
        self.taxon_to_node_id.insert(taxon.clone(), node_id);
        self.leaf_idx_to_row_idx.insert(node_id, self.leaf_idx_to_row_idx.len());
        self.row_idx_to_leaf_idx.push(node_id);
        self.nodes.push(Node::new(node_id, Some(taxon.clone())));
        return self.nodes.last().unwrap().id;
    }

    /// Add a new internal node to the tree.
    pub fn add_internal_node(&mut self) -> NodeId {
        let node_id = NodeId(self.n_nodes());
        self.nodes.push(Node::new(node_id, None));
        return self.nodes.last().unwrap().id;
    }

    /// Add a node to the tree.
    // TODO: Remove this method.
    pub fn add_node(&mut self, taxon: Option<Taxon>) -> NodeId {
        match taxon {
            Some(t) => self.add_leaf_node(t),
            None => self.add_internal_node(),
        }
    }

    /// Retrieve a node from the tree.
    #[must_use]
    pub fn get_node(&self, node_id: NodeId) -> &Node {
        &self.nodes[node_id.0]
    }

    /// Retrieve a mutable node from the tree.
    pub fn get_node_mut(&mut self, node_id: NodeId) -> &mut Node {
        &mut self.nodes[node_id.0]
    }

    /// Add an edge to the tree.
    pub fn add_edge(&mut self, parent: NodeId, child: NodeId, length: Edge) {
        self.get_node_mut(parent).add_child(child);
        self.get_node_mut(child).set_parent(parent, length);
    }

    /// Set the depth of each node in the tree.
    pub fn assign_node_depth(&mut self) {

        // Iterate over each node to make sure there is only one root node.
        let mut root = None;
        self.nodes.iter().for_each(|node| {
            if node.is_root() {
                if root.is_some() {
                    panic!("Multiple root nodes detected!");
                }
                root = Some(node.id);
            }
        });

        // Panic if no root node could be found.
        if root.is_none() {
            panic!("No root node found!");
        }


        // Set the depth of all nodes
        self.set_node_depth_dfs(root.unwrap());
    }

    fn set_node_depth_dfs(&mut self, root_id: NodeId) {

        let mut nodes_at_depth: HashMap<NodeDepth, Vec<NodeId>> = HashMap::new();

        // Seed the stack with the root node.
        let mut stack = vec![(root_id, NodeDepth(0))];

        // Iterate over the stack.
        while let Some((node_id, depth)) = stack.pop() {

            // Load the node
            let node = self.get_node_mut(node_id);

            // Add the node to the hashmap.
            node.set_depth(depth);
            nodes_at_depth.entry(depth).or_insert_with(Vec::new).push(node.id);

            // Add the children to the stack
            node.children.iter().copied().for_each(|child_id| {
                stack.push((child_id, depth + NodeDepth(1)));
            });
        }

        // Check that the root is the only node at depth 0.
        let nodes_at_depth_0 = match { nodes_at_depth.get(&NodeDepth(0)) } {
            Some(nodes) => nodes,
            None => panic!("Root node not found!"),
        };
        assert!(!(nodes_at_depth_0.len() != 1 && nodes_at_depth_0[0] != root_id), "Root node not found!");

        // Check that the deepest nodes have no children.
        let deepest_node_depth = NodeDepth(nodes_at_depth.len() - 1);
        nodes_at_depth.get(&deepest_node_depth).into_iter().for_each(|nodes| {
            for node_id in nodes.iter() {
                let node = self.get_node(*node_id);
                assert!(node.is_leaf(), "Node has children: {:?}", node_id);
            }
        });

        // Save the hashmap
        self.nodes_at_depth = nodes_at_depth;
        // self.nodes_at_depth = nodes_at_depth;
    }

    /// Wrapper method to calculate the pairwise distances at a given depth.
    fn calculate_distances_at_depth(&mut self, depth: NodeDepth, row_vec: &mut [f64]) {

        // Iterate over all nodes a this depth
        self.get_node_idxs_at_depth(depth).into_iter().for_each(|node_id| {
            let node = self.get_node(node_id);

            // This is a leaf node, set the child distances to 0
            if node.is_leaf() {
                self.get_node_mut(node_id).set_desc_distances_as_leaf();
            } else if node.children.len() == 1 {
                // This only happens with malformed trees, bring forward the distances.
                let child_node = self.get_node(node.children[0]);
                if child_node.desc_distances.is_some() {
                    let mut new_desc_distances = child_node.desc_distances.as_ref().unwrap().clone();
                    for (_k, v) in new_desc_distances.iter_mut() {
                        *v = *v + node.parent_distance.unwrap_or(Edge(0.0));
                    }
                    self.get_node_mut(node_id).set_desc_distances(&Some(new_desc_distances));
                } else {
                    panic!("Unknown error, please report this.")
                }
            } else {
                // 1. Set the descendant distances for this node.
                self.set_node_descendant_distance(node_id);

                // 2. Calculate the pairwise distances for the leaf nodes.
                self.calc_pairwise_distances_to_leaf_nodes(node_id, row_vec);

                // Free un-used memory
                self.unset_node_child_distances(node_id);
            }
        });

    }

    /// For a given node, iterate over its children and unset the descendant distances.
    /// This is mainly to free memory during distance matrix calculation.
    fn unset_node_child_distances(&mut self, node_id: NodeId) {
        self.get_node(node_id).children.clone().iter().for_each(|child_id| {
            self.get_node_mut(*child_id).desc_distances = None;
        });
    }


    /// This function sets the descendant distances for a given node.
    /// Note: This is only the distances to all descendant nodes (including leaf nodes).
    /// The pairwise leaf distances are calculated in another method.
    fn set_node_descendant_distance(&mut self, node_id: NodeId) {
        let mut node_desc_distances: HashMap<NodeId, Edge> = HashMap::new();

        let node = self.get_node(node_id);

        // Iterate over each child node and bring forward the distances.
        node.children.clone().into_iter().for_each(|child_idx| {

            // Load the child node, and the descendant distances
            // Expects the descendant distances were already initialised.
            let child_node = self.get_node(child_idx);
            let child_dist = child_node.parent_distance.unwrap();
            let child_dists = child_node.desc_distances.as_ref().unwrap();

            // Bring forward the descendant distances to the children.
            for (grandchild_idx, grandchild_dist) in child_dists {
                node_desc_distances.insert(*grandchild_idx, child_dist + *grandchild_dist);
            }
        });
        self.get_node_mut(node_id).desc_distances = Some(node_desc_distances);
    }



    /// This function calculates the pairwise distances to all leaf nodes.
    /// Assumes that the memory has not been freed for the mapping.
    fn calc_pairwise_distances_to_leaf_nodes(&self, node_id: NodeId, row_vec: &mut [f64]) {

        let node = self.get_node(node_id);
        let children_idxs = node.children.clone();

        // Iterate over each child for this node, and do a pairwise calculation for each of the leaf
        // nodes in each of the children, defined by the desc_distances.
        for i in 0..children_idxs.len() {
            let child_i_idx = children_idxs[i];
            let child_i_node = self.get_node(child_i_idx);
            let child_i_parent_distance = child_i_node.parent_distance.unwrap();
            let child_i_desc_distances = child_i_node.desc_distances.as_ref().unwrap();

            for j in 0..i {
                let child_j_idx = children_idxs[j];
                let child_j_node = self.get_node(child_j_idx);
                let child_j_parent_distance = child_j_node.parent_distance.unwrap();
                let child_j_desc_distances = child_j_node.desc_distances.as_ref().unwrap();

                // Calculate the pairwise distances between these nodes.
                for (child_i_node_id, child_i_dist) in child_i_desc_distances.iter() {
                    for (child_j_node_id, child_j_dist) in child_j_desc_distances.iter() {

                        // Find the corresponding row vector index to store this comparison
                        let row_idx = self.get_row_vec_idx_dist_between_leaf_idx(*child_i_node_id, *child_j_node_id);

                        // Calculate the new distance between these nodes
                        let new_dist =  *child_i_dist + *child_j_dist  + child_i_parent_distance + child_j_parent_distance;

                        // Save this distance in the row vector
                        row_vec[row_idx] = new_dist.0;
                    }
                }
            }
        }


    }

    /// Orders the leaf nodes for reproducibility.
    fn order_leaf_node_idx(&mut self) {
        let mut new_leaf_idx_to_row_idx: HashMap<NodeId, usize> = HashMap::new();
        let mut new_row_idx_to_leaf_idx: Vec<NodeId> = vec![NodeId::default(); self.n_leaf_nodes()];
        for (new_idx, (_taxon, node_id)) in self.taxon_to_node_id.iter().sorted_by_key(|x| x.0).enumerate() {
            new_leaf_idx_to_row_idx.insert(*node_id, new_idx);
            new_row_idx_to_leaf_idx[new_idx] = *node_id;
        }
        self.leaf_idx_to_row_idx = new_leaf_idx_to_row_idx;
        self.row_idx_to_leaf_idx = new_row_idx_to_leaf_idx;
    }

    pub fn dm(&mut self, norm: bool) -> Array2<f64> {
        // For reproducibility, order the taxa
        self.order_leaf_node_idx();

        // Create the output row matrix
        let num_leaf = self.n_leaf_nodes();
        let mut row_vec = create_row_vec_from_mat_dims(num_leaf);

        // Compute the depth of each node
        // TODO: No need to do this again if no new nodes have been added.
        self.assign_node_depth();

        // Process the deepest nodes first
        let depths = self.nodes_at_depth.keys().sorted().rev().copied().collect::<Vec<_>>();
        for cur_depth in depths {
            self.calculate_distances_at_depth(cur_depth, &mut row_vec);
        }

        let mut array = row_vec_to_symmat(&row_vec);

        if norm {
            let tree_length = self.length();
            array.mapv_inplace(|x| x / tree_length.0);
        }
        array
    }


    pub fn load_from_newick_path(&mut self, path: &str) {
        let mut lp_tree: LpTree<String> = LpTree::default();
        read_newick(path.to_string(), &mut lp_tree);
        for cur_node in &lp_tree.arena {
            if cur_node.name.is_empty() {
                self.add_node(None);
            } else {
                self.add_node(Some(Taxon(cur_node.name.clone())));
            }
        }
        for cur_node in &lp_tree.arena {
            if cur_node.parent.is_some() {
                let parent_node = cur_node.parent.unwrap();
                self.add_edge(NodeId(parent_node), NodeId(cur_node.idx), Edge(cur_node.l as f64));
            }
        }
    }
}

