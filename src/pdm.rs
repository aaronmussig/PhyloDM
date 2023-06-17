use std::collections::HashMap;
use std::panic;

use itertools::Itertools;
use light_phylogeny::ArenaTree as LpTree;
use light_phylogeny::read_newick;
use ndarray::Array2;

use crate::error::PhyloErr;
use crate::tree::{Edge, NodeDepth, NodeId, Taxon};
use crate::tree::Node;
use crate::util::{create_row_vec_from_mat_dims, row_idx_from_mat_coords, row_vec_to_symmat};

/// Create and manipulate the Phylogenetic Distance Matrix.
///
/// # Examples
///
/// ```
/// use phylodm::PDM;
///
/// let mut tree = PDM::default();
/// tree.load_from_newick_path("examples/example.tree");
///
/// let (taxa, dm) = tree.matrix(true).unwrap();
/// ```
#[derive(Default)]
pub struct PDM {
    pub(crate) nodes: Vec<Node>,
    pub(crate) taxon_to_node_id: HashMap<Taxon, NodeId>,
    pub(crate) leaf_idx_to_row_idx: HashMap<NodeId, usize>,
    pub(crate) leaf_idx_to_row_idx_vec: Vec<usize>,
    pub(crate) row_idx_to_leaf_idx: Vec<NodeId>,
    pub(crate) nodes_at_depth: HashMap<NodeDepth, Vec<NodeId>>,
    pub(crate) row_vec: Option<Vec<f64>>,
}

impl PDM {
    /// Return the number of nodes (leaf + internal) in the tree.
    #[must_use]
    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Return all leaf nodes in the tree.
    ///
    /// # Errors
    /// In the event that a leaf node has no taxon, this function will return an error.
    /// This case should never happen, but is possible due to the way the tree is constructed.
    pub fn leaf_nodes(&self) -> Result<Vec<Taxon>, PhyloErr> {
        let mut out = Vec::with_capacity(self.row_idx_to_leaf_idx.len());
        for leaf_idx in &self.row_idx_to_leaf_idx {
            let leaf = self.get_node(*leaf_idx);
            let Some(taxon) = &leaf.taxon else {
                return Err(PhyloErr("Leaf node has no taxon! Please report this error.".to_string()));
            };
            out.push(taxon.clone());
        }
        Ok(out)
    }

    /// Return the sum of all branches in the tree.
    #[must_use]
    pub fn length(&self) -> Edge {
        self.nodes
            .iter()
            .map(|n| n.parent_distance.unwrap_or(Edge(0.0)))
            .sum()
    }

    /// Returns a vector of node indices at a specific depth. Errors if depth doesn't exist.
    fn get_node_idxs_at_depth(&self, depth: NodeDepth) -> Result<Vec<NodeId>, PhyloErr> {
        if self.nodes_at_depth.get(&depth).is_none() {
            return Err(PhyloErr("No nodes at depth! Please report this error.".to_string()));
        }
        let Some(nodes_at_depth) = self.nodes_at_depth.get(&depth) else {
            return Err(PhyloErr("No nodes at depth! Please report this error.".to_string()));
        };
        let mut nodes: Vec<NodeId> = Vec::with_capacity(nodes_at_depth.len());
        for node_id in nodes_at_depth {
            nodes.push(*node_id);
        }
        Ok(nodes)
    }

    /// Return the number of leaf nodes in the tree.
    #[must_use]
    pub fn n_leaf_nodes(&self) -> usize {
        self.leaf_idx_to_row_idx.len()
    }

    /// Returns the row vector index for a given leaf node id.
    #[must_use]
    fn get_row_vec_idx_from_leaf_idx(&self, leaf_id: NodeId) -> usize {
        self.leaf_idx_to_row_idx_vec[leaf_id.0]
    }

    /// Get the row vector index for two taxa in the tree.
    #[must_use]
    fn get_row_vec_idx_dist_between_leaf_idx(&self, i: NodeId, j: NodeId) -> usize {
        let n_taxa = self.n_leaf_nodes();
        let i_idx = self.get_row_vec_idx_from_leaf_idx(i);
        let j_idx = self.get_row_vec_idx_from_leaf_idx(j);
        row_idx_from_mat_coords(n_taxa, i_idx, j_idx)
    }

    /// Add a new leaf node to the tree.
    /// Returns the ID of the new node.
    fn add_leaf_node(&mut self, taxon: &Taxon) -> Result<NodeId, PhyloErr> {
        // Panic if the taxon is already in the tree.
        if self.taxon_to_node_id.contains_key(taxon) {
            return Err(PhyloErr(format!("Taxon already exists in the tree: '{taxon:?}'")));
        }

        // Create the new node, and place it in the tree.
        let node_id = NodeId(self.n_nodes());
        self.taxon_to_node_id.insert(taxon.clone(), node_id);
        self.leaf_idx_to_row_idx.insert(node_id, self.leaf_idx_to_row_idx.len());
        self.row_idx_to_leaf_idx.push(node_id);
        self.nodes.push(Node::new(node_id, Some(taxon.clone())));
        return Ok(self.nodes.last().unwrap().id);
    }

    /// Add a new internal node to the tree.
    fn add_internal_node(&mut self) -> NodeId {
        let node_id = NodeId(self.n_nodes());
        self.nodes.push(Node::new(node_id, None));
        return self.nodes.last().unwrap().id;
    }

    /// Add a node to the tree.
    pub(crate) fn add_node(&mut self, taxon: Option<&Taxon>) -> Result<NodeId, PhyloErr> {
        match taxon {
            Some(t) => self.add_leaf_node(t),
            None => Ok(self.add_internal_node()),
        }
    }

    /// Retrieve a node from the tree.
    #[must_use]
    fn get_node(&self, node_id: NodeId) -> &Node {
        &self.nodes[node_id.0]
    }

    /// Retrieve a mutable node from the tree.
    #[must_use]
    fn get_node_mut(&mut self, node_id: NodeId) -> &mut Node {
        &mut self.nodes[node_id.0]
    }

    /// Add an edge to the tree.
    ///
    /// # Arguments
    ///
    /// * `parent`: - `NodeId` of the parent node.
    /// * `child`:  - `NodeId` of the child node.
    /// * `length`: - The branch length between these nodes.
    ///
    pub(crate) fn add_edge(&mut self, parent: NodeId, child: NodeId, length: Edge) {
        self.get_node_mut(parent).add_child(child);
        self.get_node_mut(child).set_parent(parent, length);
    }

    /// Set the depth of each node in the tree.
    fn assign_node_depth(&mut self) -> Result<(), PhyloErr> {
        // Iterate over each node to make sure there is only one root node.
        let mut root = None;
        for node in &self.nodes {
            if node.is_root() {
                if root.is_some() {
                    return Err(PhyloErr("Multiple root nodes detected!".to_string()));
                }
                root = Some(node.id);
            }
        }

        if root.is_none() {
            return Err(PhyloErr("No root node detected!".to_string()));
        }

        // Set the depth of all nodes
        self.set_node_depth_dfs(root.unwrap())?;
        Ok(())
    }

    /// Set the depth of all nodes using a depth first search.
    fn set_node_depth_dfs(&mut self, root_id: NodeId) -> Result<(), PhyloErr> {
        let mut nodes_at_depth: HashMap<NodeDepth, Vec<NodeId>> = HashMap::new();

        // Seed the stack with the root node.
        let mut stack = vec![(root_id, NodeDepth(0))];

        // Iterate over the stack.
        while let Some((node_id, depth)) = stack.pop() {
            // Load the node
            let node = self.get_node_mut(node_id);

            // Add the node to the hashmap.
            node.set_depth(depth);
            nodes_at_depth
                .entry(depth)
                .or_insert_with(Vec::new)
                .push(node.id);

            // Add the children to the stack
            node.children.iter().for_each(|child_id| {
                stack.push((*child_id, depth + NodeDepth(1)));
            });
        }

        // Check that the root is the only node at depth 0.
        if nodes_at_depth.get(&NodeDepth(0)).is_none() {
            return Err(PhyloErr("Root node not found!".to_string()));
        }
        let Some(nodes_at_depth_0) = nodes_at_depth.get(&NodeDepth(0)) else {
            return Err(PhyloErr("No nodes were found at depth 0, report this error.".to_string()));
        };
        if nodes_at_depth_0.len() != 1 && nodes_at_depth_0[0] != root_id {
            return Err(PhyloErr("Root node not found!".to_string()));
        }

        // Check that the deepest nodes have no children.
        let deepest_node_depth = NodeDepth(nodes_at_depth.len() - 1);
        let Some(deepest_nodes) = nodes_at_depth.get(&deepest_node_depth) else {
            return Err(PhyloErr("No nodes were found at depth max, report this error.".to_string()));
        };
        for node_id in deepest_nodes {
            let node = self.get_node(*node_id);
            if !node.is_leaf() {
                return Err(PhyloErr("Node has children!".to_string()));
            }
        }

        // Save the hashmap
        self.nodes_at_depth = nodes_at_depth;
        Ok(())
    }

    /// Wrapper method to calculate the pairwise distances at a given depth.
    fn calculate_distances_at_depth(&mut self, depth: NodeDepth, row_vec: &mut [f64]) -> Result<(), PhyloErr> {
        // Iterate over all nodes a this depth
        for &node_id in &self.get_node_idxs_at_depth(depth)? {
            let node = self.get_node(node_id);

            // This is a leaf node, set the child distances to 0
            if node.is_leaf() {
                self.get_node_mut(node_id).set_desc_distances_as_leaf();
            } else if node.children.len() == 1 {
                // This only happens with malformed trees, bring forward the distances.
                let child_node = self.get_node(node.children[0]);
                if child_node.desc_distances.is_some() {
                    let mut new_desc_distances =
                        child_node.desc_distances.as_ref().unwrap().clone();
                    for edge in new_desc_distances.values_mut() {
                        *edge = *edge + node.parent_distance.unwrap_or(Edge(0.0));
                    }
                    self.get_node_mut(node_id)
                        .set_desc_distances(&Some(new_desc_distances));
                } else {
                    return Err(PhyloErr("Unknown error, please report this.".to_string()));
                }
            } else {
                // 1. Set the descendant distances for this node.
                self.set_node_descendant_distance(node_id);

                // 2. Calculate the pairwise distances for the leaf nodes.
                self.calc_pairwise_distances_to_leaf_nodes(node_id, row_vec);

                // Free un-used memory
                self.unset_node_child_distances(node_id);
            }
        }
        Ok(())
    }

    /// For a given node, iterate over its children and unset the descendant distances.
    /// This is mainly to free memory during distance matrix calculation.
    fn unset_node_child_distances(&mut self, node_id: NodeId) {
        self.get_node(node_id)
            .children
            .clone()
            .into_iter()
            .for_each(|child_id| {
                self.get_node_mut(child_id).desc_distances = None;
            });
    }

    /// This function sets the descendant distances for a given node.
    /// Note: This is only the distances to all descendant nodes (including leaf nodes).
    /// The pairwise leaf distances are calculated in another method.
    fn set_node_descendant_distance(&mut self, node_id: NodeId) {
        let mut node_desc_distances: HashMap<NodeId, Edge> = HashMap::new();

        let node = self.get_node(node_id);

        // Iterate over each child node and bring forward the distances.
        for child_idx in &node.children {
            // Load the child node, and the descendant distances
            // Expects the descendant distances were already initialised.
            let child_node = self.get_node(*child_idx);
            let child_dist = child_node.parent_distance.unwrap();
            let child_dists = child_node.desc_distances.as_ref().unwrap();

            // Bring forward the descendant distances to the children.
            for (grandchild_idx, grandchild_dist) in child_dists {
                node_desc_distances.insert(*grandchild_idx, child_dist + *grandchild_dist);
            }
        }

        self.get_node_mut(node_id).desc_distances = Some(node_desc_distances);
    }

    /// This function calculates the pairwise distances to all leaf nodes.
    /// Assumes that the memory has not been freed for the mapping.
    fn calc_pairwise_distances_to_leaf_nodes(&self, node_id: NodeId, row_vec: &mut [f64]) {
        let node = self.get_node(node_id);
        let children_idxs = &node.children;

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
                for (child_i_node_id, child_i_dist) in child_i_desc_distances {
                    for (child_j_node_id, child_j_dist) in child_j_desc_distances {
                        // Find the corresponding row vector index to store this comparison
                        let row_idx = self.get_row_vec_idx_dist_between_leaf_idx(
                            *child_i_node_id,
                            *child_j_node_id,
                        );

                        // Calculate the new distance between these nodes
                        let new_dist = *child_i_dist
                            + *child_j_dist
                            + child_i_parent_distance
                            + child_j_parent_distance;

                        // Save this distance in the row vector
                        row_vec[row_idx] = new_dist.0;
                    }
                }
            }
        }
    }

    /// Orders the leaf nodes for reproducibility.
    fn order_leaf_node_idx(&mut self) {
        let mut new_leaf_idx_to_row_idx: HashMap<NodeId, usize> = HashMap::with_capacity(self.n_leaf_nodes());
        let mut new_row_idx_to_leaf_idx: Vec<NodeId> = vec![NodeId::default(); self.n_leaf_nodes()];

        // Find the maximum index for the leaf node
        let &max_leaf_idx = self.leaf_idx_to_row_idx.keys().max().unwrap();
        let mut new_row_idx_to_leaf_idx_vec: Vec<usize> = vec![0; max_leaf_idx.0 + 1];

        for (new_idx, (_taxon, node_id)) in self
            .taxon_to_node_id
            .iter()
            .sorted_by_key(|x| x.0)
            .enumerate()
        {
            new_leaf_idx_to_row_idx.insert(*node_id, new_idx);
            new_row_idx_to_leaf_idx[new_idx] = *node_id;
            new_row_idx_to_leaf_idx_vec[node_id.0] = new_idx;
        }
        self.leaf_idx_to_row_idx = new_leaf_idx_to_row_idx;
        self.row_idx_to_leaf_idx = new_row_idx_to_leaf_idx;
        self.leaf_idx_to_row_idx_vec = new_row_idx_to_leaf_idx_vec;
    }

    /// Return the symmetrical pairwise distance matrix.
    ///
    /// # Arguments
    /// * `norm` - True if the result should be normalized by the sum of all branches in the tree.
    ///
    /// # Errors
    /// If any errors are encountered due to unexpected tree structures, an error will be raised.
    pub fn matrix(&mut self, norm: bool) -> Result<(Vec<Taxon>, Array2<f64>), PhyloErr> {
        self.compute_row_vec()?;

        let mut array = row_vec_to_symmat(self.row_vec.as_ref().unwrap());

        if norm {
            let tree_length = self.length();
            array.mapv_inplace(|x| x / tree_length.0);
        }
        Ok((self.leaf_nodes()?, array))
    }

    /// Initialise the PDM from a newick file.
    ///
    /// # Errors
    /// If any errors are encountered due to unexpected tree structures, an error will be raised.
    pub fn load_from_newick_path(&mut self, path: &str) -> Result<(), PhyloErr> {
        // Catch any errors that light_phylogeny may throw
        // Suppress the stderr message
        let prev_hook = panic::take_hook();
        panic::set_hook(Box::new(|_| {}));
        let newick_ok = panic::catch_unwind(|| {
            let mut lp_tree: LpTree<String> = LpTree::default();
            read_newick(path.to_string(), &mut lp_tree);
            lp_tree
        });
        panic::set_hook(prev_hook);
        if newick_ok.is_err() {
            return Err(PhyloErr("Error reading newick file".to_string()));
        }
        let lp_tree = newick_ok.unwrap();

        // Import the tree
        for cur_node in &lp_tree.arena {
            if cur_node.name.is_empty() {
                self.add_node(None)?;
            } else {
                self.add_node(Some(&Taxon(cur_node.name.clone())))?;
            }
        }
        for cur_node in &lp_tree.arena {
            if cur_node.parent.is_some() {
                let parent_node = cur_node.parent.unwrap();
                self.add_edge(
                    NodeId(parent_node),
                    NodeId(cur_node.idx),
                    Edge(cur_node.l as f64),
                );
            }
        }
        self.compute_row_vec()?;
        Ok(())
    }

    /// Computes the row vector. Required if the PDM was manually created (i.e. not from a newick file).
    pub fn compute_row_vec(&mut self) -> Result<(), PhyloErr> {
        // For reproducibility, order the taxa
        self.order_leaf_node_idx();

        // Create the output row matrix
        let num_leaf = self.n_leaf_nodes();
        let mut row_vec = create_row_vec_from_mat_dims(num_leaf);

        // Compute the depth of each node
        // TODO: No need to do this again if no new nodes have been added.
        self.assign_node_depth()?;

        // Process the deepest nodes first
        let depths = self
            .nodes_at_depth
            .keys()
            .sorted()
            .rev()
            .copied()
            .collect::<Vec<_>>();
        for cur_depth in depths {
            self.calculate_distances_at_depth(cur_depth, &mut row_vec)?;
        }
        self.row_vec = Some(row_vec);
        Ok(())
    }

    fn get_taxon_node_idx(&self, taxon: &Taxon) -> NodeId {
        self.taxon_to_node_id[taxon]
    }

    /// Return the distance between two taxa.
    ///
    /// # Arguments
    /// * `a`: - The first taxon.
    /// * `b`: - The second taxon.
    /// * `norm` - True if the result should be normalised by the sum of all branches in the tree.
    #[must_use]
    pub fn distance(&self, a: &Taxon, b: &Taxon, norm: bool) -> f64 {
        assert!(
            !self.row_vec.is_none(),
            "The PDM has not been computed yet, call compute_row_vec() first."
        );
        let a_idx = self.get_taxon_node_idx(a);
        let b_idx = self.get_taxon_node_idx(b);
        let row_idx = self.get_row_vec_idx_dist_between_leaf_idx(a_idx, b_idx);
        let dist = self.row_vec.as_ref().unwrap()[row_idx];
        return if norm { dist / self.length().0 } else { dist };
    }
}

