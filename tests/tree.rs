#[cfg(test)]
mod tests {
    use phylodm::PDM;
    use phylodm::tree::{Taxon, Edge};

    #[test]
    fn test_tree_dm_twice() {
        let mut tree = PDM::default();
        let _ = tree.load_from_newick_path("tests/test.tree");
        let arr = tree.matrix(false).unwrap();
        let arr2 = tree.matrix(false).unwrap();
        assert_eq!(arr, arr2);
    }

    #[test]
    fn test_tree() {
        let mut tree = PDM::default();
        let _ = tree.load_from_newick_path("tests/test.tree");
        let (_, arr) = tree.matrix(false).unwrap();

        assert_eq!(arr[[0, 0]], 0.0);
        assert_eq!(arr[[1, 0]], 84.0);
        assert_eq!(arr[[1, 1]], 0.0);
        assert_eq!(arr[[2, 0]], 72.0);
        assert_eq!(arr[[2, 1]], 28.0);
        assert_eq!(arr[[2, 2]], 0.0);
        assert_eq!(arr[[3, 0]], 31.0);
        assert_eq!(arr[[3, 1]], 83.0);
        assert_eq!(arr[[3, 2]], 71.0);
        assert_eq!(arr[[3, 3]], 0.0);
        assert_eq!(arr[[4, 0]], 85.0);
        assert_eq!(arr[[4, 1]], 17.0);
        assert_eq!(arr[[4, 2]], 29.0);
        assert_eq!(arr[[4, 3]], 84.0);
        assert_eq!(arr[[4, 4]], 0.0);
        assert_eq!(arr[[5, 0]], 72.0);
        assert_eq!(arr[[5, 1]], 48.0);
        assert_eq!(arr[[5, 2]], 36.0);
        assert_eq!(arr[[5, 3]], 71.0);
        assert_eq!(arr[[5, 4]], 49.0);
        assert_eq!(arr[[5, 5]], 0.0);
        assert_eq!(arr[[6, 0]], 55.0);
        assert_eq!(arr[[6, 1]], 53.0);
        assert_eq!(arr[[6, 2]], 41.0);
        assert_eq!(arr[[6, 3]], 54.0);
        assert_eq!(arr[[6, 4]], 54.0);
        assert_eq!(arr[[6, 5]], 41.0);
        assert_eq!(arr[[6, 6]], 0.0);
        assert_eq!(arr[[7, 0]], 47.0);
        assert_eq!(arr[[7, 1]], 71.0);
        assert_eq!(arr[[7, 2]], 59.0);
        assert_eq!(arr[[7, 3]], 46.0);
        assert_eq!(arr[[7, 4]], 72.0);
        assert_eq!(arr[[7, 5]], 59.0);
        assert_eq!(arr[[7, 6]], 42.0);
        assert_eq!(arr[[7, 7]], 0.0);
        assert_eq!(arr[[8, 0]], 71.0);
        assert_eq!(arr[[8, 1]], 27.0);
        assert_eq!(arr[[8, 2]], 5.0);
        assert_eq!(arr[[8, 3]], 70.0);
        assert_eq!(arr[[8, 4]], 28.0);
        assert_eq!(arr[[8, 5]], 35.0);
        assert_eq!(arr[[8, 6]], 40.0);
        assert_eq!(arr[[8, 7]], 58.0);
        assert_eq!(arr[[8, 8]], 0.0);
        assert_eq!(arr[[9, 0]], 75.0);
        assert_eq!(arr[[9, 1]], 21.0);
        assert_eq!(arr[[9, 2]], 19.0);
        assert_eq!(arr[[9, 3]], 74.0);
        assert_eq!(arr[[9, 4]], 22.0);
        assert_eq!(arr[[9, 5]], 39.0);
        assert_eq!(arr[[9, 6]], 44.0);
        assert_eq!(arr[[9, 7]], 62.0);
        assert_eq!(arr[[9, 8]], 18.0);
        assert_eq!(arr[[9, 9]], 0.0);

        assert_eq!(
            tree.leaf_nodes().unwrap(),
            [
                Taxon("T1".to_string()),
                Taxon("T10".to_string()),
                Taxon("T2".to_string()),
                Taxon("T3".to_string()),
                Taxon("T4".to_string()),
                Taxon("T5".to_string()),
                Taxon("T6".to_string()),
                Taxon("T7".to_string()),
                Taxon("T8".to_string()),
                Taxon("T9".to_string())
            ]
        );
    }

    #[test]
    fn test_tree_norm() {
        let mut tree = PDM::default();
        let _ = tree.load_from_newick_path("tests/test.tree");
        let (_, arr) = tree.matrix(true).unwrap();

        assert_eq!(arr[[0, 0]], 0.0);
        assert_eq!(arr[[1, 0]], 0.49122807017543857);
        assert_eq!(arr[[1, 1]], 0.0);
        assert_eq!(arr[[2, 0]], 0.42105263157894735);
        assert_eq!(arr[[2, 1]], 0.16374269005847952);
        assert_eq!(arr[[2, 2]], 0.0);
        assert_eq!(arr[[3, 0]], 0.18128654970760233);
        assert_eq!(arr[[3, 1]], 0.4853801169590643);
        assert_eq!(arr[[3, 2]], 0.4152046783625731);
        assert_eq!(arr[[3, 3]], 0.0);
        assert_eq!(arr[[4, 0]], 0.49707602339181284);
        assert_eq!(arr[[4, 1]], 0.09941520467836257);
        assert_eq!(arr[[4, 2]], 0.1695906432748538);
        assert_eq!(arr[[4, 3]], 0.49122807017543857);
        assert_eq!(arr[[4, 4]], 0.0);
        assert_eq!(arr[[5, 0]], 0.42105263157894735);
        assert_eq!(arr[[5, 1]], 0.2807017543859649);
        assert_eq!(arr[[5, 2]], 0.21052631578947367);
        assert_eq!(arr[[5, 3]], 0.4152046783625731);
        assert_eq!(arr[[5, 4]], 0.28654970760233917);
        assert_eq!(arr[[5, 5]], 0.0);
        assert_eq!(arr[[6, 0]], 0.3216374269005848);
        assert_eq!(arr[[6, 1]], 0.30994152046783624);
        assert_eq!(arr[[6, 2]], 0.23976608187134502);
        assert_eq!(arr[[6, 3]], 0.3157894736842105);
        assert_eq!(arr[[6, 4]], 0.3157894736842105);
        assert_eq!(arr[[6, 5]], 0.23976608187134502);
        assert_eq!(arr[[6, 6]], 0.0);
        assert_eq!(arr[[7, 0]], 0.27485380116959063);
        assert_eq!(arr[[7, 1]], 0.4152046783625731);
        assert_eq!(arr[[7, 2]], 0.34502923976608185);
        assert_eq!(arr[[7, 3]], 0.26900584795321636);
        assert_eq!(arr[[7, 4]], 0.42105263157894735);
        assert_eq!(arr[[7, 5]], 0.34502923976608185);
        assert_eq!(arr[[7, 6]], 0.24561403508771928);
        assert_eq!(arr[[7, 7]], 0.0);
        assert_eq!(arr[[8, 0]], 0.4152046783625731);
        assert_eq!(arr[[8, 1]], 0.15789473684210525);
        assert_eq!(arr[[8, 2]], 0.029239766081871343);
        assert_eq!(arr[[8, 3]], 0.4093567251461988);
        assert_eq!(arr[[8, 4]], 0.16374269005847952);
        assert_eq!(arr[[8, 5]], 0.2046783625730994);
        assert_eq!(arr[[8, 6]], 0.23391812865497075);
        assert_eq!(arr[[8, 7]], 0.3391812865497076);
        assert_eq!(arr[[8, 8]], 0.0);
        assert_eq!(arr[[9, 0]], 0.43859649122807015);
        assert_eq!(arr[[9, 1]], 0.12280701754385964);
        assert_eq!(arr[[9, 2]], 0.1111111111111111);
        assert_eq!(arr[[9, 3]], 0.4327485380116959);
        assert_eq!(arr[[9, 4]], 0.1286549707602339);
        assert_eq!(arr[[9, 5]], 0.22807017543859648);
        assert_eq!(arr[[9, 6]], 0.2573099415204678);
        assert_eq!(arr[[9, 7]], 0.36257309941520466);
        assert_eq!(arr[[9, 8]], 0.10526315789473684);
        assert_eq!(arr[[9, 9]], 0.0);

        assert_eq!(
            tree.leaf_nodes().unwrap(),
            [
                Taxon("T1".to_string()),
                Taxon("T10".to_string()),
                Taxon("T2".to_string()),
                Taxon("T3".to_string()),
                Taxon("T4".to_string()),
                Taxon("T5".to_string()),
                Taxon("T6".to_string()),
                Taxon("T7".to_string()),
                Taxon("T8".to_string()),
                Taxon("T9".to_string())
            ]
        );
    }

    #[test]
    fn test_set_lengths() {
        
        // Load in the default tree from test.tree
        // Then calculate the distance matrix as arr_normal
        let mut tree_normal = PDM::default();
        let _ = tree_normal.load_from_newick_path("tests/test.tree");
        let (_taxon, arr_normal) = tree_normal.matrix(false).unwrap();
        let tree_normal_length = tree_normal.length();

        // Load in the modified tree with all branch lengths = 7
        // Then calculate the distance matrix as arr_bl7
        let mut tree_bl7 = PDM::default();
        let _ = tree_bl7.load_from_newick_path("tests/test_bl_7.tree");
        let (_taxon, arr_bl7) = tree_bl7.matrix(false).unwrap();
        let tree_bl7_length = tree_bl7.length();
        
        // Use the default tree and modify all branch lengths to be 7
        // Then calculate the distance matrix as arr_modified_7
        let _ = tree_normal.update_all_edge_lengths(Edge(7.0));
        let (_taxon, arr_modified_7) = tree_normal.matrix(false).unwrap();
        let tree_modified_7_length: Edge = tree_normal.length();
        
        println!("\n\nNormal");
        println!("{:?}", arr_normal);
        println!("{:?}", tree_normal_length);

        println!("Branch lengths 7 (from file)");
        println!("{:?}", arr_bl7);
        println!("{:?}", tree_bl7_length);

        println!("Branch lengths 7 (modified)");
        println!("{:?}", arr_modified_7);
        println!("{:?}", tree_modified_7_length);

        // Somehow the modified tree is 7 units longer that the one read from file.
        assert_eq!(arr_bl7, arr_modified_7);
        assert_eq!(tree_bl7_length.0, tree_modified_7_length.0);
      
    }


    #[test]
    fn test_update_edge_lengths() {
        let mut tree = PDM::default();
        
        let taxon_b = Taxon("B".to_string());
        let taxon_c = Taxon("C".to_string());
        let taxon_d = Taxon("D".to_string());
        
        let root_node = tree.add_node(None).unwrap();
        let node_a = tree.add_node(None).unwrap();
        let node_b = tree.add_node(Some(&taxon_b)).unwrap();
        let node_c = tree.add_node(Some(&taxon_c)).unwrap();
        let node_d = tree.add_node(Some(&taxon_d)).unwrap();
        
        tree.add_edge(root_node, node_a, Edge(1.0));
        tree.add_edge(root_node, node_b, Edge(2.0));
        tree.add_edge(node_a, node_c, Edge(3.0));
        tree.add_edge(node_a, node_d, Edge(7.0));
        
        tree.compute_row_vec().unwrap();
        
        let dist_b_to_c_before = tree.distance(&taxon_b, &taxon_c, false);
        
        let node_ids = vec![node_a, node_c];
        let lengths = vec![Edge(11.0), Edge(12.0)];
        
        let _ = tree.update_edge_lengths(&node_ids, &lengths);

        let dist_b_to_c_after = tree.distance(&taxon_b, &taxon_c, false);
        
        assert_eq!(dist_b_to_c_before, 6.0);
        assert_eq!(dist_b_to_c_after, 25.0);

    }


}
