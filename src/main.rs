use crate::tree::ArenaTree;

mod node;
mod tree;
mod types;
mod util;


fn main() {
    let mut tree = ArenaTree::default();
    tree.load_from_newick_path(&"/Users/aaron/tmp/500_taxa.tree".to_string());

    // tree.load_from_newick_path(&"/Users/aaron/git/PhyloDM/src/newick.txt".to_string());
    let arr = tree.dm(&true);
    println!("{:?}", arr);
    println!("Length: {:?}", tree.length());
    return
}
