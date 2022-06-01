import dendropy
import typer


def main(path_tree: str):
    tree = dendropy.Tree.get_from_path(path_tree, schema='newick')
    tree.phylogenetic_distance_matrix()


if __name__ == "__main__":
    typer.run(main)
