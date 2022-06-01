import dendropy
import typer

from phylodm import PhyloDM


def main(path_tree: str):
    tree = dendropy.Tree.get_from_path(path_tree, schema='newick')
    pdm = PhyloDM.load_from_dendropy(tree)
    pdm.dm()


if __name__ == "__main__":
    typer.run(main)
