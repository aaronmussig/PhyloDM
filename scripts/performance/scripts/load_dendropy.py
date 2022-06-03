import dendropy
import typer


def main(path_tree: str):
    dendropy.Tree.get_from_path(path_tree, schema='newick')


if __name__ == "__main__":
    typer.run(main)
