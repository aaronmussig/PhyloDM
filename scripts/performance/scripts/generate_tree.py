import random

import typer
from dendropy.simulate import treesim


def main(n_taxa: int, output_path: str):
    tree = treesim.birth_death_tree(
        birth_rate=1.0,
        death_rate=0.2,
        num_extant_tips=n_taxa,
        rng=random.Random(42)
    )
    with open(output_path, 'w') as f:
        f.write(tree.as_string(schema='newick')[5:])


if __name__ == "__main__":
    typer.run(main)
