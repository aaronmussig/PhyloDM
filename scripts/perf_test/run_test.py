import argparse
import random
import sys

from dendropy.simulate import treesim

from phylodm.pdm import PDM


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('n_tax', type=int, help='number of taxon to simulate')
    parser.add_argument('method', type=str, help='phylodm or dendropy')
    args = parser.parse_args(args)

    tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5,
                                    ntax=args.n_tax, rng=random.Random(42))

    if args.method == 'phylodm':
        PDM.get_from_dendropy(tree, method='pd')
        sys.exit(0)
    elif args.method == 'dendropy':
        tree.phylogenetic_distance_matrix()
        sys.exit(0)

    sys.exit(1)


if __name__ == '__main__':
    main()
