import argparse
import sys

import dendropy

from phylodm.pdm import PDM


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('path_tree', type=str, help='path to the tree file')
    parser.add_argument('method', type=str, help='phylodm or dendropy')
    args = parser.parse_args(args)

    tree = dendropy.Tree.get_from_path(args.path_tree, schema='newick')

    if args.method == 'phylodm':
        PDM.get_from_dendropy(tree, method='pd')
        sys.exit(0)
    elif args.method == 'dendropy':
        tree.phylogenetic_distance_matrix()
        sys.exit(0)

    sys.exit(1)


if __name__ == '__main__':
    main()
