import argparse

from project.coxeter import CoxeterGraph
from project.managers import BraidQueue


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--type',
        type=str,
        help='type of (twisted) Coxeter system',
        choices=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                 '2A', '2B', '2C', '2D', '2E', '2F', '2G',
                 'A~', 'B~', 'C~', 'D~', 'E~', 'F~', 'G~',
                 '2A~', '2B~', '2C~', '2E~',
                 'rA~', 'fA~']
    )
    parser.add_argument(
        '--rank',
        type=int,
        help='rank of Coxeter system'
    )
    parser.add_argument(
        '--verify',
        dest='verify',
        action='store_true',
        help='verify relations after constructing them (required if limit is set)'
    )
    parser.add_argument(
        '--verbosity',
        type=int,
        help='amount of output to print',
        default=2,
        choices=[0, 1, 2, 3]
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='only look for relations of length up to this (optional) limit'
    )
    return parser.parse_args()


def get_coxeter_graph(coxeter_type, rank):
    coxeter_graph_constructor_dict = {
        'A': CoxeterGraph.A,
        'B': CoxeterGraph.B,
        'C': CoxeterGraph.B,
        'D': CoxeterGraph.D,
        'E': CoxeterGraph.E,
        'F': CoxeterGraph.F,
        'G': CoxeterGraph.G,
        'H': CoxeterGraph.H,
        '2A': CoxeterGraph.A_twist,
        '2B': CoxeterGraph.B_twist,
        '2C': CoxeterGraph.B_twist,
        '2D': CoxeterGraph.D_twist,
        '2E': CoxeterGraph.E_twist,
        '2F': CoxeterGraph.F_twist,
        '2G': CoxeterGraph.G_twist,
        'A~': CoxeterGraph.A_tilde,
        'B~': CoxeterGraph.B_tilde,
        'C~': CoxeterGraph.C_tilde,
        'D~': CoxeterGraph.D_tilde,
        'E~': CoxeterGraph.E_tilde,
        'F~': CoxeterGraph.F_tilde,
        'G~': CoxeterGraph.G_tilde,
        '2A~': CoxeterGraph.A_tilde_twist,
        '2B~': CoxeterGraph.B_tilde_twist,
        '2C~': CoxeterGraph.C_tilde_twist,
        '2E~': CoxeterGraph.E_tilde_twist,
        'rA~': CoxeterGraph.A_tilde_rotate,
        'fA~': CoxeterGraph.A_tilde_flip
    }
    return coxeter_graph_constructor_dict[coxeter_type](rank)


def solve(coxeter_type, rank, verbosity, verify, limit):
    if limit is not None and not verify:
        print('Error: `--verify` is required if `--limit` is given')
        return
    try:
        g = get_coxeter_graph(coxeter_type, rank)
    except:
        print('Invalid type and rank: (%s, %s)' % (coxeter_type, rank))
    else:
        q = BraidQueue(g, verbose_level=verbosity)
        q.go(limit=limit, verify=verify)


if __name__ == '__main__':
    args = get_arguments()
    solve(args.type, args.rank, args.verbosity, args.verify, args.limit)
