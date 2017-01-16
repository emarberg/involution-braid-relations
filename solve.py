import argparse

from project.algebra import CoxeterGraph
from project.managers import BraidQueue


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--type',
        type=str,
        help='type of (twisted) Coxeter system',
        choices=['A', 'B', 'D', 'E', 'F', 'G', 'H',
                 '2A', '2B', '2D', '2E', '2F', '2G',
                 'A~', 'B~', 'C~', 'D~', 'E~', 'F~', 'G~']
    )
    parser.add_argument(
        '--rank',
        type=int,
        help='rank of Coxeter system'
    )
    parser.add_argument(
        '--verify',
        dest='do_sanity_check',
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
    # replace 'X~' with 'X_tilde'
    if coxeter_type.endswith('~'):
        coxeter_type = coxeter_type[:-1] + '_tilde'
    # reverse '2X' to be 'X2'
    elif len(coxeter_type) != 1:
        coxeter_type = ''.join(reversed(coxeter_type))

    return getattr(CoxeterGraph, coxeter_type)(rank)


def solve(coxeter_type, rank, verbosity, do_sanity_check, limit):
    try:
        g = get_coxeter_graph(coxeter_type, rank)
    except:
        print('Invalid type and rank: (%s, %s)' % (coxeter_type, rank))
        return
    else:
        q = BraidQueue(g, verbose_level=verbosity)
        q.go(do_sanity_check=do_sanity_check, limit=limit)


if __name__ == '__main__':
    args = get_arguments()
    solve(args.type, args.rank, args.verbosity, args.do_sanity_check, args.limit)
