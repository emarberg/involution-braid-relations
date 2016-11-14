import argparse

from project.algebra import CoxeterGraph
from project.managers import SolverQueue


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--type',
        type=str,
        help='type of (twisted) Coxeter system',
        choices=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', '2A', '2D', '2E', '2F', '2G']
    )
    parser.add_argument(
        '--rank',
        type=int,
        help='rank of Coxeter system'
    )
    parser.add_argument(
        '--sanity-check',
        dest='do_sanity_check',
        action='store_true',
        help='verify relations after constructing them (an optional sanity check)'
    )
    parser.add_argument(
        '--verbosity',
        type=int,
        help='amount of output to print',
        default=2,
        choices=[0, 1, 2, 3]
    )
    return parser.parse_args()


def main(coxeter_type, rank, verbosity, do_sanity_check):
    # reverse '2A', '2D', etc., to be 'A2', 'D2', and so on
    try:
        reversed_coxeter_type = ''.join(reversed(coxeter_type))
        g = getattr(CoxeterGraph, reversed_coxeter_type)(rank)
    except:
        print('Invalid type and rank: (%s, %s)' % (coxeter_type, rank))
    else:
        q = SolverQueue(g, verbose_level=verbosity)
        q.go()
        if do_sanity_check:
            q._print_status('')
            q._print_status('')
            q._print_status('Step 3. Verifying minimal relations (optional sanity check).')
            q._print_status('')
            q.sanity_check()


if __name__ == '__main__':
    args = get_arguments()
    main(args.type, args.rank, args.verbosity, args.do_sanity_check)
