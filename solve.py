import argparse
import logging

from project.coxeter import CoxeterGraph
from project.managers import BraidQueue

logger = logging.getLogger(__name__)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--type',
        type=str,
        help='type of (twisted) Coxeter system',
        choices=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                 '2A', '2B', '2C', '2D', '2E', '2F', '2G',
                 'A~', 'B~', 'C~', 'D~', 'E~', 'F~', 'G~',
                 '2A~', '2B~', '2C~', '2D~', '2E~',
                 'rA~', 'fA~', 'hD~', 'sD~']
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
        '--limit',
        type=int,
        default=None,
        help='only look for relations of length up to this (optional) limit'
    )
    parser.add_argument(
        '--log',
        dest='loglevel',
        type=str,
        help="desired level of logging ('DEBUG', 'INFO', 'WARNING', 'ERROR')",
        default='INFO'
    )
    parser.add_argument(
        '--s',
        type=str,
        help="index of generator s for abbreviated calculation",
        default=None
    )
    parser.add_argument(
        '--t',
        type=str,
        help="index of generator t for abbreviated calculation",
        default=None
    )
    return parser.parse_args()


def get_coxeter_graph(coxeter_type, rank, s, t):
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
        '2D~': CoxeterGraph.D_tilde_twist,
        '2E~': CoxeterGraph.E_tilde_twist,
        'rA~': CoxeterGraph.A_tilde_rotate,
        'fA~': CoxeterGraph.A_tilde_flip,
        'hD~': CoxeterGraph.D_tilde_half_twist,
        'sD~': CoxeterGraph.D_tilde_small_twist,
    }
    # subtract 1 in affine case since affine constructors return systems of rank n + 1
    if coxeter_type.endswith('~'):
        rank -= 1
    g = coxeter_graph_constructor_dict[coxeter_type](rank)
    # convert strings s, t to integers if necessary
    if s is not None and s not in g.generators:
        s = int(s)
    if t is not None and t not in g.generators:
        t = int(t)
    return g, s, t


def setup_logging(loglevel):
    loglevel = args.loglevel
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)


def solve(coxeter_type, rank, verify, limit, s=None, t=None):
    if limit is not None and not verify:
        print('Error: `--verify` is required if `--limit` is given')
        return
    if s == t and t is not None:
        print('Error: `inputs for `--s` and `--t` must be distinct')
        return
    try:
        g, s, t = get_coxeter_graph(coxeter_type, rank, s, t)
    except:
        print('Invalid type and rank: (%s, %s)' % (coxeter_type, rank))
    else:
        q = BraidQueue(g, s, t)
        q.go(limit=limit, verify=verify)
        summarize_parabolic_config(g, s, t, q)


def summarize_parabolic_config(coxeter_graph, s, t, q):
    if None in [s, t]:
        return
    if q.neighborhood == set(coxeter_graph.generators):
        print('')
        print('Calculation shows that (W, S, J, * | s, t) is a never parabolic system for:')
        print('  s = %s' % s)
        print('  t = %s' % t)
        print('No proper subset J of S satisfies the required conditions.')
    else:
        print('')
        print('Calculation shows that (W, S, J, * | s, t) is a parabolic system for:')
        print('  s = %s' % s)
        print('  t = %s' % t)
        print('  J = %s' % q.neighborhood)
        print('  S = %s' % set(coxeter_graph.generators))


if __name__ == '__main__':
    args = get_arguments()
    setup_logging(args.loglevel)
    solve(args.type, args.rank, args.verify, args.limit, args.s, args.t)
