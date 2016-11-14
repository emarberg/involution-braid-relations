import time
import itertools

from algebra import (
    CoxeterTransform,
    CoxeterWord,
    Monomial,
    Polynomial,
    QuadraticNumber,
    RationalNumber,
    Root,
    RootTransform
)
from utils import reverse_tuple


class ConstraintsManager:
    def __init__(self):
        self.pivots = []
        # list of linear Polynomials which must be <= 0
        self.nonpositive_constraints = set()
        # list of linear Polynomials which must be == 0
        self.zero_constraints = set()
        # list of Roots which must be != 0
        self.nonzero_constraints = set()
        # list of quadratic Polynomials which must be == 0
        self.quadratic_constraints = set()

    def __eq__(self, other):
        return \
            ConstraintsManager == type(other) and \
            self.pivots == other.pivots and \
            self.nonpositive_constraints == other.nonpositive_constraints and \
            self.zero_constraints == other.zero_constraints and \
            self.nonzero_constraints == other.nonzero_constraints and \
            self.quadratic_constraints == other.quadratic_constraints

    def copy(self):
        other = ConstraintsManager()
        other.pivots = self.pivots.copy()
        other.nonpositive_constraints = self.nonpositive_constraints.copy()
        other.zero_constraints = self.zero_constraints.copy()
        other.nonzero_constraints = self.nonzero_constraints.copy()
        other.quadratic_constraints = self.quadratic_constraints.copy()
        return other

    def add_zero_constraint(self, constraint):
        if type(constraint) == Root:
            for v in constraint.coefficients.values():
                self.add_zero_constraint(v)
            return
        elif type(constraint) in [int, RationalNumber, QuadraticNumber]:
            constraint = Polynomial(constraint)
        elif type(constraint) != Polynomial:
            raise Exception('Constraint cannot be of type `%s`' % type(constraint))

        if constraint.is_degree_one():
            self.zero_constraints.add(constraint)
        else:
            self.quadratic_constraints.add(constraint)

    def add_nonpositive_constraint(self, constraint):
        if type(constraint) == Root:
            for v in constraint.coefficients.values():
                self.add_nonpositive_constraint(v)
            return
        elif type(constraint) in [int, RationalNumber, QuadraticNumber]:
            constraint = Polynomial(constraint)
        elif type(constraint) != Polynomial:
            raise Exception('Constraint cannot be of type `%s`' % type(constraint))

        if -constraint in self.nonpositive_constraints:
            self.add_zero_constraint(constraint)
        elif not any(constraint <= f for f in self.nonpositive_constraints):
            self.nonpositive_constraints.add(constraint)

    def add_nonzero_constraint(self, root):
        self.nonzero_constraints.add(root)

    def row_reduce_zero_constraints(self):
        degenerate_constraints = set()
        while self.zero_constraints:
            row = self.zero_constraints.pop()
            variables = row.get_variables()
            if not variables:
                if row != 0:
                    degenerate_constraints.add(row)
                continue

            var = Monomial({variables.pop(): 1})
            substitution = -row / row[var]
            substitution[var] = 0

            for i in range(len(self.pivots)):
                prev_var, prev_subst = self.pivots[i]
                c = prev_subst[var]
                if c != 0:
                    new_subst = prev_subst + c*substitution
                    # note that since new_subst was assigned to newly created in previous line,
                    # the following deletion does not affect any other existing polynomial.
                    new_subst[var] = 0
                    self.pivots[i] = (prev_var, new_subst)

            self.pivots.append((var, substitution))
            self.row_reduce_by_pivot(var, substitution)

        self.zero_constraints = degenerate_constraints

    def row_reduce_by_pivot(self, var, substitution):
        self.zero_constraints = {
            c for c in self.reduced_row_generator(var, substitution, self.zero_constraints)
            if c != 0
        }
        self.nonpositive_constraints = {
            c for c in self.reduced_row_generator(var, substitution, self.nonpositive_constraints)
            if not (c <= 0)
        }

    @classmethod
    def reduced_row_generator(cls, var, substitution, old_constraints):
        for constraint in old_constraints:
            c = constraint[var]
            if c != 0:
                new_constraint = constraint + c*substitution
                # note that since new_subst was assigned to newly created in previous line,
                # the following deletion does not affect any other existing polynomial.
                new_constraint[var] = 0
                yield new_constraint
            else:
                yield constraint

    def remove_vacuous_constraints(self):
        self.zero_constraints = {
            f for f in self.zero_constraints
            if not (f == 0)
        }
        self.nonpositive_constraints = {
            f for f in self.nonpositive_constraints
            if not (f <= 0) and not any(f <= g and f != g for g in self.nonpositive_constraints)
        }
        # if 10 < len(self.nonpositive_constraints):
        # self.nonpositive_constraints = [
        #     f for f in self.nonpositive_constraints
        #     if not any(f <= g and f != g for g in self.nonpositive_constraints)
        # ]
        self.nonzero_constraints = {
            r for r in self.nonzero_constraints
            if not any(v < 0 or 0 < v for v in r.coefficients.values())
        }
        self.quadratic_constraints = {
            f for f in self.quadratic_constraints if not (f == 0)
        }

    def __repr__(self):
        s = '\nconstraints:\n'
        i = 1

        def pad(j):
            return (3 - len(str(j)))*' ' + str(j)

        for var, substitution in self.pivots:
            s += '%s. %s = %s\n' % (pad(i), var, substitution)
            i += 1

        for c in self.nonpositive_constraints:
            s += '%s. 0 >= %s\n' % (pad(i), c)
            i += 1

        for c in self.zero_constraints:
            s += '%s. 0 == %s\n' % (pad(i), c)
            i += 1

        for c in self.quadratic_constraints:
            s += '%s. 0 == %s\n' % (pad(i), c)
            i += 1

        for c in self.nonzero_constraints:
            s += '%s. 0 != %s\n' % (pad(i), c)
            i += 1

        if s == '\nconstraints:\n':
            return ''
        else:
            return s

    def is_valid(self):
        if any(0 < f for _, f in self.pivots) or \
           any(0 < f for f in self.nonpositive_constraints) or \
           any(0 < f or f < 0 for f in self.zero_constraints) or \
           any(0 < f or f < 0 for f in self.quadratic_constraints) or \
           any(r == 0 for r in self.nonzero_constraints):
            return False
        return True


class BraidSolver:
    def __init__(self, coxeter_graph, s, t):
        if s in coxeter_graph.generators and t in coxeter_graph.generators:
            self.graph = coxeter_graph
            self.s = s
            self.t = t
            self.sigma = RootTransform(coxeter_graph)
            self.word_s = CoxeterWord(coxeter_graph)
            self.word_t = CoxeterWord(coxeter_graph)
            self.constraints = ConstraintsManager()
        else:
            raise Exception('Invalid inputs `s = %s`, `t = %s` in BraidSolver' % (s, t))

    def __eq__(self, other):
        return \
            BraidSolver == type(other) and \
            self.graph == other.graph and \
            self.sigma == other.sigma and \
            self.word_s.left_action == other.word_s.left_action and \
            self.word_t.left_action == other.word_t.left_action and \
            self.constraints == other.constraints

    def __len__(self):
        return len(self.sigma)

    def __repr__(self):
        s = '\n'
        s += 'Solver State:\n'
        s += '***************\n'
        s += 's = %s, word_s = %s\n' % (self.s, self.word_s)
        s += 't = %s, word_t = %s\n' % (self.t, self.word_t)
        s += '\n'
        s += 'sigma = %s' % self.sigma
        s += '\n'
        s += '      unconditional descent: %s\n' % str(self.get_unconditional_descent())
        s += ' strong conditional descent: %s\n' % str(self.get_strong_conditional_descent())
        s += '  weak conditional descents: %s\n' % list(self.get_weak_conditional_descents())
        s += str(self.constraints)
        s += '***************\n'
        return s

    def copy(self):
        other = BraidSolver(self.graph, self.s, self.t)
        other.sigma = self.sigma.copy()
        other.word_s = self.word_s.copy()
        other.word_t = self.word_t.copy()
        other.constraints = self.constraints.copy()
        return other

    def clear_constraints(self):
        self.constraints = ConstraintsManager()

    def get_semiorder(self, is_fixer=True):
        m = self.graph.get_order(self.s, self.t)
        if m % 2 != 0:
            return (m+1)//2
        elif is_fixer:
            return m//2 + 1
        else:
            return m//2

    def get_unconditional_descent(self):
        descents_to_avoid = self.word_s.left_descents | self.word_t.left_descents
        unconditional = \
            list(self.sigma.unconditional_descents.intersection(descents_to_avoid)) + \
            list(self.sigma.unconditional_descents.difference(descents_to_avoid))
        if unconditional:
            return unconditional[0]
        return None

    def get_strong_conditional_descent(self):
        for i in self.sigma:
            root = Root(self.graph)
            for j, f in self.sigma[i]:
                if f <= 0 or any(f <= g for g in self.constraints.nonpositive_constraints):
                    root += Root(self.graph, j, f)
            if root != 0:
                return i, root
        return None

    def get_weak_conditional_descents(self):
        descents_to_avoid = self.word_s.left_descents | self.word_t.left_descents
        return self.sigma.weak_conditional_descents.difference(descents_to_avoid)

    def is_quadratic_constraint_factorable(self):
        if self.constraints.quadratic_constraints:
            c = next(iter(self.constraints.quadratic_constraints))
            return c.is_factorable()
        else:
            return False

    def branch(self, should_filter=True):
        t0 = time.time()
        children, label = self._get_children()

        t1 = time.time()
        for child in children:
            child.reduce()

        t2 = time.time()
        children = [child for child in children if not should_filter or child.is_valid()]

        t3 = time.time()
        descr = '\n'
        descr += 'BRANCHING: %s\n' % label
        descr += '  CONSTRUCT: %s seconds\n' % (t1-t0)
        descr += '  REDUCTION: %s seconds\n' % (t2-t1)
        descr += '  VALIDITY : %s seconds' % (t3-t2)
        return children, descr

    def _get_children(self):
        children = []

        if len(self.sigma) == 0:
            return self._get_children_first_iteration(), 'first iteration'

        if self.is_quadratic_constraint_factorable():
            return self._get_children_from_quadratic_constraint(), 'reducing quadratic constraint'

        unconditional = self.get_unconditional_descent()
        if unconditional:
            children = self._get_children_from_unconditional_descents(unconditional)
            return children, 'unconditional descent'

        strong = self.get_strong_conditional_descent()
        if strong:
            i, root = strong
            children = self._get_children_from_strong_conditional_descents(i, root)
            return children, 'strong conditional descent'

        weak = self.get_weak_conditional_descents()
        if weak:
            children = self._get_children_from_weak_conditional_descents(weak)
            return children, 'weak conditional descents'

        if self.sigma.is_constant() and not self.sigma.is_complete():
            return self._get_children_from_new_descent(), 'new descent'

        if self.sigma.is_constant() and not self.sigma.is_identity():
            return [], 'empty'

        raise Exception('Bad case: %s' % self)

    def _get_children_first_iteration(self):
        gens = [self.s, self.t]
        alpha = Root(self.graph, self.graph.star(self.s))
        beta = Root(self.graph, self.graph.star(self.t))

        fixer = BraidSolver(self.graph, self.s, self.t)
        fixer.sigma = RootTransform(self.graph, {self.s: alpha, self.t: beta})
        for i in range(self.get_semiorder(is_fixer=True)):
            fixer.word_s.extend_left(gens[i % 2])
            fixer.word_t.extend_left(gens[(i+1) % 2])

        transposer = BraidSolver(self.graph, self.s, self.t)
        transposer.sigma = RootTransform(self.graph, {self.s: beta, self.t: alpha})
        for i in range(self.get_semiorder(is_fixer=False)):
            transposer.word_s.extend_left(gens[i % 2])
            transposer.word_t.extend_left(gens[(i+1) % 2])

        return [fixer, transposer]

    def _get_children_from_quadratic_constraint(self):
        children = []
        constraint = next(iter(self.constraints.quadratic_constraints))
        for factor in constraint.get_factors():
            child = self.copy()
            child.constraints.add_zero_constraint(factor)
            children.append(child)
        return children

    def _get_children_from_unconditional_descents(self, descent):
        current = self._get_children_from_next_unconditional_descent(descent)
        children = []
        while current:
            new = []
            for state in current:
                state.reduce()
                if state.is_valid():
                    descent = state.get_unconditional_descent()
                    if descent:
                        new += state._get_children_from_next_unconditional_descent(descent)
                    else:
                        children.append(state)
            current = new
        return children

    def _get_children_from_next_unconditional_descent(self, descent):
        children = []
        if self.sigma[descent].is_constant():
            commutes = self.sigma[descent] == Root(self.graph, self.graph.star(descent), -1)
            children.append(self._branch_from_descent(descent, translate=commutes))
        else:
            children.append(self._branch_from_descent(descent, True))
            children.append(self._branch_from_descent(descent, False))
        return children

    def _get_children_from_strong_conditional_descents(self, descent, nonzero_root):
        constraints = [(j, f) for (j, f) in nonzero_root]

        child_a = self._branch_from_descent(descent, True)
        child_a.constraints.add_nonzero_constraint(nonzero_root)

        child_b = self._branch_from_descent(descent, False)
        child_b.constraints.add_nonzero_constraint(nonzero_root)

        child_c = self.copy()
        for _, f in constraints:
            child_c.constraints.add_zero_constraint(f)

        return [child_a, child_b, child_c]

    def _get_children_from_weak_conditional_descents(self, descents):
        children = []
        for i in descents:
            children.append(self._branch_from_descent(i, True))
            children.append(self._branch_from_descent(i, False))
            children.append(self._branch_from_descent(i, descent=False))
        return children

    def _branch_from_descent(self, i, translate=True, descent=True):
        beta = self.sigma[i]
        new = self.copy()
        if not descent:
            new.constraints.add_nonpositive_constraint(-beta)
        else:
            new.constraints.add_nonpositive_constraint(beta)
            new.word_s.extend_left(i)
            new.word_t.extend_left(i)

            alpha = Root(new.graph, new.graph.star(i))
            if translate:
                new.constraints.add_zero_constraint(beta + alpha)
                new.sigma = new.sigma * i
            else:
                new.constraints.add_nonzero_constraint(beta + alpha)
                new.sigma = self.graph.star(i) * new.sigma * i
        return new

    def _get_children_from_new_descent(self):
        g = self.graph
        children = []

        candidates = [
            i for i in g.generators
            if i not in self.sigma and any(g.get_order(i, j) != 2 for j in self.sigma)
        ]
        for i in candidates:
            # add child with new weak descent i
            child = self.copy()
            child.clear_constraints()
            child.sigma[i] = sum([-Polynomial({j: 1})*Root(g, j) for j in g.generators])
            for j in child.sigma:
                child.constraints.add_zero_constraint(
                    child.sigma[i].eval_bilinear(child.sigma[j]) - g.eval_bilinear(i, j)
                )
            children.append(child)

        # add child with no new weak descents
        child = self.copy()
        for i in g.generators:
            if i not in child.sigma:
                child.sigma[i] = Root(g, i)
        children.append(child)
        return children

    def is_valid(self):
        return \
            self.are_descents_valid() and \
            self.word_s.is_reduced and self.word_t.is_reduced and \
            self.is_sigma_valid() and \
            self.constraints.is_valid()

    def are_descents_valid(self):
        """
        Returns False if word_s and word_t share a right descent, or have (distinct)
        right descents with product of order greater than m_st. These cases may be excluded
        by induction.
        """
        if not self.word_s.right_descents.isdisjoint(self.word_t.right_descents):
            return False
        if any(self.graph.get_order(self.s, self.t) < self.graph.get_order(u, v)
           for u in self.word_s.right_descents
           for v in self.word_t.right_descents):
                return False
        return True

    def is_sigma_valid(self):
        # invalid if sigma sends any root to 0 or non-negative/positive combination of simple roots
        if any(not root.is_valid() for root in self.sigma.values()):
            return False
        # if sigma sends all roots to positive roots, and no variables remain
        if self.sigma.is_constant() and self.sigma.is_complete() and self.sigma.is_positive():
            # invalid if sigma is not trivial
            if not self.sigma.is_identity():
                return False
            # invalid if word_s and word_t are not involution words for the same element
            if not self.is_valid_involution_word_pair():
                return False
        return True

    def is_valid_involution_word_pair(self):
        x = self.word_s.to_involution_word()
        y = self.word_t.to_involution_word()
        return x.is_reduced and y.is_reduced and x.left_action == y.left_action

    def is_final(self):
        return self.is_valid() and self.sigma.is_identity()

    def reduce(self):
        self.constraints.row_reduce_zero_constraints()
        self.simplify_sigma()
        zeros = set()
        for f in self.constraints.nonpositive_constraints:
            if 0 <= f:
                zeros.update(f.get_variables())
        if zeros:
            self.set_variables_to_zero(zeros)
            self.constraints.row_reduce_zero_constraints()
            self.simplify_sigma()
        self.constraints.remove_vacuous_constraints()

    def simplify_sigma(self):
        for var, substitution in self.constraints.pivots:
            for i in self.sigma:
                self.sigma[i] = self.sigma[i].set_variable(var, substitution)

            self.constraints.nonzero_constraints = {
                r.set_variable(var, substitution) for r in self.constraints.nonzero_constraints
            }
            self.constraints.quadratic_constraints = {
                f.set_variable(var, substitution) for f in self.constraints.quadratic_constraints
            }
        self.constraints.pivots = []

    def set_variables_to_zero(self, variables):
        self.constraints.nonpositive_constraints = {
            f.set_variables_to_zero(variables) for f in self.constraints.nonpositive_constraints
        }
        self.constraints.zero_constraints = {
            f.set_variables_to_zero(variables) for f in self.constraints.zero_constraints
        }
        self.constraints.nonzero_constraints = {
            r.set_variables_to_zero(variables) for r in self.constraints.nonzero_constraints
        }
        self.constraints.quadratic_constraints = {
            f.set_variables_to_zero(variables) for f in self.constraints.quadratic_constraints
        }
        for i in self.sigma:
            self.sigma[i] = self.sigma[i].set_variables_to_zero(variables)


class SolverQueue:

    VERBOSE_LEVEL_NONE = 0
    VERBOSE_LEVEL_LOW = 1
    VERBOSE_LEVEL_MEDIUM = 2
    VERBOSE_LEVEL_HIGH = 3

    def __init__(self, coxeter_graph, s=None, t=None, verbose_level=VERBOSE_LEVEL_MEDIUM):
        self.graph = coxeter_graph
        if not s or not t:
            self.queue = [
                BraidSolver(coxeter_graph, s, t)
                for s in coxeter_graph.generators
                for t in coxeter_graph.generators
                if s < t
            ]
        else:
            self.queue = [BraidSolver(coxeter_graph, s, t)]
        self.final = set()
        self.minimal = None
        self.verbose_level = verbose_level

    def __len__(self):
        return len(self.queue)

    def _update(self, children):
        added = False
        for child in children:
            if child in self.queue:
                continue
            self._print_verbose(child)
            if child.is_final():
                self._add_final(child)
            else:
                self._insert(child)
                added = True
        if not added:
            self._print_verbose('\n(no states added)\n')

    def _insert(self, child):
        i = 0
        while i < len(self.queue) and len(self.queue[i]) < len(child):
            i += 1
        self.queue.insert(i, child)

    def _root_multiplicities(self):
        d = {}
        for state in self.queue:
            l = len(state.sigma)
            d[l] = d.get(l, 0) + 1
        if not d:
            return '0'
        else:
            return ' '.join([str(ell) + '^' + str(mul) for ell, mul in sorted(d.items())])

    def _word_multiplicities(self):
        e = {}
        for state in self.queue:
            l = len(state.word_s)
            e[l] = e.get(l, 0) + 1
        if not e:
            return '0'
        else:
            return ' '.join([str(ell) + '^' + str(mul) for ell, mul in sorted(e.items())])

    def _add_final(self, child):
        u = tuple(child.word_s.word)
        v = tuple(child.word_t.word)
        if u < v:
            self.final.add((u, v))
        else:
            self.final.add((v, u))

    def get_inverse_atoms(self, relations, start_word):
        start = CoxeterTransform.from_word(self.graph, reverse_tuple(start_word))
        relations = {(reverse_tuple(a), reverse_tuple(b)) for a, b in relations}
        return start.span_by_right_relations(relations)

    def are_atoms_connected(self, relations, start_word, target_word):
        target = CoxeterTransform.from_word(self.graph, reverse_tuple(target_word))
        return target in self.get_inverse_atoms(relations, start_word)

    def minimize_relations(self, final):
        necessary, redundant = self._get_necessary_relations(final)
        return self._finalize_necessary_relations(final, necessary, redundant)

    def _get_necessary_relations(self, final):
        self._print("")
        self._print("  1. Get relations which cannot be generated by all others combined.")
        necessary, redundant = [], []
        for i in range(len(final)):
            t0 = time.time()
            tag = '     * Relation %s/%s:' % (i+1, len(final))
            u, v = final[i]
            initial = [(a, b) for a, b in final[:i] if (a, b) not in redundant]
            complement = [(a, b) for a, b in (final[:i] + final[i+1:]) if (a, b) not in redundant]
            if self.are_atoms_connected(initial, u, v):
                redundant.append((u, v))
                label = "redundant"
            elif not self.are_atoms_connected(complement, u, v):
                necessary.append((u, v))
                label = "independent"
            else:
                label = "not redundant or independent"
            t1 = time.time()
            self._print("%s (%s), computation took %s seconds" % (tag, label, t1 - t0))
        return necessary, redundant

    def _finalize_necessary_relations(self, final, necessary, redundant):
        self._print("\n  2. Get smallest set of relations which generate all the others.")
        rest = [x for x in final if x not in necessary and x not in redundant]
        for i in range(len(rest)+1):
            for current in itertools.combinations(rest, i):
                candidate = necessary + list(current)
                self._print("     * Trying combination of relations:")
                for rel in candidate:
                    self._print("         %s <---> %s" % rel)
                complement = [x for x in rest if x not in current]
                if all(self.are_atoms_connected(candidate, u, v) for u, v in complement):
                    self._print("       Works!")
                    return self._finalize_relations(candidate)
                else:
                    self._print("       Does not span.")
        raise Exception('Error in `SolverQueue._finalize_necessary_relations`: returning None')

    def _finalize_relations(self, relations):
        finalized = []
        for u, v in relations:
            # find largest common initial prefix of u and v
            i = len(u)
            while 0 < i and u[:i] != v[:i]:
                i -= 1
            start_word = u[:i]
            inverse_atoms = self.get_inverse_atoms(relations, start_word)
            w = min(reverse_tuple(a.minimal_reduced_word) for a in inverse_atoms)
            finalized += [(w + u[i:], w + v[i:])]
        return finalized

    def _print_status(self, string):
        self._print(string, self.VERBOSE_LEVEL_LOW)

    def _print_verbose(self, string):
        self._print(string, self.VERBOSE_LEVEL_HIGH)

    def _print(self, string, level=VERBOSE_LEVEL_MEDIUM):
        if level <= self.verbose_level:
            print(string)

    def next(self):
        if len(self) == 0:
            self._print('Queue is empty')
        else:
            self._print_verbose('')
            self._print_verbose('-----------')
            self._print_verbose('Next state:')
            self._print_verbose('-----------')
            self._print_verbose(self.queue[0])

            children, descr = self.queue[0].branch()
            self.queue = self.queue[1:]
            self._print(descr)

            self._print_verbose('')
            self._print_verbose('-------------')
            self._print_verbose('Child states:')
            self._print_verbose('-------------')
            self._update(children)
            self._print('States in queue                  : %s' % len(self))
            self._print('Multiplicities by non-blank roots: %s' % self._root_multiplicities())
            self._print('Multiplicities by word length    : %s' % self._word_multiplicities())
            self._print('Final states                     : %s' % len(self.final))

    def go(self):
        self._print_status('Step 1. Finding sufficient relations.')

        t1 = time.time()
        while len(self) > 0:
            self.next()
        t2 = time.time()

        self._print_status('')
        self._print_status('Duration: %s seconds' % (t2-t1))

        final = sorted(self.final, key=lambda x: (len(x[0]), x))
        self._print_status('')
        self._print_status('-----------------------')
        self._print_status('Twisted Coxeter system:')
        self._print_status('-----------------------')
        self._print_status(self.graph)

        self._print_status('')
        self._print_status('---------------------')
        self._print_status('Sufficient relations:')
        self._print_status('---------------------')
        for u, v in final:
            self._print_status('%s <---> %s' % (u, v))

        self._print_status('')
        self._print_status('')
        self._print_status('Step 2. Finding minimal sufficient relations.')
        self.minimal = self.minimize_relations(final)
        self._print_status('')
        self._print_status('Duration: %s seconds' % (time.time()-t2))

        self._print_status('')
        self._print_status('-----------------------------')
        self._print_status('Minimal sufficient relations:')
        self._print_status('-----------------------------')
        for u, v in self.minimal:
            self._print_status('%s <---> %s' % (u, v))
