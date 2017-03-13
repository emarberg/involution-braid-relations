from collections import defaultdict
import numpy as np
import time
import itertools
import logging

from project.algebra import (
    Monomial,
    Polynomial,
    QuadraticNumber,
    RationalNumber,
    Matrix
)

from project.coxeter import (
    CoxeterTransform,
    CoxeterWord,
    CoxeterVector,
    PartialTransform
)

from project.utils import (
    reverse_tuple,
    get_braids,
    InvalidInputException,
    RecurrentStateException
)

logger = logging.getLogger(__name__)


class ConstraintsManager:

    """
    Objects of this class represent a set constraints, given by equalities and inequalities
    involving (mostly linear) polynomials. The class provides methods to simplify such
    constraints, and to determine when they are impossible to simultaneously satisfy.
    """

    def __init__(self):
        # list of linear Polynomials which must be == 0
        self.linear_constraints = set()
        # list of linear Polynomials which must be <= 0
        self.nonpositive_constraints = set()
        # list of CoxeterVectors which must be != 0
        self.nonzero_constraints = set()
        self.nonnegative_indeterminates = set()

    def copy(self):
        other = ConstraintsManager()
        other.nonpositive_constraints = self.nonpositive_constraints.copy()
        other.linear_constraints = self.linear_constraints.copy()
        other.nonzero_constraints = self.nonzero_constraints.copy()
        other.nonnegative_indeterminates = self.nonnegative_indeterminates.copy()
        return other

    def is_value_nonpositive(self, f):
        """
        If the function returns True, then current nonpositive constraints imply that f
        is nonpositive for all nonnegative indeterminates. If the function returns False,
        then this property either fails or could not be determined by our naive methods.

        TODO: current implementation only involves some simple checks, and could be improved.
        For example, whether it is feasible for the desired property to fail (the case when
        the function would return False) could be solved exactly using linear programming.
        """
        if f <= 0:
            return True
        if type(f) != Polynomial:
            f = Polynomial(f)
        for g in self.nonpositive_constraints | {-h for h in self.nonnegative_indeterminates}:
            if self._compare_to_rescaled(f, g):
                return True
        return False

    @classmethod
    def _compare_to_rescaled(cls, f, g, strict=False):
        """
        Returns True if (cg - f) has all positive coefficients and (if strict is True)
        nonzero constant term, for some real number c > 0.
        """
        upper_bounds, lower_bounds = set(), {0}
        for m in f.coefficients.keys() | g.coefficients.keys():
            coeff = g[m]
            if type(coeff) == int:
                coeff = RationalNumber(coeff)

            if coeff == 0 and f[m] > 0:
                return False
            elif coeff > 0:
                lower_bounds.add(f[m] / coeff)
            elif coeff < 0:
                upper_bounds.add(f[m] / coeff)

        if len(upper_bounds) == 0:
            # have already checked that constant term of f is nonpositive if g has no constant term.
            # in strict case, need to check that f or g has constant term.
            if not strict or f.get_constant_part() != 0 or g.get_constant_part() != 0:
                return True
            else:
                return False

        lower = max(lower_bounds)
        upper = min(upper_bounds)
        if lower <= upper:
            # in strict case, determine if constant term is nonzero for some lower <= c <= upper
            a = f.get_constant_part()
            b = g.get_constant_part()
            if not strict or lower * b - a > 0 or upper * b - a > 0:
                return True
        return False

    def is_value_negative(self, f):
        """TODO: improve along the lines of is_value_nonpositive."""
        if f < 0:
            return True
        if type(f) != Polynomial:
            f = Polynomial(f)
        for g in self.nonpositive_constraints | {-h for h in self.nonnegative_indeterminates}:
            if self._compare_to_rescaled(f, g, strict=True):
                return True
        return False

    def add_zero_constraint(self, constraint):
        """
        Add input Polynomial `constraint` to linear constraints accordingly.
        If input is a CoxeterVector, do the same for all of its coefficients.
        """
        if type(constraint) == CoxeterVector:
            for v in constraint.coefficients.values():
                self.add_zero_constraint(v)
            return
        elif type(constraint) in [int, RationalNumber, QuadraticNumber]:
            constraint = Polynomial(constraint)
        elif type(constraint) != Polynomial:
            raise InvalidInputException(self, constraint, 'add_zero_constraint')

        for index in constraint.get_variables():
            self.nonnegative_indeterminates.add(Polynomial({index: 1}))

        degree = constraint.degree()
        if degree in [0, 1]:
            self.linear_constraints.add(constraint)
        else:
            raise InvalidInputException(self, constraint, 'add_zero_constraint')

    def add_nonpositive_constraint(self, constraint):
        """
        Add input Polynomial `constraint` to nonpositive constraints. If -constraint
        is already in this set, then we call `add_zero_constraint` with `constraint` as input.
        If input is a CoxeterVector, do the same for all of its coefficients.
        """
        if type(constraint) == CoxeterVector:
            for v in constraint.coefficients.values():
                self.add_nonpositive_constraint(v)
            return
        elif type(constraint) not in [int, RationalNumber, QuadraticNumber, Polynomial]:
            raise InvalidInputException(self, constraint, 'add_nonpositive_constraint')

        constraint = self._normalize_input(constraint)
        if self.is_value_nonpositive(-constraint):
            self.add_zero_constraint(constraint)
        elif not self.is_value_nonpositive(constraint):
            self.nonpositive_constraints.add(constraint)

    def _normalize_input(self, f):
        """
        Trivially normalize input by coverting to Polynomial object.

        TODO: actually rescale `f` to avoid storing duplicate constraints.
        """
        if type(f) in [int, RationalNumber, QuadraticNumber]:
            f = Polynomial(f)
        return f

    def add_nonzero_constraint(self, root):
        if type(root) == CoxeterVector:
            # TODO: normalize `root` to avoid storing duplicate constraints.
            self.nonzero_constraints.add(root)
        else:
            raise InvalidInputException(self, root, 'add_nonzero_constraint')

    def simplify(self):
        """
        Row reduce linear constaints and eliminate resulting pivot variables from all constraints.
        Returns list of (Monomial, Polynomial) pairs giving pivot variable substitutions.
        """
        variable_substitutions = []
        degenerate_constraints = set()
        while self.linear_constraints:
            row = self.linear_constraints.pop()

            # find nonzero indeterminate x_i appearing in row, or continue if row == 0
            variables = row.get_variables()
            if not variables:
                if row != 0:
                    degenerate_constraints.add(row)
                continue
            var = Monomial({variables.pop(): 1})
            substitution = -row / row[var]
            substitution[var] = 0

            # apply new variable substitution to all previous ones
            for i in range(len(variable_substitutions)):
                prev_var, subst = variable_substitutions[i]
                new_subst = subst.set_variable(var, substitution)
                variable_substitutions[i] = (prev_var, new_subst)

            variable_substitutions.append((var, substitution))
            self.apply_variable_substitution(var, substitution)

        self.linear_constraints = degenerate_constraints
        self.remove_vacuous_constraints()
        return variable_substitutions

    def apply_variable_substitution(self, var, substitution):
        def replace(f):
            return f.set_variable(var, substitution)

        # TODO: use dedicated add_xxxx_constraint methods rather than doing in place substitutions.
        # This is slightly tricky to get right, and simple implementation is usually good enough.
        self.linear_constraints = {replace(f) for f in self.linear_constraints}
        self.nonpositive_constraints = {replace(f) for f in self.nonpositive_constraints}
        self.nonzero_constraints = {replace(r) for r in self.nonzero_constraints}
        self.nonnegative_indeterminates = {replace(f) for f in self.nonnegative_indeterminates}

    def remove_vacuous_constraints(self):
        self.linear_constraints = {
            f for f in self.linear_constraints if not (f == 0)
        }

        # TODO: improve how redundant inequalities are detected in three cases below.
        # Current implementation is very simple, and will consider scalar multiples of
        # a single constraint to be different.
        self.nonpositive_constraints = {
            f for f in self.nonpositive_constraints
            if not (f <= 0) and not any(f <= g and f != g for g in self.nonpositive_constraints)
        }
        self.nonzero_constraints = {
            r for r in self.nonzero_constraints
            if not any(v < 0 or 0 < v for v in r.coefficients.values())
        }
        self.nonnegative_indeterminates = {
            f for f in self.nonnegative_indeterminates
            if not (f.is_constant() and f >= 0)
        }

    def __repr__(self):
        s = '\nconstraints:\n'
        i = 1

        def pad(j):
            return (3 - len(str(j))) * ' ' + str(j)

        for c in self.nonpositive_constraints:
            s += '%s. 0 >= %s\n' % (pad(i), c)
            i += 1
        for c in self.nonnegative_indeterminates:
            s += '%s. 0 <= %s\n' % (pad(i), c)
            i += 1

        for c in self.linear_constraints:
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
        return not (
            any(0 < f for f in self.nonpositive_constraints) or
            # we do not check f != 0 since this would be true for any nontrivial polynomial;
            # instead, we want to determine if f has all positive/negative coeffs
            any(0 < f or f < 0 for f in self.linear_constraints) or
            any(r == 0 for r in self.nonzero_constraints) or
            any(f < 0 for f in self.nonnegative_indeterminates) or
            (len(self.nonnegative_indeterminates) > 0 and sum(self.nonnegative_indeterminates) == 0)
        )


class BraidSystem:

    """
    Class for objects containing the following data:

        1. Two CoxeterWord objects of the same length.
        2. A PartialTransform object.
        3. A ConstraintsManager object.

    Each of these objects must be defined with respect to a common CoxeterGraph.
    The class implements methods for checking whether the BraidSystem is valid or
    redundant (in an appropriate technical sense) and for expanding a given system
    into a collection of more constrained systems.
    """

    def __init__(self, coxeter_graph, s, t, is_fixer=True):
        """
        The initial value of self.sigma must define a bijection from
        {alpha_s, alpha_t} to {alpha_s^*, alpha_t^*}. Which bijection is
        chosen controlled by the input `is_fixer`.
        """
        if s in coxeter_graph.generators and t in coxeter_graph.generators:
            self.constraints = ConstraintsManager()
            self.graph = coxeter_graph
            self.s = s
            self.t = t
            self.is_fixer = is_fixer
            self.word_s = CoxeterWord(coxeter_graph)
            self.word_t = CoxeterWord(coxeter_graph)
            self._extend_words()

            if self.is_fixer:
                alpha = CoxeterVector(self.graph, self.graph.star(self.s))
                beta = CoxeterVector(self.graph, self.graph.star(self.t))
            else:
                alpha = CoxeterVector(self.graph, self.graph.star(self.t))
                beta = CoxeterVector(self.graph, self.graph.star(self.s))
            self.sigma = PartialTransform(self.graph, {self.s: alpha, self.t: beta})
        else:
            raise InvalidInputException(self, (s, t))

    def _extend_words(self):
        """Helper method for initializing self.word_s and self.word_t."""
        gens = [self.s, self.t]
        for i in range(self.graph.get_semiorder(self.s, self.t, self.is_fixer)):
            self.word_s.extend_left(gens[i % 2])
            self.word_t.extend_left(gens[(i + 1) % 2])

    def __len__(self):
        return len(self.word_s)

    def __repr__(self):
        unconditional = self.get_unconditional_descent()
        conditional = self.get_conditional_descent()
        if self.is_fixer:
            a, b = 's', 't'
        else:
            a, b = 't', 's'

        s = ''
        s += 'BraidSystem data:\n'
        s += '----------------------------------------------------------------\n'
        s += 's = %s, word_s = %s\n' % (self.s, self.word_s)
        s += 't = %s, word_t = %s\n' % (self.t, self.word_t)
        s += '\n'
        s += 'initial sigma: alpha_s -> alpha_%s^*, alpha_t -> alpha_%s^*\n' % (a, b)
        s += '\n'
        s += 'sigma = %s' % self.sigma
        s += '\n'
        s += 'unconditional descent: %s\n' % unconditional
        s += 'conditional descent  : %s\n' % conditional
        s += str(self.constraints)
        s += '----------------------------------------------------------------\n'
        return s

    def copy(self):
        other = BraidSystem(self.graph, self.s, self.t, self.is_fixer)
        other.sigma = self.sigma.copy()
        other.word_s = self.word_s.copy()
        other.word_t = self.word_t.copy()
        other.constraints = self.constraints.copy()
        return other

    def _clear_constraints(self):
        self.constraints = ConstraintsManager()

    def get_children(self):
        children = []
        queue = self._get_initialized_children()
        while queue:
            child = queue.pop(0)
            if child.sigma.is_constant():
                for reduced in child.eliminate_descents():
                    if not reduced.is_redundant():
                        children.append(reduced)
                    else:
                        logger.debug("Redundant system: %s" % reduced)
            else:
                queue.extend(child.branch())
        return children

    def _get_initialized_children(self):
        """
        If self.sigma is not constant, or if self.sigma has a descent, raises an Exception.
        Otherwise, returns list of children constructed from current BraidSystem
        by replacing a single undefined value of self.sigma[i] by a generic
        CoxeterVector of the form - sum_i X_i alpha_i, so that i becomes a descent.
        If the BraidSystem has no undefined values, returns the empty list.
        """
        if len(self.sigma) == 0 or not self.sigma.is_constant() or not self.sigma.is_positive():
            raise Exception(self.sigma)  # pragma: no cover
        if self.sigma.is_complete():
            return []

        logger.debug("Constructing children with new descent.")
        g = self.graph
        # exclude the descents which will not give rise to valid states
        candidate_descents = [
            i for i in g.generators
            if i not in self.sigma and any(g.get_order(i, j) != 2 for j in self.sigma)
        ]

        children = []
        for i in candidate_descents:
            # add child with new descent i
            child = self.copy()
            child._clear_constraints()
            child.sigma[i] = sum([-Polynomial({j: 1}) * CoxeterVector(g, j) for j in g.generators])
            for j in self.sigma:
                child.constraints.add_zero_constraint(
                    child.sigma[i].eval_bilinear(child.sigma[j]) - g.eval_bilinear(i, j)
                )
            child = child.reduce()
            children.append(child)
        return children

    def branch(self):
        def reduce_children(child_list):
            reduced = [child.reduce() for child in child_list]
            return [child for child in reduced if child.is_valid()]

        children = self._get_children_from_determinant_constraint()
        if children is not None:
            return reduce_children(children)

        children = self._get_children_from_unconditional_descent()
        if children is not None:
            return reduce_children(children)

        children = self._get_children_from_conditional_descent()
        if children is not None:
            return reduce_children(children)

        raise Exception('Current state does not match any branching rule: %s' % self)  # pragma: no cover

    def eliminate_descents(self):
        """
        Given constant BraidSystem, successively conjugate by its descents until no descents
        remain, or we can determine that the system is invalid, redundant, or unrealizable.
        """
        logger.debug("Eliminating descents from constant system.")
        history = [self]
        try:
            while True:
                new = history[-1]

                descent = new.get_unconditional_descent()
                if descent is None:
                    return [new]

                commutes = new.sigma[descent] == -CoxeterVector(new.graph, new.graph.star(descent))
                next_new = new._branch_from_descent(descent, commutes=commutes)
                if not next_new.is_valid():
                    return []
                new = next_new

                logger.debug("  descent depth: %s", len(history))
                history.append(new)

                # periodically check whether state is recurrent; 32 is an arbitrary number
                if len(history) % 32 == 0 and new.is_unrealizable(history):
                    return []

        except KeyboardInterrupt:  # pragma: no cover
            logger.warning('Could not compute children for possibly recurrent state:\n%s' % self)  # pragma: no cover
            raise RecurrentStateException(new)  # pragma: no cover

    def get_unconditional_descent(self):
        """
        A generator index i is an unconditional descent if self.sigma[i] is a
        CoxeterVector of the form sum_j f_j alpha_j and for some j it holds that
        f_j is a polynomial with nonpositive coefficients and negative constant term.
        Returns such an index i if one exists, and otherwise None.
        """
        for i in sorted(self.sigma):
            if any(
                self.constraints.is_value_negative(f)
                for f in self.sigma[i].coefficients.values()
            ):
                return i

    def get_conditional_descent(self):
        """
        A generator index i is a conditional descent if self.sigma[i] is a
        CoxeterVector of the form sum_j f_j alpha_j and for some j it holds
        that f_j <= 0 or f_j <= g for some nonpositive constraint g.
        Returns such an index i if one exists, and otherwise None.
        """
        for i in sorted(self.sigma):
            if any(
                self.constraints.is_value_nonpositive(f)
                for f in self.sigma[i].coefficients.values()
            ):
                return i

    def _get_children_from_determinant_constraint(self):
        if not self.sigma.is_complete():
            return None

        determinant = self.sigma.determinant()
        if determinant is None or determinant in [-1, 1]:
            return None

        logger.debug("Constructing children from determinant constraint.")
        positive_child = self.copy()
        negative_child = self.copy()
        try:
            positive_child.constraints.add_zero_constraint(determinant - 1)
            negative_child.constraints.add_zero_constraint(determinant + 1)
        except InvalidInputException:  # pragma: no cover
            return None  # pragma: no cover
        else:
            return [positive_child, negative_child]

    def _get_children_from_unconditional_descent(self):
        descent = self.get_unconditional_descent()
        if descent is None:
            return None

        logger.debug("Constructing children from unconditional descents.")
        children = [
            self._branch_from_descent(descent, commutes=True),
            self._branch_from_descent(descent, commutes=False)
        ]
        return children

    def _get_children_from_conditional_descent(self):
        """Returns list of three children constructed from a conditional descent, or None."""
        descent = self.get_conditional_descent()
        if descent is None:
            return None

        logger.debug("Constructing children from conditional descent.")

        ascent = self.copy()
        ascent.constraints.add_nonpositive_constraint(-self.sigma[descent])

        children = [
            ascent,
            self._branch_from_descent(descent, commutes=True),
            self._branch_from_descent(descent, commutes=False)
        ]
        return children

    def _branch_from_descent(self, i, commutes=True):
        """
        Returns BraidSystem derived from self with respect to generator index i,
        which we define as the BraidSystem given by replacing self.sigma by its
        Demazure conjugate by s_i, extending self.word_s and self.word_t by i on
        the left, and adding the constraint that self.sigma[i] is a negative root.
        If input `commutes` is True, we assume self.sigma s_i = s_i^* self.sigma
        and include the relavent constraint in the derived BraidSystem.
        """
        alpha = CoxeterVector(self.graph, self.graph.star(i))
        beta = self.sigma[i]

        new = self.copy()
        new.constraints.add_nonpositive_constraint(beta)
        new.word_s.extend_left(i)
        new.word_t.extend_left(i)

        if commutes:
            new.constraints.add_zero_constraint(alpha + beta)
            new.sigma = new.sigma * i
        else:
            new.constraints.add_nonzero_constraint(alpha + beta)
            new.sigma = self.graph.star(i) * new.sigma * i

        return new

    def is_unrealizable(self, history):
        """
        Returns True if we can determine automatically that sigma_j has a descent i_j
        for all j (and therefore is never the identity), where we define

            sigma_0 = self.sigma
            sigma_{j+1} = s_{i_j}^* o sigma_j o s_{i_j}.

        We achieve this by trying to interpolate a polynomial formula for sigma_j over
        certain periods of j and then checking this formula by induction. If a BraidSystem
        has this property, then self.sigma cannot represent a Coxeter group element,
        so the BraidSystem is invalid and has no children.
        """
        if not self.sigma.is_constant():
            return False

        logger.debug("Checking recurrence of system:\n%s", self)
        patterns = self._find_patterns(len(history))
        variable = 'x'
        for pattern, repetitions in patterns:
            logger.debug("Found pattern with %i repetitions: %s", repetitions, pattern)
            n = len(pattern) * repetitions
            expected = self._get_generic_sequence(history[-n:], pattern, repetitions, variable)
            if self._confirm_generic_sequence(pattern, expected, variable):
                logger.debug("System is recurrent.")
                return True

        logger.debug("Could not determine if system is recurrent.")
        return False

    def _find_patterns(self, search_length):
        """
        Looks for repeated patterns of the form (i, j, k, ..., i, j, k, ... , i, j, k, ...)
        in self.word_s and self.word_t. Returns list of pairs (sequence, n) where sequence is
        the reversed pattern (..., k, j, i) and n is the number of times it repeatedly occurs.
        """
        patterns = []
        word_s = self.word_s.word
        word_t = self.word_t.word
        n = 2
        while 2 * n <= search_length:
            i = 1
            while True:
                a = word_s[(i - 1) * n:i * n]
                b = word_s[i * n:(i + 1) * n]
                c = word_t[(i - 1) * n:i * n]
                d = word_t[i * n:(i + 1) * n]
                if a == b == c == d and (i + 1) * n <= search_length:
                    i += 1
                else:
                    if i > 1:
                        pattern = reverse_tuple(word_s[:n])
                        patterns += [(pattern, i)]
                    break
            n += 1
        # sort so that largest patterns appear first
        return sorted(patterns, key=lambda x: -len(x[0]))

    def _get_generic_sequence(self, history, pattern, repetitions, variable):
        """
        If the pattern (i, j, ...) has length n and is repeated r times,
        and history is the list of PartialTransforms

            (f_1, f_2, ..., f_nr) = (s_i^* o f_0 o s_i, s_j^* o f_1 o s_j, ... )

        where f_0 = self.sigma, then returns the sequence of (non-constant) PartialTransforms
        whose coefficients are polynomial interpolations of the f_i's. We hope that
        one can then check by induction that this sequence is repeating.
        """
        n = len(pattern)
        matrix = Matrix.vandermonde_inverse(repetitions)

        vectors = [{
            i: [history[j * n + k].sigma[i] for j in range(repetitions)]
            for i in self.sigma
        } for k in range(n)]
        solved = [{
            i: matrix * vectors[k][i]
            for i in self.sigma
        } for k in range(n)]

        x = Polynomial(variable)
        base = Matrix([[x**i for i in range(repetitions)]])
        induction = Matrix([[(x + 1)**i for i in range(repetitions)]])

        return [
            PartialTransform(self.graph, {i: (base * sol[i])[0] for i in self.sigma})
            for sol in solved
        ] + [
            PartialTransform(self.graph, {i: (induction * sol[i])[0] for i in self.sigma})
            for sol in solved
        ]

    def _confirm_generic_sequence(self, pattern, expected, variable):
        """
        Given outputs from _find_pattern and _get_generic_sequence, try to check
        by induction that polynomial interpolation of sigma values describes recurrent
        pattern. Returns True if this succeeds.
        """
        sigma, n = expected[0], len(pattern)
        for index in range(2 * n):
            i = pattern[(index + 1) % n]

            # check that inductive hypothesis holds
            if sigma != expected[index] or not sigma[i].is_negative():
                return False

            # s_i commutes with sigma if and only if the following is zero
            root = sigma[i] + CoxeterVector(self.graph, self.graph.star(i))

            # check that if root is not zero for all values of X, then it is never zero.
            # return False if this fails, as then the next value of sigma depends on X.
            if not any(f < 0 or f > 0 for _, f in root):
                return False

            if root == 0:
                sigma = sigma * i
            else:
                sigma = self.graph.star(i) * sigma * i
        return True

    def is_valid(self):
        """Return True if current state could be parent of final state for a braid relation."""
        if self._are_words_valid() and self._is_sigma_valid() and self.constraints.is_valid():
            return True
        else:
            logger.debug("Invalid or redundant system: %s" % self)
            return False

    def _are_words_valid(self):
        """
        Returns False if word_s and word_t are not reduced, or share a right descent,
        or have (distinct) right descents with product of order greater than m_st.
        These cases may be excluded by induction.
        """
        if not (self.word_s.is_reduced and self.word_t.is_reduced):
            return False
        if not self.word_s.right_descents.isdisjoint(self.word_t.right_descents):
            return False
        if any(self.graph.get_order(self.s, self.t) < self.graph.get_order(u, v)
           for u in self.word_s.right_descents
           for v in self.word_t.right_descents):
                return False
        return True

    def _is_sigma_valid(self):
        """Returns False if sigma sends any root to 0 or vector which is not negative/positive."""
        return not any(not root.is_valid() for root in self.sigma.values())

    def is_redundant(self):
        """
        Returns True if prospective relation encoded by BraidSystem reduces to a shorter
        relation that holds by induction, or is rendered invalid by such a relation.
        """
        if not self.sigma.is_constant() or self.sigma.is_complete():
            return False

        # filter relations to retain only those shorter than current words
        relations = self.sigma.get_relations()
        relations = {(a, b) for a, b in relations if len(a) < len(self.word_s.word)}

        def extract_descents(word):
            """
            Generates equivalence class of group elements from given relations,
            starting from the product s_1s_2...s_n of simple generators corresponding
            to the input word. It is assumed that every relations has the form
            (s, t, s, ...) ~ (t, s, t, ...). Returns union of left descent sets over class.
            """
            descents = set()
            for w in self.graph.get_inverse_atoms(relations, word):
                descents |= w.get_inverse().right_descents
            return descents

        descents_start = extract_descents(self.word_s.word)
        descents_target = extract_descents(self.word_t.word)

        # return True if we can pull out descents s', t' in each word with m(s',t') > m(s,t)
        return any(
            i == j or
            self.graph.get_order(self.s, self.t) < self.graph.get_order(i, j)
            for i in descents_start for j in descents_target
        )

    def is_realized(self):
        """
        Returns True if self.sigma acts trivially, and self.word_s and self._word_t
        are involution words for the same twisted involution.
        """
        if all(self.sigma[i] == CoxeterVector(self.graph, i) for i in self.sigma):
            x = self.word_s.to_involution()
            y = self.word_t.to_involution()
            if x.is_reduced and y.is_reduced and x.left_action == y.left_action:
                return True
        return False

    def reduce(self):
        """Reduce constraints, apply resulting simplifications to self.sigma, and return self."""
        variable_substitutions = self.constraints.simplify()
        for var, substitution in variable_substitutions:
            for i in self.sigma:
                self.sigma[i] = self.sigma[i].set_variable(var, substitution)
        return self


class BraidQueue:

    """
    Class implementing a queue for storing BraidSystem objects.

    The queue is processed by popping the first BraidSystem, computing its children
    and extracting any possible involution braid relations, and then adding all child
    systems to the back of the queue. See methods `next` and `go`.

    Once all processing is done, the class has methods which can be used to compute
    the minimal spanning subset of a sufficient set of involution braid relations,
    and to do some sanity checks.
    """

    def __init__(self, coxeter_graph, s=None, t=None):
        """Inputs `s` and `t` should be elements of `coxeter_graph.generators`."""
        self.graph = coxeter_graph
        self.recurrent_states = []
        self.minimal_relations = []
        self.sufficient_relations = coxeter_graph.get_half_braid_relations()

        # if s or t is not provided, initialize queue with all non-commuting pairs (s, t)
        if s is None or t is None:
            self.queue = [
                BraidSystem(coxeter_graph, u, v, b)
                for u in coxeter_graph.generators for v in coxeter_graph.generators if u < v
                for b in [True, False]
            ]
            self.neighborhood = set(self.graph.generators)
        else:
            assert s in self.graph.generators and t in self.graph.generators
            assert coxeter_graph.get_order(s, t) not in [1, np.infty]
            self.queue = [
                BraidSystem(coxeter_graph, s, t, True),
                BraidSystem(coxeter_graph, s, t, False)
            ]
            self.neighborhood = {s, t, self.graph.star(s), self.graph.star(t)}

        self._filter_initial_queue()

    def _filter_initial_queue(self):
        """
        Removes BraidSystems from self.queue involving commuting generators s and t,
        since we can predict ahead of time that all such systems will become redundant
        after one iteration. Adds relevant half-braid relation when s^* = t.
        """
        filtered = []
        for state in self.queue:
            s, t = tuple(sorted([state.s, state.t]))
            if 2 < self.graph.get_order(s, t) < np.infty:
                filtered += [state]
            elif self.graph.get_order(s, t) == 2 and not state.is_fixer and self.graph.star(s) == t:
                filtered += [state]
        self.queue = filtered

    def __len__(self):
        return len(self.queue)

    def _get_multiplicities(self, state_to_length_fn):
        multiplicities = defaultdict(int)
        for state in self.queue:
            multiplicities[state_to_length_fn(state)] += 1
        if not multiplicities:
            return '0'
        else:
            return ' '.join(['%s^%s' % (ell, mul) for ell, mul in sorted(multiplicities.items())])

    def root_multiplicities(self):
        """Returns string with multiplicities of non-blank roots in sigma field in queue states."""
        return self._get_multiplicities(lambda state: len(state.sigma))

    def word_multiplicities(self):
        """Returns string with multiplicities of lengths of word_s/t fields in queue states."""
        return self._get_multiplicities(lambda state: len(state.word_s))

    def go(self, limit=None, verify=False):
        """
        Process all states in queue to find sufficient spanning relations,
        then reduce these to a minimal set. If integer input `limit` is provided,
        then we only look for relations of length less than or equal to this limit.
        """
        self.find_minimal_relations(limit)
        self.summarize()
        if verify:
            self.verify_relations()

    def find_minimal_relations(self, limit=None):
        """
        Compute spanning set of involution braid relations, then reduces this
        to a minimal subset, which is stored in the `minimal_relation` field.
        """
        self.find_sufficient_relations(limit)
        self.minimize_relations()

    def find_sufficient_relations(self, limit=None):
        logger.info('Step 1: Finding sufficient relations.')
        t0 = time.time()
        while len(self) > 0 and (limit is None or len(self.queue[0]) <= limit):
            self.next()
        t1 = time.time()

        logger.info('')
        logger.info('Duration: %s seconds' % (t1 - t0))

        # sort by word length, then lexicographically
        sufficient = sorted(self.sufficient_relations, key=lambda x: (len(x[0]), x))
        logger.info('')
        logger.info('---------------------')
        logger.info('Sufficient relations:')
        logger.info('---------------------')
        for u, v in sufficient:
            logger.info('%s <---> %s' % (u, v))

    def next(self):
        """Pop first state from queue, compute its children, then append these to the queue."""
        if len(self) == 0:
            logger.info('Queue is empty')
        else:
            next_state = self._get_next_state()
            self._update_neighborhood(next_state)
            self._update_sufficient_relations(next_state)

            try:
                t0 = time.time()
                children = next_state.get_children()
                t1 = time.time()
                logger.info("")
                logger.info("Time: %s seconds" % (t1 - t0))
            except RecurrentStateException as e:
                self.recurrent_states += [e.state]
            else:
                self._update(children)

    def _get_next_state(self):
        next_state = self.queue.pop(0)
        logger.debug('')
        logger.debug('Next state. %s', next_state)
        return next_state

    def _update(self, children):
        """Add new states from input list `children` to self.queue or to self.final."""
        for i, child in enumerate(children):
            logger.debug('%s. %s' % (i + 1, child))
            self._insert(child)

        if len(children) == 0:
            logger.debug('(no new children)')

        logger.info('Systems in queue                 : %s' % len(self))
        logger.info('Multiplicities by word length    : %s' % self.word_multiplicities())
        logger.info('Multiplicities by non-blank roots: %s' % self.root_multiplicities())
        logger.info('Relations found                  : %s' % len(self.sufficient_relations))
        logger.info('Unresolved systems               : %s' % len(self.recurrent_states))

    def _update_neighborhood(self, child):
        self.neighborhood |= set(child.sigma)
        self.neighborhood |= {
            j for j in (set(self.graph.generators) - set(child.sigma))
            if any(self.graph.get_order(i, j) > 2 for i in child.sigma)
        }
        self.neighborhood |= {j for i in child.sigma for j, _ in child.sigma[i]}
        self.neighborhood |= {self.graph.star(j) for j in self.neighborhood}

    def _insert(self, child):
        """Insert child into queue in position preserving ordering by length."""
        i = 0
        while i < len(self.queue) and len(self.queue[i]) < len(child):
            i += 1
        self.queue.insert(i, child)

    def _update_sufficient_relations(self, child):
        """Convert BraidSystem `child` to pair of words and add to self.sufficient_relations."""
        if child.is_realized():
            u = tuple(child.word_s.word)
            v = tuple(child.word_t.word)
            if u < v:
                self.sufficient_relations.add((u, v))
            else:
                self.sufficient_relations.add((v, u))

    def minimize_relations(self):
        """
        Computes minimal subset of self.sufficient_relations which
        generates the same equivalence classes.
        """
        logger.info('')
        logger.info('')
        logger.info('Step 2: Finding minimal sufficient relations.')

        t1 = time.time()
        sufficient = sorted(self.sufficient_relations, key=lambda x: (len(x[0]), x))
        necessary, redundant = self._get_necessary_relations(sufficient)
        rest = [x for x in sufficient if x not in necessary and x not in redundant]
        self.minimal_relations = self._finalize_necessary_relations(necessary, rest)
        t2 = time.time()

        logger.info('')
        logger.info('Duration: %s seconds' % (t2 - t1))

    def summarize(self):
        print('')
        print('-----------------------')
        print('Twisted Coxeter system:')
        print('-----------------------')
        print(self.graph)

        print('')
        print('-----------------------------')
        print('Minimal sufficient relations:')
        print('-----------------------------')
        for u, v in self.minimal_relations:
            print('%s <---> %s' % (u, v))

        if len(self.recurrent_states) > 0:
            print('')
            print('------------------')
            print('Unresolved states:')
            print('------------------')
            print('')
            for i, state in enumerate(self.recurrent_states):
                print('%s. %s' % (i + 1, state))

    def are_atoms_connected(self, relations, start_word, target_word):
        """
        Returns True if `start_word` and `target_word` are connected via input `relations`
        plus the ordinary braid relations. Here, `relations` should be a list of pairs
        of integer tuples indicating consecutive subsequence substitutions which can be
        made at the start of a word.
        """
        target = CoxeterTransform.from_word(self.graph, reverse_tuple(target_word))
        return target in self.graph.get_inverse_atoms(relations, start_word)

    def _get_necessary_relations(self, final):
        """
        Given a sufficent set of relations `final` (consisting of pairs of integer tuples),
        returns pair `(necessary, redundant)` where `necessary` is the list of relations
        in `final` which cannot be generated by all others combined, and the `redundant`
        is the list of relations in `final` which are generated by earlier relations.
        """
        logger.info("")
        logger.info("  1. Get relations which cannot be generated by all others combined.")

        necessary, redundant = [], []
        for i in range(len(final)):
            t0 = time.time()
            tag = '     * Relation %s/%s:' % (i + 1, len(final))

            u, v = final[i]
            initial = [(a, b) for a, b in final[:i] if (a, b) not in redundant]
            complement = [(a, b) for a, b in (final[:i] + final[i + 1:]) if (a, b) not in redundant]

            if self.are_atoms_connected(initial, u, v):
                redundant.append((u, v))
                label = "redundant"
            elif not self.are_atoms_connected(complement, u, v):
                necessary.append((u, v))
                label = "independent"
            else:
                label = "not redundant or independent"

            t1 = time.time()
            logger.info("%s (%s), computation took %s seconds" % (tag, label, t1 - t0))
        return necessary, redundant

    def _finalize_necessary_relations(self, necessary, rest):
        """
        Given list of `necessary` relations as returned by `_get_necessary_relations`
        and the list `rest` of the remaining relations which are not redundant,
        searches through all combinations of relations to return smallest set which spans.
        """
        logger.info("")
        logger.info("  2. Get smallest set of relations which generate all the others.")

        for i in range(len(rest) + 1):
            for current in itertools.combinations(rest, i):
                candidate = necessary + list(current)
                logger.info("     * Trying combination of relations:")
                for rel in candidate:
                    logger.info("         %s <---> %s" % rel)

                complement = [x for x in rest if x not in current]
                if all(self.are_atoms_connected(candidate, u, v) for u, v in complement):
                    logger.info("       Works!")
                    return self._finalize_relations(candidate)
                else:
                    logger.info("       Does not span.")

        raise Exception('Error in BraidQueue._finalize_necessary_relations: returning None')  # pragma: no cover

    def _finalize_relations(self, relations):
        """
        Rewrites relations of the form

            (a_1,...,a_k,s,t,s,...) <---> (a_1,...,a_k,t,s,t,...)

        so that the prefix (a_1,...,a_k) is lexicographically minimal,
        while still spanning the same equivalence classes.
        """
        finalized = []
        for u, v in relations:
            # find largest common initial prefix of u and v
            i = len(u)
            while 0 < i and u[:i] != v[:i]:
                i -= 1
            start_word = u[:i]
            inverse_atoms = self.graph.get_inverse_atoms(relations, start_word)
            w = min(reverse_tuple(a.minimal_reduced_word) for a in inverse_atoms)
            finalized += [(w + u[i:], w + v[i:])]
        return finalized

    @classmethod
    def _get_next_level_of_involutions_to_atoms(cls, graph, current_level=None):
        """
        Returns dictionary mapping involutions of a fixed rank to their sets of atoms.

        The input `current_level` should be a dictionary whose keys are the CoxeterTransforms
        giving all (twisted) involutions of some rank k, and whose values are the sets of
        CoxeterTransforms which are the atoms of each involution.

        The method returns a dictionary of the same type but whose keys are the
        involutions of rank k+1.
        """
        if current_level is None:
            return {CoxeterTransform(graph): {CoxeterTransform(graph)}}

        next_level = defaultdict(set)
        for involution, atoms in current_level.items():
            for i in set(graph.generators) - involution.right_descents:
                next_involution = involution.demazure_conjugate(i)
                next_level[next_involution] |= {w * i for w in atoms if i not in w.right_descents}
        return dict(next_level)

    def verify_relations(self, upper_length=100):
        """
        Checks whether self.minimal_relations span the set of atoms for any twisted involution
        whose involution length is at most upper_length, and prints out the results.
        """
        logger.info('')
        logger.info('')
        logger.info('Step 3: Verifying minimal relations.')
        logger.info('')

        t2 = time.time()
        max_length = self.graph.get_max_involution_word_length(upper_length)

        logger.info('Checking that minimal relations generate the atoms of any twisted involution.')
        logger.info('')

        if max_length == upper_length:
            logger.info('\u26A0 Coxeter group appears to be very large or infinite.')
            logger.info('  Only checking atoms up to length %s.' % max_length)
            logger.info('')
            logger.info('  (If this still takes forever, you can quit with CONTROL+C.)')
            logger.info('')

        self._check_preservation()
        self._check_spanning(max_length)
        t3 = time.time()

        logger.info('')
        logger.info('Verifying relations took %s seconds' % (t3 - t2))

    def _check_preservation(self):
        """Checks that minimal relations are involution words for same twisted involution."""
        logger.info('  * Relations preserve atoms: ', end='')
        g = self.graph
        for a, b in self.minimal_relations:
            y, z = CoxeterWord(g, a).to_involution(), CoxeterWord(g, b).to_involution()
            if y.left_action != z.left_action:
                logger.info('no')
                raise Exception('Error: minimal relations do not preserve all sets of atoms.')
        logger.info('yes')

    def _check_spanning(self, max_length):
        """
        Checks that each atom's equivalence class under minimal relations gives
        all atoms for the corresponding involution.
        """
        g = self.graph
        next_level = self._get_next_level_of_involutions_to_atoms(g)
        for length in range(max_length + 1):
            logger.info('  * Relations span atoms of length %s/%s: ' % (length, max_length), end='')
            for involution, atoms in next_level.items():
                atom = next(iter(atoms))
                word = atom.minimal_reduced_word
                if len(self.graph.get_inverse_atoms(self.minimal_relations, word)) < len(atoms):
                    logger.info('no')
                    raise Exception('Error: minimal relations fail to span all sets of atoms.')
            logger.info('yes')
            next_level = self._get_next_level_of_involutions_to_atoms(g, next_level)
