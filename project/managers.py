from collections import defaultdict
import time
import itertools

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
    InvalidInputException,
    CannotFactorException,
    RecurrentStateException
)


class ConstraintsManager:

    """
    Objects of this class represent a set constraints, given by equalities and inequalities
    involving (mostly linear) polynomials. The class provides methods to simplify such
    constraints, and to determine when they are impossible to simultaneously satisfy.
    """

    def __init__(self):
        # list of linear Polynomials which must be == 0
        self.linear_constraints = set()
        # list of quadratic Polynomials which must be == 0
        self.quadratic_constraints = set()
        # list of linear Polynomials which must be <= 0
        self.nonpositive_constraints = set()
        # list of CoxeterVectors which must be != 0
        self.nonzero_constraints = set()
        self.nonnegative_indeterminates = set()

    def __eq__(self, other):
        return \
            ConstraintsManager == type(other) and \
            self.linear_constraints == other.linear_constraints and \
            self.quadratic_constraints == other.quadratic_constraints and \
            self.nonpositive_constraints == other.nonpositive_constraints and \
            self.nonzero_constraints == other.nonzero_constraints and \
            self.quadratic_constraints == other.quadratic_constraints and \
            self.nonnegative_indeterminates == other.nonnegative_indeterminates

    def copy(self):
        other = ConstraintsManager()
        other.nonpositive_constraints = self.nonpositive_constraints.copy()
        other.linear_constraints = self.linear_constraints.copy()
        other.nonzero_constraints = self.nonzero_constraints.copy()
        other.quadratic_constraints = self.quadratic_constraints.copy()
        other.nonnegative_indeterminates = self.nonnegative_indeterminates.copy()
        return other

    def add_zero_constraint(self, constraint):
        """
        Add input Polynomial `constraint` to linear or quadratic constraints accordingly.
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
        elif degree == 2:
            self.quadratic_constraints.add(constraint)
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
        elif type(constraint) in [int, RationalNumber, QuadraticNumber]:
            constraint = Polynomial(constraint)
        elif type(constraint) != Polynomial:
            raise InvalidInputException(self, constraint, 'add_nonpositive_constraint')

        if -constraint in self.nonpositive_constraints:
            self.add_zero_constraint(constraint)
        elif not any(constraint <= f for f in self.nonpositive_constraints):
            self.nonpositive_constraints.add(constraint)

    def add_nonzero_constraint(self, root):
        if type(root) == CoxeterVector:
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
        self.linear_constraints = {
            f.set_variable(var, substitution) for f in self.linear_constraints
        }
        self.nonpositive_constraints = {
            f.set_variable(var, substitution) for f in self.nonpositive_constraints
        }
        self.nonzero_constraints = {
            r.set_variable(var, substitution) for r in self.nonzero_constraints
        }
        self.quadratic_constraints = {
            f.set_variable(var, substitution) for f in self.quadratic_constraints
        }
        self.nonnegative_indeterminates = {
            f.set_variable(var, substitution) for f in self.nonnegative_indeterminates
        }

    def remove_vacuous_constraints(self):
        self.linear_constraints = {
            f for f in self.linear_constraints if not (f == 0)
        }
        self.quadratic_constraints = {
            f for f in self.quadratic_constraints if not (f == 0)
        }
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
        return not (
            any(0 < f for f in self.nonpositive_constraints) or
            # we do not check f != 0 since this would be true for any nontrivial polynomial;
            # instead, we want to determine if f has all positive/negative coeffs
            any(0 < f or f < 0 for f in self.linear_constraints) or
            any(0 < f or f < 0 for f in self.quadratic_constraints) or
            any(r == 0 for r in self.nonzero_constraints) or
            any(f < 0 for f in self.nonnegative_indeterminates) or
            (len(self.nonnegative_indeterminates) > 0 and sum(self.nonnegative_indeterminates) == 0)
        )


class BraidSystem:
    def __init__(self, coxeter_graph, s, t, is_fixer=True):
        if s in coxeter_graph.generators and t in coxeter_graph.generators:
            self.graph = coxeter_graph
            self.s = s
            self.t = t
            self.is_fixer = is_fixer
            self.sigma = PartialTransform(coxeter_graph)
            self.word_s = CoxeterWord(coxeter_graph)
            self.word_t = CoxeterWord(coxeter_graph)
            self.constraints = ConstraintsManager()
        else:
            raise InvalidInputException(self, (s, t))

    def __eq__(self, other):
        return \
            BraidSystem == type(other) and \
            self.graph == other.graph and \
            self.sigma == other.sigma and \
            self.word_s.left_action == other.word_s.left_action and \
            self.word_t.left_action == other.word_t.left_action and \
            self.constraints == other.constraints

    def __len__(self):
        return len(self.word_s)

    def __repr__(self):
        unconditional = self.get_unconditional_descent()
        conditional, _ = self.get_conditional_descent()
        s = ''
        s += 'BraidSystem data:\n'
        s += '----------------------------------------------------------------\n'
        s += 's = %s, word_s = %s\n' % (self.s, self.word_s)
        s += 't = %s, word_t = %s\n' % (self.t, self.word_t)
        if self.is_fixer:
            s += '\ninitial sigma: alpha_s -> alpha_s^*, alpha_t -> alpha_t^*\n'
        else:
            s += '\ninitial sigma: alpha_s -> alpha_t^*, alpha_t -> alpha_s^*\n'
        s += '\n'
        s += 'sigma = %s' % self.sigma
        s += '\n'
        s += 'unconditional descent: %s\n' % unconditional
        s += 'conditional descent  : %s\n' % conditional
        s += str(self.constraints)
        s += '----------------------------------------------------------------\n'
        return s

    def copy(self):
        other = BraidSystem(self.graph, self.s, self.t)
        other.sigma = self.sigma.copy()
        other.word_s = self.word_s.copy()
        other.word_t = self.word_t.copy()
        other.constraints = self.constraints.copy()
        return other

    def _clear_constraints(self):
        self.constraints = ConstraintsManager()

    def get_children(self):
        if len(self.sigma) == 0:
            return self._get_first_children()

        children = []
        queue = [self.copy().reduce()]
        while queue:
            child = queue.pop(0)
            if child.sigma.is_constant() and child.sigma.is_positive():
                if len(child.word_s) > len(self.word_s) and not child.is_redundant():
                    children.extend(child._get_children_with_new_descent())
            else:
                queue.extend(child.branch())
        return children

    def _get_first_children(self):
        child = BraidSystem(self.graph, self.s, self.t, self.is_fixer)
        child._extend_words()
        if not child._are_words_valid():
            return []

        if self.is_fixer:
            i, j = self.s, self.t
        else:
            i, j = self.t, self.s
        alpha = CoxeterVector(self.graph, self.graph.star(i))
        beta = CoxeterVector(self.graph, self.graph.star(j))
        child.sigma = PartialTransform(self.graph, {self.s: alpha, self.t: beta})
        return child._get_children_with_new_descent()

    def _extend_words(self):
        gens = [self.s, self.t]
        for i in range(self.graph.get_semiorder(self.s, self.t, self.is_fixer)):
            self.word_s.extend_left(gens[i % 2])
            self.word_t.extend_left(gens[(i + 1) % 2])

    def _get_children_with_new_descent(self):
        """
        If self.sigma is not constant, or if self.sigma[i] is defined for all
        generator indices i, or if self.sigma has a descent, returns None.
        Otherwise, returns list of children constructed from current BraidSystem
        by replacing a single undefined value of self.sigma[i] by a generic
        CoxeterVector of the form - sum_i X_i alpha_i, so that i becomes a descent.
        We also include in the returned list one child with no new descents.
        """
        if len(self.sigma) == 0 or not self.sigma.is_constant() or not self.sigma.is_positive():
            raise Exception  # pragma: no cover
        if self.sigma.is_complete():
            return [self]

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
            for j in child.sigma:
                child.constraints.add_zero_constraint(
                    child.sigma[i].eval_bilinear(child.sigma[j]) - g.eval_bilinear(i, j)
                )
            children.append(child)

        # add child with no (new) descents
        child = self.copy()
        for i in g.generators:
            # since child has no descents, it must act as identity element
            if i in self.sigma:
                child.constraints.add_zero_constraint(self.sigma[i] - CoxeterVector(g, i))
            elif i in candidate_descents:
                child.sigma[i] = CoxeterVector(g, i)
        if child.is_valid():
            children.append(child.reduce())

        return children

    def branch(self):
        child_getters = [
            self._get_children_from_determinant_constraint,
            self._get_children_from_quadratic_constraint,
            self._get_children_from_unconditional_descent,
            self._get_children_from_conditional_descent
        ]
        for getter in child_getters:
            childen = getter()
            if childen is not None:
                for child in childen:
                    child = child.reduce()
                    if child.is_valid():
                        yield child
                return
        raise Exception('Current state does not match any branching rule: %s' % self)

    def get_unconditional_descent(self):
        """
        A generator index i is an unconditional descent if self.sigma[i] is a
        CoxeterVector of the form sum_j f_j alpha_j and for some j it holds that
        f_j is a polynomial with nonpositive coefficients and negative constant term.
        Returns such an index i if one exists, and otherwise None.
        """
        descents_to_avoid = self.word_s.left_descents | self.word_t.left_descents
        unconditional = {
            i for i in self.sigma
            if any(f < 0 for f in self.sigma[i].coefficients.values())
        }
        intersection = sorted(unconditional & descents_to_avoid)
        if intersection:
            return intersection[0]
        if unconditional:
            return sorted(unconditional)[0]

    def get_conditional_descent(self):
        """
        A generator index i is a conditional descent if self.sigma[i] is a
        CoxeterVector of the form sum_j f_j alpha_j and for some j it holds that

        (*)    f_j <= 0 or f_j <= g for some nonpositive constraint g.

        If such an index i exists, then returns the pair (i, sum_j f_j alpha_j)
        where the sum is over j such that (*) holds. Otherwise returns (None, None).
        """
        for i in self.sigma:
            nonpositive_part = CoxeterVector(self.graph)
            for j, f in self.sigma[i]:
                if f <= 0 or any(f <= g for g in self.constraints.nonpositive_constraints):
                    nonpositive_part += CoxeterVector(self.graph, j, f)
            if nonpositive_part != 0:
                return i, nonpositive_part
        return None, None

    def _get_children_from_determinant_constraint(self):
        determinant = self.sigma.determinant()
        if determinant is None or determinant in [-1, 1]:
            return None
        positive_child = self.copy()
        negative_child = self.copy()
        try:
            positive_child.constraints.add_zero_constraint(determinant - 1)
            negative_child.constraints.add_zero_constraint(determinant + 1)
        except InvalidInputException:  # pragma: no cover
            return None  # pragma: no cover
        else:
            return [positive_child, negative_child]

    def _get_children_from_quadratic_constraint(self):
        if len(self.constraints.quadratic_constraints) == 0:
            return None
        constraint = next(iter(self.constraints.quadratic_constraints))
        try:
            factors = constraint.get_real_quadratic_factors()
        except CannotFactorException:  # pragma: no cover
            return None  # pragma: no cover
        children = []
        for factor in factors:
            child = self.copy()
            child.constraints.add_zero_constraint(factor)
            children.append(child)
        return children

    def _get_children_from_unconditional_descent(self):
        descent = self.get_unconditional_descent()
        if descent is None:
            return None

        if not self.sigma[descent].is_constant():
            children = [
                self._branch_from_descent(descent, commutes=True),
                self._branch_from_descent(descent, commutes=False)
            ]
            return [c for c in children if c is not None]

        return self._iterate_descents(descent)

    def _iterate_descents(self, descent):
        history = [self]
        try:
            while True:
                new = history[-1]
                commutes = new.sigma[descent] == -CoxeterVector(new.graph, new.graph.star(descent))
                new = new._branch_from_descent(descent, commutes=commutes)
                if new is None:
                    return []

                new = new.reduce()
                if not new.is_valid():
                    return []

                history.append(new)
                print(len(history))

                # periodically whether state is recurrent; 32 is an arbitrary number
                if len(history) % 32 == 0 and new.is_recurrent(history):
                    return []  # pragma: no cover

                descent = new.get_unconditional_descent()
                if descent is None or not new.sigma[descent].is_constant():
                    return [new]
        except KeyboardInterrupt:  # pragma: no cover
            print('\nCould not compute children for possibly recurrent state:\n%s' % self)  # pragma: no cover
            raise RecurrentStateException(new)  # pragma: no cover

    def is_recurrent(self, history):
        """
        Returns True if we can determine automatically that sigma_j has a descent i_j
        for all j (and therefore is never the identity), where we define

            sigma_0 = self.sigma
            sigma_{j+1} = s_{i_j}^* o sigma_j o s_{i_j}.

        We achieve this by trying to interpolate a polynomial formula for sigma_j over
        certain periods of j and then checking this formula by induction. If a BraidSystem
        is recurrent in this sense, then self.sigma cannot represent a Coxeter group element,
        so the BraidSystem is invalid and has no children.
        """
        if not self.sigma.is_constant():
            return False
        print(history[-1])
        patterns = self._find_patterns(len(history))
        variable = 'x'
        for pattern, repetitions in patterns:
            print(pattern, repetitions)
            n = len(pattern) * repetitions
            expected = self._get_generic_sequence(history[-n:], pattern, repetitions, variable)
            if self._confirm_generic_sequence(pattern, expected, variable):
                return True

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
                return False  # pragma: no cover

            # s_i commutes with sigma if and only if the following is zero
            root = sigma[i] + CoxeterVector(self.graph, self.graph.star(i))

            # check that if root is not zero for all values of X, then it is never zero.
            # return False if this fails, as then the next value of sigma depends on X.
            if not any(f < 0 or f > 0 for _, f in root):
                return False  # pragma: no cover

            if root == 0:
                sigma = sigma * i
            else:
                sigma = self.graph.star(i) * sigma * i
        return True

    def _get_children_from_conditional_descent(self):
        """Returns list of three children constructed from a conditional descent, or None."""
        descent, nonpositive_root = self.get_conditional_descent()
        if descent is None:
            return None

        children = []

        # in this child, `descent` is a commuting descent of sigma
        child_a = self._branch_from_descent(descent, commutes=True)
        if child_a is not None:
            child_a.constraints.add_nonzero_constraint(nonpositive_root)
            children.append(child_a)

        # in this child, `descent` is a non-commuting descent of sigma
        child_b = self._branch_from_descent(descent, commutes=False)
        if child_b is not None:
            child_b.constraints.add_nonzero_constraint(nonpositive_root)
            children.append(child_b)

        # in last child, `descent` is not a descent of sigma, so nonpositive_root must be 0
        child_c = self.copy()
        for _, f in nonpositive_root:
            child_c.constraints.add_zero_constraint(f)
        children.append(child_c)

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
        if not new._are_words_valid():
            return None

        if commutes:
            new.constraints.add_zero_constraint(alpha + beta)
            new.sigma = new.sigma * i
        else:
            new.constraints.add_nonzero_constraint(alpha + beta)
            new.sigma = self.graph.star(i) * new.sigma * i
        return new

    def is_valid(self):
        """Return True if current state could be parent of final state for a braid relation."""
        return \
            self._is_sigma_valid() and \
            self.constraints.is_valid()

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

    def is_redundant(self):
        """
        Returns True if prospective relation encoded by BraidSystem reduces to a shorter
        relation that holds by induction, or is rendered invalid by such a relation.
        """
        sigma = self.sigma
        if len(self.word_s.word) > 2 and sigma.is_constant() and not sigma.is_complete():
            if self._words_can_be_reduced(self.word_s.word, self.word_t.word):
                return True
        return False

    def _words_can_be_reduced(self, start_word, target_word):
        """
        Returns True if relation between start_word and target_word reduces to a shorter
        relation that holds by induction, or is rendered invalid by such a relation.
        """
        relations = self.sigma.get_relations()

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

        descents_start = extract_descents(start_word)
        descents_target = extract_descents(target_word)

        # if descents intersect then can pull out common descent to reduce to shorter relation.
        if descents_start & descents_target:
            return True

        # return True if we can pull out descents s', t' in each word with m(s',t') > m(s,t)
        return any(
            self.graph.get_order(self.s, self.t) < self.graph.get_order(i, j)
            for i in descents_start for j in descents_target
        )

    def _is_sigma_valid(self):
        # invalid if sigma sends any root to 0 or non-negative/positive combination of simple roots
        if any(not root.is_valid() for root in self.sigma.values()):
            return False

        # if sigma sends all roots to positive roots, and no variables remain
        if self.sigma.is_constant() and self.sigma.is_complete() and self.sigma.is_positive():
            # invalid if sigma is not trivial
            if not self.sigma.is_identity():
                return False
            # invalid if word_s and word_t are not involution words for the same element
            x = self.word_s.to_involution()
            y = self.word_t.to_involution()
            if not (x.is_reduced and y.is_reduced and x.left_action == y.left_action):
                return False
        return True

    def is_leaf(self):
        return self.sigma.is_constant()

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

    The queue is processed by popping the first BraidSystem, computing its children,
    and then either adding each of these back into the queue or, when a child is a leaf,
    converting it to a involution braid relation. See methods `next` and `go`.

    Once all processing is done, the class has methods which can be used to compute
    the minimal spanning subset of a sufficient set of involution braid relations,
    and to do some sanity checks.
    """

    VERBOSE_LEVEL_NONE = 0
    VERBOSE_LEVEL_LOW = 1
    VERBOSE_LEVEL_MEDIUM = 2
    VERBOSE_LEVEL_HIGH = 3

    def __init__(self, coxeter_graph, s=None, t=None, is_fixer=True,
                 verbose_level=VERBOSE_LEVEL_MEDIUM):
        """Inputs `s` and `t` should be elements of `coxeter_graph.generators`."""
        self.graph = coxeter_graph

        # if s or t is not provided, initialize queue with all pairs of generators (s, t)
        if s is None or t is None:
            self.queue = [
                BraidSystem(coxeter_graph, s, t, b)
                for s in coxeter_graph.generators for t in coxeter_graph.generators if s < t
                for b in [True, False]
            ]
            self.neighborhood = set(self.graph.generators)
        else:
            self.queue = [BraidSystem(coxeter_graph, s, t, is_fixer)]
            self.neighborhood = {s, t}

        self.recurrent_states = []
        self.sufficient_relations = set()
        self.leaf_states = []
        self.minimal_relations = []
        self.verbose_level = verbose_level

    def __len__(self):
        return len(self.queue)

    def _print(self, string, end=None, level=VERBOSE_LEVEL_MEDIUM):
        if level <= self.verbose_level:
            print(string, end=end)  # noqa

    def _print_status(self, string, end=None):
        self._print(string, end=end, level=self.VERBOSE_LEVEL_LOW)

    def _print_verbose(self, string, end=None):
        self._print(string, end=end, level=self.VERBOSE_LEVEL_HIGH)

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
        self.find_relations(limit)
        self.minimize_relations()
        if verify:
            self.verify_relations()

    def find_relations(self, limit=None):
        self._print_status('Step 1: Finding sufficient relations.')
        t0 = time.time()
        while len(self) > 0 and (limit is None or len(self.queue[0]) <= limit):
            self.next()
        t1 = time.time()

        self._print_status('')
        self._print_status('Duration: %s seconds' % (t1 - t0))
        self._print_status('')
        self._print_status('-----------------------')
        self._print_status('Twisted Coxeter system:')
        self._print_status('-----------------------')
        self._print_status(self.graph)

        # sort by word length, then lexicographically
        sufficient = sorted(self.sufficient_relations, key=lambda x: (len(x[0]), x))
        self._print_status('')
        self._print_status('---------------------')
        self._print_status('Sufficient relations:')
        self._print_status('---------------------')
        for u, v in sufficient:
            self._print_status('%s <---> %s' % (u, v))

    def next(self):
        """Pop first state from queue, compute its children, then append these to the queue."""
        if len(self) == 0:
            self._print('Queue is empty')
        else:
            next_state = self._get_next_state()
            try:
                t0 = time.time()
                children = next_state.get_children()
                t1 = time.time()
                description = "\nTime: %s seconds" % (t1 - t0)
            except RecurrentStateException as e:
                self.recurrent_states += [e.state]
            else:
                self._update(children, description)

    def _get_next_state(self):
        next_state = self.queue.pop(0)
        self._print_verbose('')
        self._print_verbose('-----------')
        self._print_verbose('Next state:')
        self._print_verbose('-----------')
        self._print_verbose('')
        self._print_verbose(next_state)
        return next_state

    def _update(self, children, description):
        """Add new states from input list `children` to self.queue or to self.final."""
        self._print(description)
        self._print_verbose('')
        self._print_verbose('-------------')
        self._print_verbose('Child states:')
        self._print_verbose('-------------')
        self._print_verbose('')

        i = 0
        for child in children:
            self.neighborhood |= set(child.sigma)
            if child.is_leaf():
                self._add_to_sufficient_relations(child)
            else:
                self._insert(child)
                self._print_verbose('%s. %s' % (i + 1, child))
                i += 1
        if i == 0:
            self._print_verbose('\n(no states added)\n')

        self._print('States in queue                  : %s' % len(self))
        self._print('Multiplicities by word length    : %s' % self.word_multiplicities())
        self._print('Multiplicities by non-blank roots: %s' % self.root_multiplicities())
        self._print('Relations found                  : %s' % len(self.sufficient_relations))
        self._print('Unresolved states                : %s' % len(self.recurrent_states))

    def _insert(self, child):
        """Insert child into queue in position preserving ordering by length."""
        i = 0
        while i < len(self.queue) and len(self.queue[i]) < len(child):
            i += 1
        self.queue.insert(i, child)

    def _add_to_sufficient_relations(self, child):
        """Convert BraidSystem `child` to pair of words and add to self.sufficient_relations."""
        u = tuple(child.word_s.word)
        v = tuple(child.word_t.word)
        if u < v:
            self.sufficient_relations.add((u, v))
        else:
            self.sufficient_relations.add((v, u))
        self.leaf_states.append(child)

    def minimize_relations(self):
        """
        Computes minimal subset of self.sufficient_relations which
        generates the same equivalence classes.
        """
        self._print_status('')
        self._print_status('')
        self._print_status('Step 2: Finding minimal sufficient relations.')

        t1 = time.time()
        sufficient = sorted(self.sufficient_relations, key=lambda x: (len(x[0]), x))
        necessary, redundant = self._get_necessary_relations(sufficient)
        rest = [x for x in sufficient if x not in necessary and x not in redundant]
        self.minimal_relations = self._finalize_necessary_relations(necessary, rest)
        t2 = time.time()

        self._print_status('')
        self._print_status('Duration: %s seconds' % (t2 - t1))
        self._print_status('')
        self._print_status('-----------------------------')
        self._print_status('Minimal sufficient relations:')
        self._print_status('-----------------------------')
        for u, v in self.minimal_relations:
            self._print_status('%s <---> %s' % (u, v))

        if len(self.recurrent_states) > 0:
            self._print_status('')
            self._print_status('------------------')
            self._print_status('Unresolved states:')
            self._print_status('------------------')
            self._print_status('')
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
        self._print("")
        self._print("  1. Get relations which cannot be generated by all others combined.")

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
            self._print("%s (%s), computation took %s seconds" % (tag, label, t1 - t0))
        return necessary, redundant

    def _finalize_necessary_relations(self, necessary, rest):
        """
        Given list of `necessary` relations as returned by `_get_necessary_relations`
        and the list `rest` of the remaining relations which are not redundant,
        searches through all combinations of relations to return smallest set which spans.
        """
        self._print("\n  2. Get smallest set of relations which generate all the others.")

        for i in range(len(rest) + 1):
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
        self._print_status('')
        self._print_status('')
        self._print_status('Step 3: Verifying minimal relations.')
        self._print_status('')

        t2 = time.time()
        max_length = self.graph.get_max_involution_word_length(upper_length)

        self._print('Checking that minimal relations generate the atoms of any twisted involution.')
        self._print('')

        if max_length == upper_length:
            self._print('\u26A0 Coxeter group appears to be very large or infinite.')
            self._print('  Only checking atoms up to length %s.' % max_length)
            self._print('')
            self._print('  (If this still takes forever, you can quit with CONTROL+C.)')
            self._print('')

        self._check_preservation()
        self._check_spanning(max_length)
        t3 = time.time()

        self._print('')
        self._print_status('Verifying relations took %s seconds' % (t3 - t2))

    def _check_preservation(self):
        """Checks that minimal relations are involution words for same twisted involution."""
        self._print('  * Relations preserve atoms: ', end='')
        g = self.graph
        for a, b in self.minimal_relations:
            y, z = CoxeterWord(g, a).to_involution(), CoxeterWord(g, b).to_involution()
            if y.left_action != z.left_action:
                self._print('no')
                raise Exception('Error: minimal relations do not preserve all sets of atoms.')
        self._print('yes')

    def _check_spanning(self, max_length):
        """
        Checks that each atom's equivalence class under minimal relations gives
        all atoms for the corresponding involution.
        """
        g = self.graph
        next_level = self._get_next_level_of_involutions_to_atoms(g)
        for length in range(max_length + 1):
            self._print('  * Relations span atoms of length %s/%s: ' % (length, max_length), end='')
            for involution, atoms in next_level.items():
                atom = next(iter(atoms))
                word = atom.minimal_reduced_word
                if len(self.graph.get_inverse_atoms(self.minimal_relations, word)) < len(atoms):
                    self._print('no')
                    raise Exception('Error: minimal relations fail to span all sets of atoms.')
            self._print('yes')
            next_level = self._get_next_level_of_involutions_to_atoms(g, next_level)
