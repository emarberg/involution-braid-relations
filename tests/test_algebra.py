import pytest

from ..braids import RationalNumber, QuadraticNumber, Polynomial, PrimeFactorization, Monomial


class TestRationalNumbers:

    @pytest.mark.parametrize("p, q, numerator, denominator", [
        (1, 2, 1, 2),
        (-1, 2, -1, 2),
        (1, -2, -1, 2),
        (48, 96, 1, 2),
        (-48, -96, 1, 2),
        (RationalNumber(1, 2), 3, 1, 6),
        (RationalNumber(2, 3), RationalNumber(4, 7), 7, 6),
        (-1, RationalNumber(15, 6), -2, 5),
        (-4 * 9 * 5 * 7, -2 * 3 * 25 * 49, 2 * 3, 5 * 7)
    ])
    def test_constructor(self, p, q, numerator, denominator):
        """Test constructor for RationalNumber with a variety of valid inputs."""
        x = RationalNumber(p, q)
        assert x.numerator == numerator and x.denominator == denominator

    def test_repr(self):
        """Tests for conversion of RationalNumber to str."""
        assert str(RationalNumber()) == '0'
        assert str(RationalNumber(1)) == '1'
        assert str(RationalNumber(-5, 7)) == '-5/7'

    def test_hash(self):
        """Test that hashes for RationalNumbers which are integers match the hashes of the ints."""
        N = 100
        for p in range(-N, N):
            for q in range(-N, N):
                try:
                    x = RationalNumber(p, q)
                except:
                    assert q == 0
                else:
                    if p % q == 0:
                        assert hash(x) == hash(p//q)
                        assert x.is_integer()

        assert len({1, RationalNumber(1, 1)}) == 1

    @pytest.mark.parametrize("a, b, expected", [
    ])
    def test_le(self, a, b, expected):
        """Test <= operator for RationalNumbers."""
        assert (a <= b) == (-b <= -a) == expected

    @pytest.mark.parametrize("a, b, expected", [
    ])
    def test_lt(self, a, b, expected):
        """Test <= operator for RationalNumbers."""
        assert (a < b) == (-b < -a) == expected

    @pytest.mark.parametrize("a, b, c", [
        (RationalNumber(3, 7), 0, RationalNumber(3, 7)),
        (RationalNumber(3, 7), 1, RationalNumber(10, 7)),
        (RationalNumber(3, 7), -2, RationalNumber(-11, 7)),
        (RationalNumber(3, 7), RationalNumber(-3, 7), 0),
        (RationalNumber(3, 7), RationalNumber(4, 9), RationalNumber(3*9 + 4*7, 7*9)),
    ])
    def test_addition(self, a, b, c):
        """Tests for addition and substraction of RationalNumbers."""
        assert a + b == b + a == c
        assert -a - b == -b - a == -c

    @pytest.mark.parametrize("a, b, c", [
        (RationalNumber(3, 7), 0, 0),
        (RationalNumber(3, 7), 1, RationalNumber(3, 7)),
        (RationalNumber(3, 7), -2, RationalNumber(-6, 7)),
        (RationalNumber(3, 7), RationalNumber(-3, 7), RationalNumber(-9, 49)),
        (RationalNumber(3, 7), RationalNumber(4, 9), RationalNumber(4, 21)),
    ])
    def test_multiplication(self, a, b, c):
        """Tests for multiplication of RationalNumbers."""
        assert a * b == b * a == c

    @pytest.mark.parametrize("a, b, c", [
        (RationalNumber(3, 7), 1, RationalNumber(3, 7)),
        (1, RationalNumber(3, 7), RationalNumber(7, 3)),
        (RationalNumber(3, 7), -2, RationalNumber(-3, 14)),
        (-2, RationalNumber(3, 7), RationalNumber(-14, 3)),
        (RationalNumber(3, 7), RationalNumber(-3, 7), -1),
        (RationalNumber(3, 7), RationalNumber(4, 9), RationalNumber(27, 28)),
    ])
    def test_division(self, a, b, c):
        """Tests for  division of RationalNumbers."""
        assert a / b == c

    @pytest.mark.parametrize("p, q, e, num, den", [
        (3, 2, 4, 81, 16),
        (3, 2, -4, 16, 81),
        (3, 2, 0, 1, 1)
    ])
    def test_power(self, p, q, e, num, den):
        """Tests for exponentiation of RationalNumbers."""
        q = RationalNumber(p, q)**e
        assert q.numerator == num and q.denominator == den

    @pytest.mark.parametrize("p, q", [
        (1, 0),
        (1.0, 20),
        (1, 20.0),
        (0, 0),
        (QuadraticNumber(), 1),
        (Polynomial(), 1),
        (RationalNumber(1, 2), None),
        (None, RationalNumber(1, 2)),
        (1, RationalNumber(0, 1)),
        (RationalNumber(7, 8), RationalNumber(0, -100))
    ])
    def test_constructor_errors(self, p, q):
        """Test error handling of invalid arguments passed to RationalNumber constructor."""
        try:
            # denominator must be nonzero
            RationalNumber(p, q)
        except Exception as e:
            assert str(e).startswith('Invalid input types to RationalNumber: ')
        else:
            assert False

    @pytest.mark.parametrize("a, b", [
        (RationalNumber(3, 7), 1.2),
        (1.2, RationalNumber(3, 7)),
        (RationalNumber(3, 7), None),
        (None, RationalNumber(3, 7)),
    ])
    def test_operator_errors(self, a, b):
        """Tests for error handling of invalid operations involving RationalNumbers."""
        try:
            a + b
        except:
            pass
        else:
            assert False

        try:
            a - b
        except:
            pass
        else:
            assert False

        try:
            a * b
        except:
            pass
        else:
            assert False

        try:
            a / b
        except:
            pass
        else:
            assert False

    @pytest.mark.parametrize("zero", [
        0,
        0.0,
        RationalNumber(0, 1),
        QuadraticNumber(0),
        Polynomial(),
    ])
    def test_division_errors(self, zero):
        """Test error handling for attempted division of RationalNumber by zero."""
        a = RationalNumber(3, 7)
        try:
            a / zero
        except:
            pass
        else:
            assert False

    def test_power_errors(self):
        # cannot compute 0**0
        try:
            RationalNumber(0, 1)**0
        except:
            pass
        else:
            assert False

        # cannot compute RationlNumber**RationalNumber, even if exponent is integer
        try:
            RationalNumber(5, 2)**RationalNumber(6, 2)
        except:
            pass
        else:
            assert False


class TestPrimeFactorization:

    @pytest.mark.parametrize("i, expected", [
        (1, {}),
        (-1, {-1: 1}),
        (2*3*5*7, {2: 1, 3: 1, 5: 1, 7: 1}),
        (8*27*125, {2: 3, 3: 3, 5: 3}),
        (-3, {-1: 1, 3: 1}),
        (-2*3*5*7, {-1: 1, 2: 1, 3: 1, 5: 1, 7: 1}),
        (-8*27*125, {-1: 1, 2: 3, 3: 3, 5: 3}),
    ])
    def test_constructor(self, i, expected):
        pf = PrimeFactorization(i)
        assert pf.factorization == expected

    def test_repr(self):
        """Test conversion of PrimeFactorization to str."""
        assert str(PrimeFactorization(2*3*5*7)) == '{2: 1, 3: 1, 5: 1, 7: 1}'

    def test_lt(self):
        """Test < operator for PrimeFactorizations."""
        assert PrimeFactorization(-2*3) < PrimeFactorization(-2) < PrimeFactorization(-1)
        assert PrimeFactorization(-1) < PrimeFactorization(1) < PrimeFactorization(2*3)

    @pytest.mark.parametrize("m, n, expected", [
        (1, 1, {}),
        (1, -1, {-1: 1}),
        (-1, -1, {}),
        (2*3*5, -7*9*11, {-1: 1, 2: 1, 3: 3, 5: 1, 7: 1, 11: 1}),
        (-2*3*5, -7*9*11, {2: 1, 3: 3, 5: 1, 7: 1, 11: 1})
    ])
    def test_multiplication(self, m, n, expected):
        """Tests for multiplication of PrimeFactorizations."""
        m = PrimeFactorization(m)
        n = PrimeFactorization(n)
        assert (n * m).factorization == (m * n).factorization == expected

    @pytest.mark.parametrize("i, expected", [
        (1, {}),
        (-1, {-1: 1}),
        (2, {2: 1}),
        (-4, {-1: 1}),
        (4*3*25, {3: 1})
    ])
    def test_get_square_free_part(self, i, expected):
        pf = PrimeFactorization(i)
        assert pf.get_square_free_part().factorization == expected

    @pytest.mark.parametrize("i, expected", [
        (1, 1),
        (-1, 1),
        (2, 1),
        (-4, 2),
        (4*3*125, 10)
    ])
    def get_truncated_square_root(self, i, expected):
        pf = PrimeFactorization(i)
        assert pf.get_truncated_square_root().factorization == expected

    @pytest.mark.parametrize("i", [
        0,
        0.0,
        QuadraticNumber(),
        RationalNumber()
    ])
    def test_constructor_errors(self, i):
        """Test error handling of invalid arguments passed to QuadraticNumber constructor."""
        try:
            PrimeFactorization(i)
        except Exception as e:
            assert str(e).startswith('Invalid input to PrimeFactorization: ')
        else:
            assert False

    @pytest.mark.parametrize("i", [
        12,
        3.0,
        QuadraticNumber(),
        RationalNumber()
    ])
    def test_multiplication_errors(self, i):
        """Test error handling for invalid multiplication of PrimeFactorizations."""
        try:
            PrimeFactorization(12) * i
        except Exception as e:
            assert str(e).startswith('Cannot multiply PrimeFactorization with ')
        else:
            assert False


class TestQuadraticNumbers:

    @pytest.mark.parametrize("i, expected", [
        (0, {}),
        (1, {PrimeFactorization(1): 1}),
        (-1, {PrimeFactorization(1): -1}),
        (RationalNumber(7, 8), {PrimeFactorization(1): RationalNumber(7, 8)}),
        (RationalNumber(0), {})
    ])
    def test_constructor(self, i, expected):
        """Test constructor for QuadraticNumber with a variety of valid inputs."""
        n = QuadraticNumber(i)
        assert n.coefficients == expected

    def test_repr(self):
        """Tests for conversion of QuadraticNumber to str."""
        assert str(QuadraticNumber()) == '0'
        assert str(QuadraticNumber(1)) == '1'
        assert str(QuadraticNumber(RationalNumber(-5, 7))) == '-5/7'
        assert str(QuadraticNumber.sqrt((RationalNumber(-25, 3)))) == '(5/3*sqrt(-3))'
        assert str(-QuadraticNumber.sqrt((RationalNumber(-3)))) == '(-sqrt(-3))'
        assert str(-3 * QuadraticNumber.sqrt((RationalNumber(-25, 3)))) == '(-5*sqrt(-3))'

    IRRATIONAL_SQUARE_1 = 3 + QuadraticNumber.sqrt(5)
    IRRATIONAL_SQUARE_2 = 7 + 3*QuadraticNumber.sqrt(5)

    @pytest.mark.parametrize("i, expected", [
        (0, {}),
        (RationalNumber(), {}),
        (QuadraticNumber(), {}),
        (1, {1: 1}),
        (-4, {-1: 2}),
        (8*27, {6: 6}),
        (125, {5: 5}),
        (RationalNumber(7, 8), {14: RationalNumber(1, 4)}),
        (RationalNumber(-4, 3), {-3: RationalNumber(2, 3)}),
        (RationalNumber(125), {5: 5}),
        (QuadraticNumber(125), {5: 5}),
        (QuadraticNumber(RationalNumber(125)), {5: 5}),
        (QuadraticNumber(RationalNumber(7, 8)), {14: RationalNumber(1, 4)}),
        (QuadraticNumber(RationalNumber(-4, 3)), {-3: RationalNumber(2, 3)}),
        (IRRATIONAL_SQUARE_1, {2: RationalNumber(1, 2), 10: RationalNumber(1, 2)}),
        (5 * IRRATIONAL_SQUARE_1 / 7, {70: RationalNumber(1, 14), 14: RationalNumber(5, 14)}),
        (IRRATIONAL_SQUARE_2, {2: RationalNumber(3, 2), 10: RationalNumber(1, 2)}),
        (RationalNumber(-4, 3) * IRRATIONAL_SQUARE_2, {-6: 1, -30: RationalNumber(1, 3)}),
    ])
    def test_sqrt(self, i, expected):
        """Test QuadraticNumber.sqrt with a variety of valid inputs."""
        n = QuadraticNumber.sqrt(i)
        assert {f.n: v for f, v in n.coefficients.items()} == expected

    def test_hash(self):
        """Test that hashes for QuadraticNumbers are consistent."""
        hash(QuadraticNumber(0)) == hash(0)
        hash(QuadraticNumber(-1)) == hash(-1)
        hash(QuadraticNumber(12)) == hash(12)
        hash(QuadraticNumber(2000)) == hash(2000)

        hash(QuadraticNumber(RationalNumber(3, 7))) == hash(RationalNumber(3, 7))
        hash(QuadraticNumber(RationalNumber(-3, 7))) == hash(RationalNumber(-3, 7))
        hash(QuadraticNumber(RationalNumber(3, 700))) == hash(RationalNumber(3, 700))
        hash(QuadraticNumber(RationalNumber(700))) == hash(700)

        hash(QuadraticNumber.sqrt(2))

    @pytest.mark.parametrize("a, b, expected", [
    ])
    def test_le(self, a, b, expected):
        """Test <= operator for QuadraticNumbers."""
        assert (a <= b) == (-b <= -a) == expected

    @pytest.mark.parametrize("a, b, expected", [
    ])
    def test_lt(self, a, b, expected):
        """Test <= operator for QuadraticNumbers."""
        assert (a < b) == (-b < -a) == expected

    @pytest.mark.parametrize("a, b, expected", [
        (QuadraticNumber(0), 0, {}),
        (QuadraticNumber(0), RationalNumber(0), {}),
        (QuadraticNumber(0), QuadraticNumber(0), {}),
        (QuadraticNumber(5), -5, {}),
        (QuadraticNumber(5), RationalNumber(-5), {}),
        (QuadraticNumber(5), QuadraticNumber(-5), {}),
        (QuadraticNumber(-2), QuadraticNumber.sqrt(5), {1: -2, 5: 1}),
        (QuadraticNumber(-2), QuadraticNumber.sqrt(-5), {1: -2, -5: 1}),
        (QuadraticNumber.sqrt(125), QuadraticNumber.sqrt(5), {5: 6}),
        (QuadraticNumber.sqrt(125), QuadraticNumber.sqrt(2), {2: 1, 5: 5}),
    ])
    def test_addition(self, a, b, expected):
        """Tests for addition and substraction of QuadraticNumbers."""
        assert {f.n: v for f, v in (a + b).coefficients.items()} == expected
        assert {f.n: v for f, v in (b + a).coefficients.items()} == expected
        assert {f.n: v for f, v in (-(-a - b)).coefficients.items()} == expected
        assert {f.n: v for f, v in (-(-b - a)).coefficients.items()} == expected

    @pytest.mark.parametrize("a, b, expected", [
        (QuadraticNumber(0), 0, {}),
        (QuadraticNumber(0), RationalNumber(0), {}),
        (QuadraticNumber(0), QuadraticNumber(0), {}),
        (QuadraticNumber(1), 1, {1: 1}),
        (QuadraticNumber(1), RationalNumber(1), {1: 1}),
        (QuadraticNumber(1), QuadraticNumber(1), {1: 1}),
        (QuadraticNumber(-1), 1, {1: -1}),
        (QuadraticNumber(-1), RationalNumber(1), {1: -1}),
        (QuadraticNumber(-1), QuadraticNumber(1), {1: -1}),
        (QuadraticNumber(-1), -1, {1: 1}),
        (QuadraticNumber(-1), RationalNumber(-1), {1: 1}),
        (QuadraticNumber(-1), QuadraticNumber(-1), {1: 1}),
        (QuadraticNumber(5), QuadraticNumber.sqrt(5), {5: 5}),
        (QuadraticNumber.sqrt(5), QuadraticNumber.sqrt(5), {1: 5}),
        (QuadraticNumber.sqrt(125), QuadraticNumber.sqrt(3), {15: 5}),
    ])
    def test_multiplication(self, a, b, expected):
        """Tests for multiplication and division of QuadraticNumbers."""
        assert {f.n: v for f, v in (a * b).coefficients.items()} == expected
        assert {f.n: v for f, v in (b * a).coefficients.items()} == expected

    @pytest.mark.parametrize("a, b, expected", [
        (QuadraticNumber(0), 1, {}),
        (QuadraticNumber(1), 12, {1: RationalNumber(1, 12)}),
        (QuadraticNumber(1), RationalNumber(1, 12), {1: 12}),
        (QuadraticNumber(1), QuadraticNumber(1), {1: 1}),
        (QuadraticNumber(-1), 1, {1: -1}),
        (QuadraticNumber(-1), RationalNumber(1), {1: -1}),
        (QuadraticNumber(-1), QuadraticNumber(1), {1: -1}),
        (QuadraticNumber(-1), -1, {1: 1}),
        (QuadraticNumber(-1), RationalNumber(-1), {1: 1}),
        (QuadraticNumber(-1), QuadraticNumber(-1), {1: 1}),
        (QuadraticNumber(5), QuadraticNumber.sqrt(5), {5: 1}),
        (QuadraticNumber.sqrt(5), QuadraticNumber.sqrt(5), {1: 1}),
        (QuadraticNumber.sqrt(125), QuadraticNumber.sqrt(3), {15: RationalNumber(5, 3)}),
    ])
    def test_division(self, a, b, expected):
        """Tests for multiplication and division of QuadraticNumbers."""
        assert {f.n: v for f, v in (a / b).coefficients.items()} == expected

    @pytest.mark.parametrize("n, e, expected", [
        (QuadraticNumber.sqrt(12), 0, {1: 1}),
        (QuadraticNumber(0), 3, {}),
        (QuadraticNumber(1), 10, {1: 1}),
        (QuadraticNumber.sqrt(2), -3, {2: RationalNumber(1, 4)}),
        ((1 - QuadraticNumber.sqrt(3))/2, 2, {1: 1, 3: -RationalNumber(1, 2)}),
        ((1 - QuadraticNumber.sqrt(3))/2, -2, {1: 4, 3: 2}),
    ])
    def test_power(self, n, e, expected):
        """Tests for exponentiation of QuadraticNumbers."""
        assert {f.n: v for f, v in (n**e).coefficients.items()} == expected

    @pytest.mark.parametrize("i", [
        0.0,
        QuadraticNumber(),
        Polynomial()
    ])
    def test_constructor_errors(self, i):
        """Test error handling of invalid arguments passed to QuadraticNumber constructor."""
        try:
            QuadraticNumber(i)
        except Exception as e:
            assert str(e).startswith('Invalid input type to QuadraticNumber: ')
        else:
            assert False

    @pytest.mark.parametrize("i", [
        (0.0),
        (QuadraticNumber.sqrt(2)),
        (Polynomial(4))
    ])
    def test_sqrt_errors(self, i):
        """Test error handling of invalid arguments passed to QuadraticNumber constructor."""
        try:
            QuadraticNumber.sqrt(i)
        except Exception as e:
            assert str(e).startswith('Cannot compute square root of ')
        else:
            assert False

    @pytest.mark.parametrize("zero", [
        0,
        0.0,
        RationalNumber(0, 1),
        QuadraticNumber(0),
        Polynomial(),
    ])
    def test_division_errors(self, zero):
        """Test error handling for attempted division of QuadraticNumber by zero."""
        a = QuadraticNumber(RationalNumber(3, 7)) + QuadraticNumber.sqrt(2)
        try:
            a / zero
        except:
            pass
        else:
            assert False

    @pytest.mark.parametrize("i", [
        0.0,
        1.0,
        QuadraticNumber(),
        RationalNumber(),
        QuadraticNumber(1),
        RationalNumber(1)
    ])
    def test_power_errors(self, i):
        """Test error handling for invalid exponentiation of Monomials."""
        try:
            QuadraticNumber(12) ** i
        except Exception as e:
            assert str(e).startswith(
                '** not implemented when exponent is non-positive or non-integer'
            )
        else:
            assert False

    def test_indeterminate_power(self):
        try:
            QuadraticNumber(0)**0
        except Exception as e:
            assert str(e).startswith('Cannot compute indeterminate power 0**0')
        else:
            assert False

    @pytest.mark.parametrize("a, b", [
        (QuadraticNumber(3), 1.2),
        (1.2, QuadraticNumber(3)),
        (QuadraticNumber(7), None),
        (None, QuadraticNumber(7)),
    ])
    def test_operator_errors(self, a, b):
        """Tests for error handling of invalid operations involving QuadraticNumbers."""
        try:
            a + b
        except:
            pass
        else:
            assert False

        try:
            a - b
        except:
            pass
        else:
            assert False

        try:
            a * b
        except:
            pass
        else:
            assert False

        try:
            a / b
        except:
            pass
        else:
            assert False


class TestMonomial:

    @pytest.mark.parametrize("i, expected", [
        (None, {}),
        ('x', {ord('x'): 1}),
        (12, {12: 1}),
        ({1: 1, 3: 1}, {1: 1, 3: 1}),
        ({1: 1, 3: 1, 4: 0}, {1: 1, 3: 1}),
    ])
    def test_constructor(self, i, expected):
        """Test constructor for Monomial with a variety of valid inputs."""
        assert Monomial(i).exponents == expected

    def test_repr(self):
        """Test conversion of Monomial to str."""
        assert str(Monomial()) == '1'
        assert str(Monomial('x')) == 'x'
        assert str(Monomial(-1)) == 'y_1'
        assert str(Monomial({0: 1, 1: 2, 2: 3})) == 'x_0x_1^2x_2^3'

    def test_eq(self):
        """Test == operator for Monomials."""
        assert Monomial('x') == Monomial(ord('x'))
        assert Monomial() == 1 == RationalNumber(1) == QuadraticNumber(1) == Polynomial(1)

    def test_lt(self):
        """Test < operator for Monomials."""
        assert Monomial() < Monomial(0) < Monomial({0: 1, 1: 1}) < Monomial({0: 2})

    def test_getitem(self):
        """Test [] operator for Monomials."""
        x = Monomial({0: 3})
        assert x[0] == 3
        assert x[1] == 0

    @pytest.mark.parametrize("m, n, expected", [
        ('x', 'y', {ord('x'): 1, ord('y'): 1}),
        ('x', 'x', {ord('x'): 2}),
        ('x', None, {ord('x'): 1}),
        ({0: 1, 1: 2, 2: 3}, {1: 1, 2: 2, 3: 3}, {0: 1, 1: 3, 2: 5, 3: 3}),
        ({0: 1, 1: 2}, {0: -1, 1: -1}, {1: 1})
    ])
    def test_multiplication(self, m, n, expected):
        """Tests for multiplication of Monomials."""
        m = Monomial(m)
        n = Monomial(n)
        assert (n * m).exponents == (m * n).exponents == expected

    @pytest.mark.parametrize("n, e, expected", [
        ('x', 0, {}),
        ('x', 5, {ord('x'): 5}),
        ({0: 1, 1: 2}, 3, {0: 3, 1: 6}),
        ({0: 1, 1: 2}, -3, {0: -3, 1: -6}),
        (None, 12, {}),
    ])
    def test_power(self, n, e, expected):
        n = Monomial(n)
        assert (n**e).exponents == expected

    @pytest.mark.parametrize("i", [
        {0: 1, 'x': 2},
        0.0,
        'abc',
    ])
    def test_constructor_errors(self, i):
        """Test error handling of invalid arguments passed to Monomial constructor."""
        try:
            Monomial(i)
        except Exception as e:
            assert str(e).startswith('Invalid input to Monomial: ')
        else:
            assert False

    @pytest.mark.parametrize("i", [
        1,
        3.0,
        QuadraticNumber(),
        RationalNumber(),
        Polynomial()
    ])
    def test_multiplication_errors(self, i):
        """Test error handling for invalid multiplication of Monomials."""
        try:
            Monomial(12) * i
        except Exception as e:
            assert str(e).startswith('Cannot multiply Monomial with ')
        else:
            assert False

    @pytest.mark.parametrize("i", [
        0.0,
        1.0,
        QuadraticNumber(),
        RationalNumber(),
        QuadraticNumber(1),
        RationalNumber(1)
    ])
    def test_power_errors(self, i):
        """Test error handling for invalid exponentiation of Monomials."""
        try:
            Monomial('x') ** i
        except Exception as e:
            assert str(e).startswith('Cannot exponentiate Monomial by ')
        else:
            assert False


class TestPolynomial:

    @pytest.mark.parametrize("i, expected", [
        (None, {}),
        (0, {}),
        (RationalNumber(0), {}),
        (QuadraticNumber(0), {}),
        ('x', {Monomial('x'): 1}),
        ({0: 1, 1: 2}, {Monomial({0: 1, 1: 2}): 1}),
        (Monomial('x'), {Monomial('x'): 1}),
        (-3, {Monomial(): -3}),
        (RationalNumber(3, 7), {Monomial(): RationalNumber(3, 7)}),
        (QuadraticNumber.sqrt(3)/7, {Monomial(): QuadraticNumber.sqrt(3)/7}),
    ])
    def test_constructor(self, i, expected):
        """Test constructor for Polynomial with a variety of valid inputs."""
        assert Polynomial(i).coefficients == expected

    def test_eq(self):
        assert 0 == RationalNumber(0) == QuadraticNumber(0) == Polynomial(0)
        assert 1 == RationalNumber(1) == QuadraticNumber(1) == Polynomial(1)
        assert -7 == RationalNumber(-7) == QuadraticNumber(-7) == Polynomial(-7)

        q = RationalNumber(3, 4)
        assert q == QuadraticNumber(q) == Polynomial(q)

        q = RationalNumber(-3, 4)
        assert q == QuadraticNumber(q) == Polynomial(q)

        assert QuadraticNumber.sqrt(3) == Polynomial(QuadraticNumber.sqrt(3))
        assert QuadraticNumber.sqrt(-3) == Polynomial(QuadraticNumber.sqrt(-3))

    @pytest.mark.parametrize("a, b, expected", [
        (0, Polynomial('x'), True),
        (RationalNumber(0), Polynomial('x'), True),
        (QuadraticNumber(0), Polynomial('x'), True),
        (Polynomial(0), Polynomial('x'), True),
        (1, Polynomial('x'), False),
        (0, 2*Polynomial('x'), True),
        (RationalNumber(3, 5), Polynomial('x'), False),
        (QuadraticNumber.sqrt(3), Polynomial('x'), False),
        (0, Polynomial('x') * Polynomial('y'), True),
        (Polynomial('x'), Polynomial('x')**2 + 2*Polynomial('x'), True),
        (0, Polynomial('x') - 1, False),
        (0, Polynomial('x')**2 - 2*Polynomial('x') + 1, False),
        (-2*Polynomial('x'), Polynomial('x')**2 - 2*Polynomial('x') + 1, True),
        (Polynomial('x') + Polynomial('x')**2, 2*Polynomial('x') + 3*Polynomial('x')**2, True),
    ])
    def test_le(self, a, b, expected):
        """Test <= operator for Polynomials."""
        assert (a <= b) == (-b <= -a) == expected

        q = RationalNumber(127, 128)
        assert (q*a <= q*b) == (-q*b <= -q*a) == expected

        r = 1 + QuadraticNumber.sqrt(127)
        assert (r*a <= r*b) == (-r*b <= -r*a) == expected

        if (r-q)*b >= 0:
            assert (q*a <= r*b) == (-r*b <= -q*a) == expected

    @pytest.mark.parametrize("a, b, expected", [
        (0, Polynomial('x'), False),
        (0, 2*Polynomial('x'), False),
        (1, Polynomial('x'), False),
        (RationalNumber(3, 5), Polynomial('x'), False),
        (QuadraticNumber.sqrt(3), Polynomial('x'), False),
        (0, Polynomial('x') * Polynomial('y'), False),
        (1, Polynomial('x') * Polynomial('y') + 2, True),
        (RationalNumber(1), Polynomial('x') * Polynomial('y') + 2, True),
        (QuadraticNumber(1), Polynomial('x') * Polynomial('y') + 2, True),
        (Polynomial(1), Polynomial('x') * Polynomial('y') + 2, True),
        (Polynomial('x'), Polynomial('x')**2 + 2*Polynomial('x'), False),
        (0, Polynomial('x') + 1, True),
        (-2*Polynomial('x') - 3, Polynomial('x') - 1, True),
        (-2*Polynomial('x'), Polynomial('x')**2 - 2*Polynomial('x') + 1, True)
    ])
    def test_lt(self, a, b, expected):
        """Test < operator for Polynomials."""
        assert (a < b) == (-b < -a) == expected

        q = RationalNumber(127, 128)
        assert (q*a < q*b) == (-q*b < -q*a) == expected

        r = 1 + QuadraticNumber.sqrt(127)
        assert (r*a < r*b) == (-r*b < -r*a) == expected

        if (r-q)*b > 0:
            assert (q*a < r*b) == (-r*b < -q*a) == expected

    def test_hash(self):
        """Test that hashes for Polynomials are consistent."""
        assert hash(Polynomial(0)) == hash(0)
        assert hash(Polynomial(RationalNumber(7, 8))) == hash(RationalNumber(7, 8))
        assert hash(Polynomial(QuadraticNumber.sqrt(7))) == hash(QuadraticNumber.sqrt(7))

        x = Polynomial('x')
        assert hash(4*x) == hash(RationalNumber(4)*x) == hash(QuadraticNumber(4)*x)
        assert hash(4*x/3) == hash(RationalNumber(4, 3)*x) \
                           == hash(QuadraticNumber(RationalNumber(4, 3))*x)

    def test_getitem(self):
        """Test [] operator for Polynomials."""
        x = Polynomial('x')
        f = 3*x**3 + 2*x - 7

        assert f[1] == -7
        assert f['x'] == 2
        assert f[{ord('x'): 3}] == 3
        assert f[Monomial('x')] == 2

    def test_addition(self):
        """Tests for addition and substraction of Polynomials."""
        pass

    def test_multiplication(self):
        """Tests for multiplication and division of Polynomials."""
        pass

    def test_power(self):
        """Tests for exponentiation of Polynomials."""
        pass

    @pytest.mark.parametrize("i", [
        Polynomial(0),
        1.0,
        'abc',
        {'a': 1},
        Polynomial(10)
    ])
    def test_constructor_errors(self, i):
        """Test error handling of invalid arguments passed to Polynomial constructor."""
        try:
            Polynomial(i)
        except Exception as e:
            assert str(e).startswith('Invalid input to Polynomial: ')
        else:
            assert False

    def test_addition_errors(self):
        """Test error handling for invalid addition of Polynomials by zero."""
        pass

    def test_multiplication_errors(self):
        """Test error handling for invalid multiplication of Polynomials."""
        pass

    @pytest.mark.parametrize("divisor", [
        0,
        0.0,
        RationalNumber(0),
        QuadraticNumber(0),
        Polynomial(0),
        Polynomial('x')
    ])
    def test_division_errors(self, divisor):
        """Test error handling for attempted division of Polynomials by zero."""
        try:
            Polynomial('x') / divisor
        except:
            pass
        else:
            assert False

    @pytest.mark.parametrize("a, b", [
        (Polynomial(3), 1.2),
        (1.2, Polynomial(3)),
        (Polynomial(7), None),
        (None, Polynomial(7)),
    ])
    def test_operator_errors(self, a, b):
        """Tests for error handling of invalid operations involving Polynomials."""
        try:
            a + b
        except:
            pass
        else:
            assert False

        try:
            a - b
        except:
            pass
        else:
            assert False

        try:
            a * b
        except:
            pass
        else:
            assert False

        try:
            a / b
        except:
            pass
        else:
            assert False

    def test_power_errors(self):
        """Test error handling for invalid exponentiation of Polynomials."""
        pass
