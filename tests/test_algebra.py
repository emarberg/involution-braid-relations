import pytest

from ..braids import RationalNumber, QuadraticNumber, Polynomial, PrimeFactorization


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

    @pytest.mark.parametrize("p, q", [
        (1, 0),
        (1.0, 20),
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

    def test_addition(self):
        """Tests for addition and substraction of RationalNumbers."""
        a = RationalNumber(3, 7)
        b = -a
        c = RationalNumber(4, 9)

        assert a + 0 == 0 + a == a
        assert a + b == b + a == 0 == RationalNumber()
        assert a + 1 == 1 + a == RationalNumber(10, 7)
        assert a - 2 == -2 + a == RationalNumber(-11, 7)
        assert a + c == c + a == RationalNumber(3*9 + 4*7, 7*9)
        assert a - c == -c + a == RationalNumber(3*9 - 4*7, 7*9)

        # cannot add RationalNumber and float
        try:
            a + 1.2
        except:
            pass
        else:
            assert False

    def test_multiplication(self):
        """Tests for multiplication and division of RationalNumbers."""
        a = RationalNumber(3, 7)
        b = -a
        c = RationalNumber(4, 9)

        assert a * 1 == 1 * a == a
        assert a * b == b * a == RationalNumber(-9, 49)
        assert a * -2 == -2 * a == RationalNumber(-6, 7)
        assert a * c == c * a == RationalNumber(4, 21)

        assert a / 2 == RationalNumber(3, 14)
        assert 2 / a == RationalNumber(14, 3)
        assert a / c == RationalNumber(27, 28)
        assert c / a == RationalNumber(28, 27)

        # cannot multiply RationalNumber and float
        try:
            a * -3.8
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

    def test_power(self):
        """Tests for exponentiation of RationalNumbers."""
        q = RationalNumber(3, 2)**4
        assert q.numerator == 81 and q.denominator == 16

        q = RationalNumber(3, 2)**-4
        assert q.numerator == 16 and q.denominator == 81

        q = RationalNumber(3, 2)**0
        assert q.numerator == 1 and q.denominator == 1

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
        (8*27*125*343, {2: 3, 3: 3, 5: 3, 7: 3}),
        (-3, {-1: 1, 3: 1}),
        (-2*3*5*7, {-1: 1, 2: 1, 3: 1, 5: 1, 7: 1}),
        (-8*27*125*343, {-1: 1, 2: 3, 3: 3, 5: 3, 7: 3}),
    ])
    def test_constructor(self, i, expected):
        pf = PrimeFactorization(i)
        assert pf.factorization == expected

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

    def test_multiplication(self):
        pass

    def test_get_square_free_part(self):
        pass

    def get_truncated_square_root(self):
        pass

    def test_get_divisor_exponent(self):
        pass

    def test_get_prime_factorization(self):
        pass


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

    @pytest.mark.parametrize("i", [
        (0.0),
        (QuadraticNumber()),
        (Polynomial())
    ])
    def test_constructor_errors(self, i):
        """Test error handling of invalid arguments passed to QuadraticNumber constructor."""
        try:
            QuadraticNumber(i)
        except Exception as e:
            assert str(e).startswith('Invalid input type to QuadraticNumber: ')
        else:
            assert False

    IRRATIONAL_SQUARE_1 = 3 + QuadraticNumber.sqrt(5)
    IRRATIONAL_SQUARE_2 = 7 + 3*QuadraticNumber.sqrt(5)

    @pytest.mark.parametrize("i, expected", [
        (0, {}),
        (RationalNumber(), {}),
        (QuadraticNumber(), {}),
        (1, {PrimeFactorization(1): 1}),
        (-4, {PrimeFactorization(-1): 2}),
        (8*27, {PrimeFactorization(6): 6}),
        (125, {PrimeFactorization(5): 5}),
        (RationalNumber(7, 8), {PrimeFactorization(14): RationalNumber(1, 4)}),
        (RationalNumber(-4, 3), {PrimeFactorization(-3): RationalNumber(2, 3)}),
        (RationalNumber(125), {PrimeFactorization(5): 5}),
        (QuadraticNumber(125), {PrimeFactorization(5): 5}),
        (QuadraticNumber(RationalNumber(125)), {PrimeFactorization(5): 5}),
        (QuadraticNumber(RationalNumber(7, 8)), {PrimeFactorization(14): RationalNumber(1, 4)}),
        (QuadraticNumber(RationalNumber(-4, 3)), {PrimeFactorization(-3): RationalNumber(2, 3)}),
        (IRRATIONAL_SQUARE_1, {
            PrimeFactorization(2): RationalNumber(1, 2),
            PrimeFactorization(10): RationalNumber(1, 2)
        }),
        (5 * IRRATIONAL_SQUARE_1 / 7, {
            PrimeFactorization(70): RationalNumber(1, 14),
            PrimeFactorization(14): RationalNumber(5, 14)
        }),
        (IRRATIONAL_SQUARE_2, {
            PrimeFactorization(2): RationalNumber(3, 2),
            PrimeFactorization(10): RationalNumber(1, 2)
        }),
        (RationalNumber(-4, 3) * IRRATIONAL_SQUARE_2, {
            PrimeFactorization(-6): RationalNumber(1),
            PrimeFactorization(-30): RationalNumber(1, 3)
        }),
    ])
    def test_sqrt(self, i, expected):
        """Test QuadraticNumber.sqrt with a variety of valid inputs."""
        n = QuadraticNumber.sqrt(i)
        assert n.coefficients == expected

    def test_hash(self):
        """Test that hashes for QuadraticNumbers are consistent."""
        pass

    def test_addition(self):
        """Tests for addition and substraction of QuadraticNumbers."""
        pass

    def test_multiplication(self):
        """Tests for multiplication and division of QuadraticNumbers."""
        pass

    def test_division_errors(self):
        """Test error handling for attempted division of QuadraticNumbers by zero."""
        pass

    def test_power(self):
        """Tests for exponentiation of QuadraticNumbers."""
        pass


class TestMonomial:

    def test_constructor(self):
        pass

    def test_eq(self):
        pass

    def test_multiplication(self):
        pass

    def test_power(self):
        pass


class TestPolynomial:

    def test_constructor(self):
        """Test constructor for Polynomial with a variety of valid inputs."""
        pass

    def test_constructor_errors(self):
        """Test error handling of invalid arguments passed to Polynomial constructor."""
        pass

    def test_hash(self):
        """Test that hashes for Polynomials are consistent."""
        pass

    def test_addition(self):
        """Tests for addition and substraction of Polynomials."""
        pass

    def test_multiplication(self):
        """Tests for multiplication and division of Polynomials."""
        pass

    def test_division_errors(self):
        """Test error handling for attempted division of Polynomials by zero."""
        pass

    def test_power(self):
        """Tests for exponentiation of Polynomials."""
        pass


class TestComparisons:
    pass
