#!/usr/bin/env python3
"""Rigorous certificate for the uniform degree-10 theta-moment bound.

Required dependency (the certificate checks this exact Python wrapper version):

    python3 -m pip install --target /tmp/theta-cert-deps python-flint==0.9.0
    PYTHONPATH=/tmp/theta-cert-deps python3 tmp/theta_moment_certificate.py

python-flint delegates real Taylor-series arithmetic to Arb.  Arb operations
return outward-rounded midpoint-radius balls.  Every comparison below is an
interval comparison: it returns True only when the complete ball lies on the
claimed side of the exact rational number.

Mathematical conventions
------------------------

Let

    H_{2k}(s) = sum_{j=0}^{2k} (-1)^j binom(2k,j)
                xi^(j)(s) xi^(2k-j)(s).

For s = 1/2-y and y>0, expansion of

    C_y(x) = -(1/y) d/ds [xi(s+ix) xi(s-ix)]

gives

    A_{2k}(y) = (-1)^k C_y^(2k)(0) = -H_{2k}'(s)/y.

The functional equation makes H_{2k}(1-s)=H_{2k}(s), so H'_{2k}(1/2)=0.
Consequently the endpoint formulas, including their signs and factors, are

    A_{2k}(0)   = H_{2k}''(1/2),
    A_{2k}(1/2) = -2 H_{2k}'(0).

The script constructs xi Taylor balls directly and applies these formulas.
"""

from __future__ import annotations

import math
import platform
from fractions import Fraction

import flint
from flint import arb, arb_series, ctx


REQUIRED_PYTHON_FLINT = "0.9.0"
DECIMAL_DIGITS = 80
SERIES_CAP = 18  # coefficients x^0,...,x^17; derivatives through order 12 used


def derivative(series: arb_series, order: int) -> arb_series:
    """Return the indicated formal derivative."""
    for _ in range(order):
        series = series.derivative()
    return series


def xi_taylor(center: str) -> arb_series:
    """Taylor ball for xi(center+x), with no pole cancellation in zeta/gamma.

    We use s*Gamma(s/2)/2 = Gamma(1+s/2), hence

        xi(s) = (s-1) pi^(-s/2) Gamma(1+s/2) zeta(s).

    This expression is analytic at both centers used here, 0 and 1/2.
    """
    s = arb_series([arb(center), arb(1)])
    log_pi = arb.pi().log()
    return (s - 1) * (-s * log_pi / 2).exp() * (1 + s / 2).gamma() * s.zeta()


def endpoint_moments(center: str) -> list[arb]:
    """Return A_0,...,A_10 at y=0 (center 1/2) or y=1/2 (center 0)."""
    xi = xi_taylor(center)
    moments: list[arb] = []
    for k in range(6):
        h = arb_series([arb(0)])
        for j in range(2 * k + 1):
            h += (
                (-1) ** j
                * math.comb(2 * k, j)
                * derivative(xi, j)
                * derivative(xi, 2 * k - j)
            )
        if center == "0.5":
            # y -> 0: -H'(1/2-y)/y -> H''(1/2).
            value = derivative(h, 2).coeffs()[0]
        elif center == "0":
            # y = 1/2: -H'(0)/(1/2) = -2 H'(0).
            value = -2 * derivative(h, 1).coeffs()[0]
        else:
            raise ValueError(center)
        moments.append(value)
    return moments


def rational_ball(numerator: int, denominator: int) -> arb:
    return arb(numerator) / denominator


def sharp_cross_endpoint_certificate(
    at_zero: list[arb], at_half: list[arb]
) -> tuple[list[arb], list[arb], arb, arb]:
    """Certify the optimal degree-10 cross-endpoint polynomial radius.

    The five coefficients are the exact endpoint ratios supplied by the
    analytic moment formulas, enclosed by Arb.  We prove that the derivative
    in t=x^2 is negative on 0<=t<=51 using its degree-four Bernstein balls.
    After shifting t=51+r, all five power-basis coefficient balls are negative,
    proving negativity for r>=0 as well.  We then certify opposite signs at the
    displayed decimal-rational endpoints.
    """
    bounds = [
        at_half[1] / at_zero[0],
        at_zero[2] / at_half[0],
        at_half[3] / at_zero[0],
        at_zero[4] / at_half[0],
        at_half[5] / at_zero[0],
    ]
    polynomial = [
        arb(1),
        -bounds[0] / 2,
        bounds[1] / math.factorial(4),
        -bounds[2] / math.factorial(6),
        bounds[3] / math.factorial(8),
        -bounds[4] / math.factorial(10),
    ]
    derivative_coefficients = [
        (degree + 1) * polynomial[degree + 1] for degree in range(5)
    ]
    power_on_unit_interval = [
        coefficient * 51**degree
        for degree, coefficient in enumerate(derivative_coefficients)
    ]
    bernstein = []
    for j in range(5):
        coefficient = arb(0)
        for k in range(j + 1):
            coefficient += (
                power_on_unit_interval[k]
                * math.comb(j, k)
                / math.comb(4, k)
            )
        assert coefficient < 0
        bernstein.append(coefficient)

    # For t >= 51, write Q'(t) = Q'(51+r).  If every power coefficient in
    # r is negative, then Q'(51+r) < 0 for every r >= 0.
    shifted_tail = []
    for degree in range(5):
        coefficient = arb(0)
        for original_degree in range(degree, 5):
            coefficient += (
                derivative_coefficients[original_degree]
                * math.comb(original_degree, degree)
                * 51 ** (original_degree - degree)
            )
        assert coefficient < 0
        shifted_tail.append(coefficient)

    root_lower = rational_ball(70362433, 10_000_000)
    root_upper = rational_ball(70362434, 10_000_000)

    def evaluate(x: arb) -> arb:
        return sum(
            (coefficient * x ** (2 * degree) for degree, coefficient in enumerate(polynomial)),
            arb(0),
        )

    assert evaluate(root_lower) > 0
    assert evaluate(root_upper) < 0
    return bernstein, shifted_tail, root_lower, root_upper


def bernstein_derivative_certificate() -> tuple[
    list[Fraction], Fraction, Fraction, Fraction
]:
    """Exact-fraction check through the first positive zero of P_*.

    Write t=x^2 and express d/dt P_*(sqrt(t)) in the degree-four Bernstein
    basis on t=51*s, 0<=s<=1. Negative Bernstein coefficients prove the
    derivative is negative throughout an interval containing the certified
    root bracket below.
    """
    polynomial = [
        Fraction(1),
        -Fraction(837, 10000) / 2,
        Fraction(199, 10000) / math.factorial(4),
        -Fraction(403, 50000) / math.factorial(6),
        Fraction(429, 100000) / math.factorial(8),
        -Fraction(1497, 500000) / math.factorial(10),
    ]

    # Power coefficients of q(t) = d/dt P_*(sqrt(t)).
    derivative_coefficients = [
        (degree + 1) * polynomial[degree + 1] for degree in range(5)
    ]
    # Substitute t=51*s.
    power_on_unit_interval = [
        coefficient * 51**degree
        for degree, coefficient in enumerate(derivative_coefficients)
    ]
    # Convert power coefficients to degree-four Bernstein coefficients:
    # s^k = sum_{j=k}^4 C(j,k)/C(4,k) B_{j,4}(s).
    bernstein = []
    for j in range(5):
        bernstein.append(
            sum(
                power_on_unit_interval[k]
                * Fraction(math.comb(j, k), math.comb(4, k))
                for k in range(j + 1)
            )
        )
    assert all(coefficient < 0 for coefficient in bernstein)

    x = Fraction(7)
    endpoint = sum(
        coefficient * x ** (2 * degree)
        for degree, coefficient in enumerate(polynomial)
    )
    assert endpoint == Fraction(299909357, 86400000000)
    assert endpoint > 0

    root_lower = Fraction(70165148, 10_000_000)
    root_upper = Fraction(70165149, 10_000_000)
    value_lower = sum(
        coefficient * root_lower ** (2 * degree)
        for degree, coefficient in enumerate(polynomial)
    )
    value_upper = sum(
        coefficient * root_upper ** (2 * degree)
        for degree, coefficient in enumerate(polynomial)
    )
    assert value_lower > 0
    assert value_upper < 0
    return bernstein, endpoint, root_lower, root_upper


def main() -> None:
    if flint.__version__ != REQUIRED_PYTHON_FLINT:
        raise RuntimeError(
            f"python-flint {REQUIRED_PYTHON_FLINT} required; "
            f"found {flint.__version__}"
        )

    # Setting dps changes the underlying binary precision.  Arb performs all
    # elementary/special-function and series operations with rigorous ball
    # enclosures at this precision.
    ctx.dps = DECIMAL_DIGITS
    ctx.cap = SERIES_CAP

    print("CERTIFICATE uniform degree-10 theta-moment bound")
    print("python", platform.python_version())
    print("python-flint", flint.__version__)
    print("FLINT", flint.__FLINT_VERSION__)
    print("ctx.dps", ctx.dps)
    print("ctx.prec_bits", ctx.prec)
    print("ctx.cap", ctx.cap)
    print("rounding Arb midpoint-radius balls with rigorous outward enclosures")
    print(
        "series calls arb_series.exp, gamma, zeta, log; "
        "xi coefficients through x^17, derivatives through order 12"
    )

    at_zero = endpoint_moments("0.5")
    at_half = endpoint_moments("0")

    print("\nA_{2k}(0) = H_{2k}''(1/2)")
    for k, value in enumerate(at_zero):
        print(f"k={k} {value}")
    print("\nA_{2k}(1/2) = -2 H_{2k}'(0)")
    for k, value in enumerate(at_half):
        print(f"k={k} {value}")

    # Cross-endpoint bounds.  Raw A_{2k}(y) are pointwise increasing in y,
    # hence A_{2k}(0)/A_0(1/2) <= m_{2k}(y)
    #       <= A_{2k}(1/2)/A_0(0).  No ratio monotonicity is assumed.
    comparisons = [
        (
            "m2 upper",
            at_half[1] / at_zero[0],
            rational_ball(837, 10000),
            "<",
        ),
        (
            "m4 lower",
            at_zero[2] / at_half[0],
            rational_ball(199, 10000),
            ">",
        ),
        (
            "m6 upper",
            at_half[3] / at_zero[0],
            rational_ball(403, 50000),
            "<",
        ),
        (
            "m8 lower",
            at_zero[4] / at_half[0],
            rational_ball(429, 100000),
            ">",
        ),
        (
            "m10 upper",
            at_half[5] / at_zero[0],
            rational_ball(1497, 500000),
            "<",
        ),
    ]

    print("\nInterval-vs-exact-rational comparisons")
    for name, value, bound, relation in comparisons:
        proved = value < bound if relation == "<" else value > bound
        print(f"{name}: {value} {relation} {bound} PROVED={proved}")
        assert proved

    sharp_bernstein, sharp_tail, sharp_root_lower, sharp_root_upper = (
        sharp_cross_endpoint_certificate(at_zero, at_half)
    )
    print("\nSharp cross-endpoint polynomial derivative Bernstein balls on t=51*s")
    for coefficient in sharp_bernstein:
        print(coefficient)
    print("\nSharp derivative power-basis balls after t=51+r")
    for coefficient in sharp_tail:
        print(coefficient)
    print("sharp globally unique positive root bracket", sharp_root_lower, sharp_root_upper)

    bernstein, endpoint, root_lower, root_upper = bernstein_derivative_certificate()
    print("\nExact Bernstein coefficients of d/dt P_*(sqrt(t)) on t=51*s")
    for coefficient in bernstein:
        print(coefficient)
    print("P_*(7)", endpoint)
    print("P_*(7) decimal", float(endpoint))
    print("unique positive root bracket", root_lower, root_upper)
    print("root bracket decimal", float(root_lower), float(root_upper))
    print("\nALL CERTIFICATE ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
