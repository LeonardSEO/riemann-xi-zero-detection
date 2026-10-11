#!/usr/bin/env python3
"""Rigorous Arb/Flint upper bound for the N=1 F1 tail constant K0.

Run without changing project dependencies:
    uv run --offline --with python-flint python3 \
      synthesis/route_reports/n1_tail_bound.py \
      --certificate synthesis/route_reports/n1_tail_bound_certificate.json

The Riemann sum uses Arb interval polynomial evaluation on [3,100].
The replacement pi >= 3 is a majorant.  The tail [100,infinity) uses
absolute polynomial coefficients and Arb upper incomplete gamma values.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction

from flint import arb, ctx


def poly_add(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for j, value in enumerate(a):
        out[j] += value
    for j, value in enumerate(b):
        out[j] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_scale(a, c):
    return [value * c for value in a]


def poly_deriv(a):
    return [a[j + 1] * (j + 1) for j in range(len(a) - 1)] or [Fraction(0)]


def poly_mul_t(a):
    return [Fraction(0)] + a


def recurrence(p):
    # P_{r+1}=(1/2-2T)P_r+2T P_r'
    return poly_add(
        poly_scale(p, Fraction(1, 2)),
        poly_add(poly_scale(poly_mul_t(p), -2),
                 poly_scale(poly_mul_t(poly_deriv(p)), 2)),
    )


def eval_poly(p, t):
    value = arb(0)
    for coefficient in reversed(p):
        value = value * t + arb(coefficient.numerator) / coefficient.denominator
    return value


def rational_arb(num, den):
    return arb(f"{num}/{den}")


def interval_upper_sum(n, p4):
    """Upper bound for integral_3^100 T^-1/2 |P4(T)| exp(-T) dT."""
    width = arb(97) / n
    total = arb(0)
    for j in range(n):
        lo = rational_arb(3 * n + 97 * j, n)
        hi = rational_arb(3 * n + 97 * (j + 1), n)
        midpoint = (lo + hi) / 2
        radius = (hi - lo) / 2
        cell = arb(midpoint, radius)
        polynomial_bound = abs(eval_poly(p4, cell))
        inverse_sqrt_bound = 1 / lo.sqrt()
        exponential_bound = (-lo).exp()
        total += width * polynomial_bound * inverse_sqrt_bound * exponential_bound
    return total


def tail_upper_100(p4):
    """Absolute-coefficient upper bound for integral_100^inf."""
    total = arb(0)
    for j, coefficient in enumerate(p4):
        # Integral T^(j-1/2) exp(-T) dT = Gamma(j+1/2,100).
        total += abs(arb(coefficient.numerator) / coefficient.denominator) * (
            arb(100).gamma_upper(arb(j) + arb(1) / 2)
        )
    return total


def upper_string(value):
    return str(value.upper())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100_000)
    parser.add_argument("--certificate", type=str)
    args = parser.parse_args()

    if args.n <= 0:
        raise SystemExit("--n must be positive")
    ctx.dps = 80

    # P0=4T^2-6T and the stated derivative recurrence.
    polys = [[Fraction(0), Fraction(-6), Fraction(4)]]
    for _ in range(4):
        polys.append(recurrence(polys[-1]))
    # Exact normalization checks: P0=4T^2-6T and P1=-8T^3+30T^2-15T.
    assert polys[0] == [Fraction(0), Fraction(-6), Fraction(4)]
    assert polys[1] == [Fraction(0), Fraction(-15), Fraction(30), Fraction(-8)]
    # The signed P1 value gives f'(0+) = pi*exp(-pi)*(30*pi-15-8*pi^2)
    # because f'=exp(u/2-T)P1 and P1(pi)=pi*(-15+30*pi-8*pi^2).
    p3, p4 = polys[3], polys[4]

    pi = arb.pi()
    a1_formula = pi * (-pi).exp() * (30 * pi - 15 - 8 * pi * pi)
    p1_pi = eval_poly(polys[1], pi)
    a1_from_p1 = (-pi).exp() * p1_pi
    # Overlap of the two Arb balls is required, not just numerical agreement.
    if not (a1_formula.contains(a1_from_p1)) and not (a1_from_p1.contains(a1_formula)):
        raise AssertionError('P1(pi) and analytical a1 balls do not overlap')
    a1 = a1_formula
    if not (a1 > 0):
        raise AssertionError("a1 must have a strictly positive Arb lower bound")
    f3 = (-pi).exp() * eval_poly(p3, pi)
    f3_abs_upper = abs(f3)

    core_integral_upper = interval_upper_sum(args.n, p4)
    tail_integral_upper = tail_upper_100(p4)
    weighted_f4_upper = (
        core_integral_upper + tail_integral_upper
    ) / (2 * pi.sqrt())
    K0_upper = 2 * (f3_abs_upper + weighted_f4_upper)

    ratio_upper = K0_upper / (2 * a1)
    # Choose the first integer R with a strict interval proof R^2 > ratio_upper.
    R = 1
    while not (arb(R * R) > ratio_upper):
        R += 1
    if not (arb(R * R) > ratio_upper):
        raise AssertionError("strict R^2 > K0/(2*a1) interval proof failed")

    result = {
        "precision_bits_or_dps": 80,
        "n_cells": args.n,
        "polynomial_P0_coefficients_low_to_high": [str(c) for c in polys[0]],
        "polynomial_P1_coefficients_low_to_high": [str(c) for c in polys[1]],
        "polynomial_P3_coefficients_low_to_high": [str(c) for c in p3],
        "polynomial_P4_coefficients_low_to_high": [str(c) for c in p4],
        "a1_ball": str(a1),
        "a1_from_P1_ball": str(a1_from_p1),
        "a1_overlap_checked": True,
        "a1_lower": str(a1.lower()),
        "abs_f3_upper": str(f3_abs_upper),
        "integral_3_100_upper_before_factor": str(core_integral_upper.upper()),
        "integral_100_inf_upper_before_factor": str(tail_integral_upper.upper()),
        "weighted_f4_upper": str(weighted_f4_upper.upper()),
        "K0_upper": str(K0_upper.upper()),
        "K0_over_2a1_upper": str(ratio_upper.upper()),
        "R": R,
        "R_squared": R * R,
        "claim": (
            "Using |F1(z)+2*a1/z^2| <= K0*|z|^-4 on "
            "|Im z|<=1/2, all F1 zeros in that strip have |Re z|<R."
        ),
        "hypotheses": [
            "T=pi*exp(2u), pi>3; [pi,infinity] is majorized by [3,infinity].",
            "The [3,100] integral is bounded by interval-cell suprema.",
            "The [100,infinity] integral uses absolute polynomial coefficients.",
            "The asymptotic remainder bound with this K0 is the separately proved IBP bound.",
        ],
    }
    print(json.dumps(result, indent=2))
    if args.certificate:
        with open(args.certificate, "w", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")


if __name__ == "__main__":
    main()
