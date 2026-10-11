#!/usr/bin/env python3
"""Full-theta Arb certificate for a positive N=3 Möbius cross increment.

Uses the already available python-flint 0.9.0, without installing a runtime:
  PYTHONPATH=/tmp/riemann-flint python3 synthesis/certificates/296_theta_recursion_arb.py

Pass --output PATH only to save a new result. Default execution never edits JSON.
All quadrature cells are exact dyadic intervals. The omitted theta n>=2
terms and the complete two-sided physical integral tails are bounded below.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import flint
from flint import arb, ctx


def phi_bounds(u):
    """Full Phi(|u|) lower/upper on an Arb real interval."""
    u = abs(u)
    tt = arb.pi() * (2 * u).exp()
    prefactor = 4 * tt * tt - 6 * tt
    assert prefactor.lower() > 0
    head = prefactor * (u / 2 - tt).exp()
    # Every n>=2 term is positive. Ratios of n^4 exp(-n^2 T)
    # are <=(3/2)^4 exp(-5T)<1/2, so their sum <=32 exp(-4T).
    theta_tail = 128 * tt * tt * (u / 2 - 4 * tt).exp()
    # A very wide exponential enclosure may round its positive lower
    # endpoint to zero (or slightly negative). True positivity permits zero.
    lower = head.lower() if head.lower() > 0 else arb(0)
    return lower, head.upper() + theta_tail.upper()


def global_tail(x):
    """Upper for integral_{|r|>x} Phi(r)Phi(r+s)dr, for ANY real s."""
    pi = arb.pi()
    tt = pi * (2 * x).exp()
    assert tt > arb(5) / 4
    # Phi<=8 pi^2 exp(9|r|/2-pi exp(2|r|)) and sup Phi<4.
    # Gamma(9/4,T)<=T^(5/4)exp(-T)/(1-5/(4T)); elementary log bound.
    return (
        32
        / pi.sqrt().sqrt()
        * tt
        * tt.sqrt().sqrt()
        * (-tt).exp()
        / (1 - arb(5) / (4 * tt))
    ).upper()


def gram_bounds(s, cells, x=arb(2)):
    width = 2 * x / cells
    lower = arb(0)
    upper = arb(0)
    for k in range(cells):
        lo = -x + k * width
        cell = arb(lo + width / 2, width / 2)
        p_lo, p_hi = phi_bounds(cell)
        q_lo, q_hi = phi_bounds(cell + s)
        lower += width * p_lo * q_lo
        upper += width * p_hi * q_hi
    upper += global_tail(x)
    return lower.lower(), upper.upper()


def enclosing(lo, hi):
    assert hi > lo
    return arb((lo + hi) / 2, (hi - lo) / 2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cells", type=int, default=1024)
    args = parser.parse_args()
    if args.cells < 1024 or args.cells & (args.cells - 1):
        raise ValueError("cells must be a power of two at least 1024")
    ctx.prec = 192
    started = time.perf_counter()
    pi = arb.pi()
    assert arb(81) / 16 * (-5 * pi).exp() < arb(1) / 2
    assert 32 * (-3 * pi).exp() < 1
    assert 8 * pi * pi * (-pi).exp() < 4
    assert pi > arb(9) / 4
    grams = {}
    for label, s in [
        ("zero", arb(0)),
        ("log_3", arb(3).log()),
        ("log_3_over_2", (arb(3) / 2).log()),
    ]:
        lo, hi = gram_bounds(s, args.cells)
        grams[label] = enclosing(lo, hi)
    cross = (
        grams["log_3_over_2"] / arb(2).sqrt() - grams["log_3"]
    ) / arb(3).sqrt()
    diagonal = grams["zero"] / 3
    increment = 2 * cross + diagonal
    assert cross > arb(1) / 32
    assert increment - diagonal > arb(1) / 16
    result = {
        "status": "PASS",
        "scope": "raw full Phi, N=3 signed cross increment; no RH assertion",
        "fourier_convention": "F_f(z)=integral f(x) exp(-izx) dx",
        "precision_bits": ctx.prec,
        "python_flint_version": flint.__version__,
        "python_version": platform.python_version(),
        "quadrature_cells": args.cells,
        "quadrature_interval": "[-2,2]",
        "gram_balls": {k: str(v) for k, v in grams.items()},
        "cross_ball": str(cross),
        "cross_strict_lower": "1/32",
        "diagonal_ball": str(diagonal),
        "energy_increment_ball": str(increment),
        "increment_minus_diagonal_strict_lower": "1/16",
        "complete_physical_tail_upper": str(global_tail(arb(2))),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "measured_runtime_seconds": time.perf_counter() - started,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
