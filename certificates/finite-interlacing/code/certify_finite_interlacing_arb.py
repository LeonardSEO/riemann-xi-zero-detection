#!/usr/bin/env python3
"""Certify a genuine two-step finite counting shift.

The calculation uses only the finite semilocal formulas frozen in
``../FROZEN_SPECIFICATION.md``.  At ``lambda^2=6`` and ``N=2`` it:

* includes exactly the prime powers 2, 3, 4=2^2, and 5;
* encloses the five independent Weil-form entries with Arb/Acb;
* isolates a simple inversion-even ground eigenvalue;
* proves the Dirichlet functional does not vanish on its line;
* encloses the two positive ZD-P active frequencies from the exact
  relative-determinant numerator; and
* certifies ``X_{sqrt(6),2}(8)=2``.

No zeta zero, limiting determinant, or cross-lambda input is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import flint
from flint import acb, arb, arb_mat, ctx


HERE = Path(__file__).resolve()
PHASE_ROOT = HERE.parents[1]
FROZEN_SPEC = PHASE_ROOT / "FROZEN_SPECIFICATION.md"
DEFAULT_OUTPUT = PHASE_ROOT / "certificates" / "finite_interlacing_arb.json"

PRECISION_BITS = 384
DISPLAY_DIGITS = 36
EPSILON_POWER = 40
TAIL_INTEGRAND_BOUND = 20
PRIME_POWERS = ((2, 2, 1), (3, 3, 1), (4, 2, 2), (5, 5, 1))
UNIQUE_ENTRIES = {
    "a0": (0, 0),
    "a1": (1, 1),
    "a2": (2, 2),
    "c1": (1, 0),
    "c2": (2, 0),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"finite interlacing certification failed: {message}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_record(value: arb) -> dict[str, object]:
    return {
        "ball": value.str(DISPLAY_DIGITS),
        "lower": value.lower().str(DISPLAY_DIGITS, radius=False),
        "upper": value.upper().str(DISPLAY_DIGITS, radius=False),
        "contains_zero": bool(value.contains(0)),
    }


def q_zero(a: int, b: int) -> arb:
    return arb(2) if a == b else arb(0)


def q_arb(a: int, b: int, y: arb, length: arb) -> arb:
    r = y / length
    if a == b:
        return arb(2) * (arb(1) - r) * (arb(2 * a) * r).cos_pi()
    return (
        (arb(2 * b) * r).sin_pi() - (arb(2 * a) * r).sin_pi()
    ) / (arb.pi() * arb(a - b))


def q_acb(a: int, b: int, y: acb, length: arb) -> acb:
    r = y / length
    if a == b:
        return 2 * (1 - r) * (2 * arb.pi() * a * r).cos()
    return (
        (2 * arb.pi() * b * r).sin()
        - (2 * arb.pi() * a * r).sin()
    ) / (arb.pi() * (a - b))


def archimedean_integral(
    a: int,
    b: int,
    length: arb,
    epsilon: arb,
) -> tuple[arb, acb]:
    q0 = q_zero(a, b)

    def integrand(y: acb, analytic: bool) -> acb:
        del analytic
        return (
            (y / 2).exp() * q_acb(a, b, y, length) - q0
        ) / (y.exp() - (-y).exp())

    middle = acb.integral(
        integrand,
        epsilon,
        length,
        rel_tol=arb(2) ** -250,
        abs_tol=arb(2) ** -250,
        deg_limit=300,
        eval_limit=150000,
        depth_limit=35,
        use_heap=True,
    )
    require(middle.is_finite(), f"nonfinite integral for {(a, b)}")
    require(middle.imag.contains(0), f"integral not certified real for {(a, b)}")

    # For |a|,|b| <= 2, |q| <= 2 and |q'| <= 34 on [0,epsilon].
    # Since exp(y/2)-1 <= y and 2*sinh(y) >= 2y there,
    # the removable integrand is bounded by 18.  Radius 20*epsilon is used.
    omitted = arb(0, arb(TAIL_INTEGRAND_BOUND) * epsilon)
    return middle.real + omitted, middle


def pole_term(a: int, b: int, length: arb) -> arb:
    pi = arb.pi()
    return (
        arb(32)
        * length
        * (length / 4).sinh() ** 2
        * (length**2 - arb(16) * pi**2 * a * b)
        / (
            (length**2 + arb(16) * pi**2 * a**2)
            * (length**2 + arb(16) * pi**2 * b**2)
        )
    )


def constant_archimedean(a: int, b: int, length: arb) -> arb:
    exponential = length.exp()
    constant = arb.const_euler() + (
        arb(4)
        * arb.pi()
        * (exponential - 1)
        / (exponential + 1)
    ).log()
    return -q_zero(a, b) * constant / 2


def prime_terms(a: int, b: int, length: arb) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for prime_power, prime, exponent in PRIME_POWERS:
        y = arb(exponent) * arb(prime).log()
        weight = arb(prime).log() / (arb(prime) ** exponent).sqrt()
        q_value = q_arb(a, b, y, length)
        contribution = -weight * q_value
        result.append(
            {
                "n": prime_power,
                "p": prime,
                "m": exponent,
                "weight": ball_record(weight),
                "q_value": ball_record(q_value),
                "signed_contribution": ball_record(contribution),
                "_value": contribution,
            }
        )
    return result


def qw_entry(
    a: int,
    b: int,
    length: arb,
    epsilon: arb,
) -> tuple[arb, dict[str, object]]:
    pole = pole_term(a, b, length)
    constant = constant_archimedean(a, b, length)
    integral, middle = archimedean_integral(a, b, length, epsilon)
    primes = prime_terms(a, b, length)
    prime_sum = sum((item["_value"] for item in primes), arb(0))
    value = pole + constant - integral + prime_sum
    for item in primes:
        del item["_value"]
    return value, {
        "indices": [a, b],
        "pole": ball_record(pole),
        "constant_archimedean": ball_record(constant),
        "integral": ball_record(integral),
        "middle_integral_imaginary": middle.imag.str(DISPLAY_DIGITS),
        "prime_terms": primes,
        "value": ball_record(value),
    }


def sorted_real_eigenvalues(matrix: arb_mat) -> list[acb]:
    values = matrix.eig(algorithm="rump")
    for value in values:
        require(value.imag.contains(0), "real symmetric eigenvalue missed real axis")
    return sorted(values, key=lambda value: float(value.real.mid()))


def build_certificate() -> dict[str, object]:
    ctx.prec = PRECISION_BITS
    length = arb(6).log()
    epsilon = arb(2) ** -EPSILON_POWER

    require(length.lower() > 1, "log(6)>1")
    require(arb.pi().upper() < 4, "pi<4")
    require(epsilon.lower() > 0 and epsilon.upper() < 1, "0<epsilon<1")
    require((arb(1) / 2).exp().upper() < 2, "exp(1/2)<2")

    entries: dict[str, arb] = {}
    entry_audits: dict[str, object] = {}
    for name, indices in UNIQUE_ENTRIES.items():
        entries[name], entry_audits[name] = qw_entry(
            *indices, length, epsilon
        )

    a0, a1, a2 = entries["a0"], entries["a1"], entries["a2"]
    c1, c2 = entries["c1"], entries["c2"]
    b1, b2 = c1, 2 * c2
    tau12 = (b1 - b2) / (1 - 2)
    tau1m2 = (b1 + b2) / (1 + 2)
    root2 = arb(2).sqrt()

    even = arb_mat(
        [
            [a0, root2 * c1, root2 * c2],
            [root2 * c1, a1 + c1, tau12 + tau1m2],
            [root2 * c2, tau12 + tau1m2, a2 + c2],
        ]
    )
    odd = arb_mat(
        [
            [a1 - c1, tau12 - tau1m2],
            [tau12 - tau1m2, a2 - c2],
        ]
    )

    even_values = sorted_real_eigenvalues(even)
    odd_values = sorted_real_eigenvalues(odd)
    ground = even_values[0].real
    other_lowers = [
        value.real.lower() for value in even_values[1:] + odd_values
    ]
    require(
        ground.upper() < min(other_lowers),
        "simple even ground interval is not separated",
    )

    raw_values, right_vectors = even.eig(right=True, algorithm="rump")
    ground_index = min(
        range(len(raw_values)),
        key=lambda index: float(raw_values[index].real.mid()),
    )
    first = right_vectors[0, ground_index]
    require(not first.contains(0), "ground V0 component may vanish")
    ratio1 = right_vectors[1, ground_index] / first
    ratio2 = right_vectors[2, ground_index] / first
    require(
        ratio1.imag.contains(0) and ratio2.imag.contains(0),
        "real ground ratios not enclosed",
    )

    xi0 = arb(1)
    xi1 = ratio1.real / root2
    xi2 = ratio2.real / root2
    coefficient_sum = xi0 + 2 * xi1 + 2 * xi2
    require(coefficient_sum.lower() > 0, "Dirichlet overlap may vanish")

    d1 = 2 * arb.pi() / length
    d2 = 2 * d1
    qa = -coefficient_sum
    qb = (
        xi0 * (d1**2 + d2**2)
        + 2 * xi1 * d2**2
        + 2 * xi2 * d1**2
    )
    qc = -xi0 * d1**2 * d2**2
    discriminant = qb**2 - 4 * qa * qc
    require(discriminant.lower() > 0, "active quadratic discriminant not positive")
    square_root = discriminant.sqrt()
    roots = sorted(
        [
            (-qb - square_root) / (2 * qa),
            (-qb + square_root) / (2 * qa),
        ],
        key=lambda value: float(value.mid()),
    )
    require(all(root.lower() > 0 for root in roots), "active x root not positive")
    active = [root.sqrt() for root in roots]
    require(
        active[0].lower() > d2.upper(),
        "smallest active frequency not beyond largest positive free frequency",
    )
    require(d2.upper() < 8 < active[0].lower(), "test point t=8 not separated")

    return {
        "certificate": "finite_interlacing_arb_lambda2_6_N2",
        "certified": True,
        "classification": (
            "actual semilocal counterexample to one-step interlacing: "
            "X_{sqrt(6),2}(8)=2"
        ),
        "implementation": {
            "python": sys.version.split()[0],
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
            "precision_bits": PRECISION_BITS,
            "script": HERE.relative_to(PHASE_ROOT).as_posix(),
            "script_sha256": sha256(HERE),
            "frozen_spec_sha256": sha256(FROZEN_SPEC),
        },
        "parameters": {
            "lambda_squared": 6,
            "N": 2,
            "L": ball_record(length),
            "epsilon": f"2^-{EPSILON_POWER}",
            "prime_power_endpoint": "p^m <= lambda^2, inclusive, no half weight",
            "retained_prime_powers": [
                {"n": n, "p": p, "m": m} for n, p, m in PRIME_POWERS
            ],
            "test_t": 8,
        },
        "scope": {
            "zeta_zeros_used": False,
            "RH_used": False,
            "cross_lambda_input_used": False,
            "determinant_convergence_used": False,
        },
        "improper_integral": {
            "split": f"[0,2^-{EPSILON_POWER}] union [2^-{EPSILON_POWER},log(6)]",
            "q_abs_bound": 2,
            "q_derivative_abs_bound": 34,
            "integrand_abs_bound": 18,
            "omitted_radius_used": f"{TAIL_INTEGRAND_BOUND}*2^-{EPSILON_POWER}",
        },
        "form_entries": entry_audits,
        "parity_eigenvalues": {
            "even": [ball_record(value.real) for value in even_values],
            "odd": [ball_record(value.real) for value in odd_values],
        },
        "ground_coefficients_scaled_xi0_1": {
            "xi0": ball_record(xi0),
            "xi_plus_minus_1": ball_record(xi1),
            "xi_plus_minus_2": ball_record(xi2),
            "sum_coefficients": ball_record(coefficient_sum),
        },
        "relative_determinant_numerator": {
            "variable": "x=z^2",
            "formula": (
                "-xi0*(d1^2-x)*(d2^2-x)"
                "+2*xi1*x*(d2^2-x)+2*xi2*x*(d1^2-x)"
            ),
            "discriminant": ball_record(discriminant),
            "x_roots": [ball_record(root) for root in roots],
            "positive_active_frequencies": [
                ball_record(frequency) for frequency in active
            ],
            "largest_positive_free_frequency_d2": ball_record(d2),
        },
        "counting_shift": {
            "t": 8,
            "N_reference_at_t": 4,
            "N_active_at_t": 2,
            "X_reference_minus_active": 2,
            "s_rank_C_minus_C0": 2,
            "rank_A_minus_B": 4,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    certificate = build_certificate()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check:
        require(args.output.exists(), f"missing certificate {args.output}")
        require(args.output.read_text() == encoded, "certificate is not reproducible")
        print(f"PASS: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(args.output)


if __name__ == "__main__":
    main()
