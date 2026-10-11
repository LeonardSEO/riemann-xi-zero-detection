#!/usr/bin/env python3
"""Independent interval certificate for the two-step finite count.

This implementation deliberately does not diagonalize an Arb matrix and does
not solve the active-frequency quadratic.  It instead

* evaluates nine form entries directly (including every retained prime power);
* constructs the even and odd parity blocks without the displacement shortcut;
* certifies their inertia by an interval LDL congruence;
* encloses the projective ground line from a two-by-two linear solve; and
* proves both active squared frequencies exceed 64 from certified signs of
  their characteristic polynomial and its derivative.

The conclusion is the finite statement

    X_{sqrt(6),2}(8) = N_reference(8) - N_active(8) = 2.

It uses no Riemann zero, RH assumption, limiting determinant, or cross-lambda
input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import flint
from flint import acb, arb, arb_mat, ctx


SCRIPT = Path(__file__).resolve()
PHASE_ROOT = SCRIPT.parents[1]
FROZEN_SPECIFICATION = PHASE_ROOT / "FROZEN_SPECIFICATION.md"
DEFAULT_CERTIFICATE = (
    PHASE_ROOT / "certificates" / "implementation_c_finite_replica.json"
)

PRECISION_BITS = 448
PRINT_DIGITS = 42
ENDPOINT_POWER = 42
OMITTED_INTEGRAND_BOUND = 20

# p^m <= 6, inclusive.  Six itself is not a prime power.
PRIME_POWERS = (
    {"n": 2, "p": 2, "m": 1},
    {"n": 3, "p": 3, "m": 1},
    {"n": 4, "p": 2, "m": 2},
    {"n": 5, "p": 5, "m": 1},
)

# These are evaluated directly.  In particular, the off-axis entries are not
# reconstructed from the displacement identity used by implementation A.
DIRECT_ENTRY_PAIRS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 1),
    (1, -1),
    (1, 2),
    (1, -2),
    (2, 2),
    (2, -2),
)

GROUND_LOWER = arb("0.0000000138")
GROUND_UPPER = arb("0.0000000141")
SEPARATION_LEVEL = arb("0.0000001")
COUNTING_POINT = arb(8)
SQUARED_COUNTING_POINT = arb(64)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"implementation C certification failed: {message}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_sign(value: arb) -> str:
    if value.lower() > 0:
        return "positive"
    if value.upper() < 0:
        return "negative"
    if value == 0:
        return "zero"
    return "indeterminate"


def ball(value: arb) -> dict[str, Any]:
    return {
        "ball": value.str(PRINT_DIGITS),
        "lower": value.lower().str(PRINT_DIGITS, radius=False),
        "upper": value.upper().str(PRINT_DIGITS, radius=False),
        "contains_zero": bool(value.contains(0)),
        "strict_sign": strict_sign(value),
    }


def overlap_real(j: int, k: int, y: arb, length: arb) -> arb:
    """The exact Fourier overlap q_{jk}(y)."""

    ratio = y / length
    if j == k:
        return 2 * (1 - ratio) * (2 * j * ratio).cos_pi()
    return (
        (2 * k * ratio).sin_pi() - (2 * j * ratio).sin_pi()
    ) / (arb.pi() * (j - k))


def overlap_complex(j: int, k: int, y: acb, length: arb) -> acb:
    ratio = y / length
    if j == k:
        return 2 * (1 - ratio) * (2 * arb.pi() * j * ratio).cos()
    return (
        (2 * arb.pi() * k * ratio).sin()
        - (2 * arb.pi() * j * ratio).sin()
    ) / (arb.pi() * (j - k))


def pole_piece(j: int, k: int, length: arb) -> arb:
    pi = arb.pi()
    return (
        32
        * length
        * (length / 4).sinh() ** 2
        * (length**2 - 16 * pi**2 * j * k)
        / (
            (length**2 + 16 * pi**2 * j**2)
            * (length**2 + 16 * pi**2 * k**2)
        )
    )


def archimedean_constant(j: int, k: int, length: arb) -> arb:
    q_at_zero = arb(2) if j == k else arb(0)
    exponential = length.exp()
    normalization = arb.const_euler() + (
        4 * arb.pi() * (exponential - 1) / (exponential + 1)
    ).log()
    return -q_at_zero * normalization / 2


def archimedean_integral(
    j: int,
    k: int,
    length: arb,
    endpoint: arb,
) -> tuple[arb, list[dict[str, Any]]]:
    """Integrate in three independent segments and enclose the removed tip."""

    q_at_zero = arb(2) if j == k else arb(0)

    def integrand(y: acb, analytic: bool) -> acb:
        del analytic
        numerator = (y / 2).exp() * overlap_complex(j, k, y, length)
        return (numerator - q_at_zero) / (y.exp() - (-y).exp())

    cuts = (endpoint, length / 3, 2 * length / 3, length)
    pieces: list[acb] = []
    for left, right in zip(cuts[:-1], cuts[1:]):
        value = acb.integral(
            integrand,
            left,
            right,
            rel_tol=arb(2) ** -300,
            abs_tol=arb(2) ** -300,
            deg_limit=320,
            eval_limit=180000,
            depth_limit=36,
            use_heap=True,
        )
        require(value.is_finite(), f"finite integral segment {(j, k)}")
        require(value.imag.contains(0), f"real integral segment {(j, k)}")
        pieces.append(value)

    # For |j|,|k| <= 2, the direct elementary estimates
    # |q_{jk}| <= 2 and |q'_{jk}| <= 34 on [0,2^-42], together with
    # exp(y/2)-1 <= y and 2*sinh(y) >= 2y, bound the removable
    # integrand by 18.  We use the strictly larger radius 20*2^-42.
    omitted = arb(0, OMITTED_INTEGRAND_BOUND * endpoint)
    total = sum((piece.real for piece in pieces), arb(0)) + omitted
    ledger = [
        {
            "interval": [
                cuts[index].str(PRINT_DIGITS),
                cuts[index + 1].str(PRINT_DIGITS),
            ],
            "real_enclosure": ball(piece.real),
            "imaginary_enclosure": piece.imag.str(PRINT_DIGITS),
        }
        for index, piece in enumerate(pieces)
    ]
    ledger.append(
        {
            "interval": ["0", endpoint.str(PRINT_DIGITS)],
            "analytic_bound": f"|integrand| <= 18; radius {OMITTED_INTEGRAND_BOUND}*2^-{ENDPOINT_POWER}",
            "real_enclosure": ball(omitted),
        }
    )
    return total, ledger


def semilocal_entry(
    j: int,
    k: int,
    length: arb,
    endpoint: arb,
) -> tuple[arb, dict[str, Any]]:
    pole = pole_piece(j, k, length)
    constant = archimedean_constant(j, k, length)
    integral, integral_ledger = archimedean_integral(
        j, k, length, endpoint
    )

    arithmetic_total = arb(0)
    prime_ledger: list[dict[str, Any]] = []
    for item in PRIME_POWERS:
        p = item["p"]
        m = item["m"]
        location = m * arb(p).log()
        weight = arb(p).log() / (arb(p) ** m).sqrt()
        overlap = overlap_real(j, k, location, length)
        contribution = -weight * overlap
        arithmetic_total += contribution
        prime_ledger.append(
            {
                **item,
                "location_m_log_p": ball(location),
                "coefficient_log_p_over_sqrt_p_to_m": ball(weight),
                "overlap": ball(overlap),
                "signed_contribution": ball(contribution),
            }
        )

    value = pole + constant - integral + arithmetic_total
    return value, {
        "indices": [j, k],
        "pole": ball(pole),
        "archimedean_constant": ball(constant),
        "archimedean_integral": ball(integral),
        "archimedean_segments": integral_ledger,
        "prime_power_terms": prime_ledger,
        "prime_power_total": ball(arithmetic_total),
        "form_entry": ball(value),
    }


def matrix_ledger(matrix: arb_mat) -> list[list[dict[str, Any]]]:
    return [
        [ball(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ]


def shifted_ldl_inertia(
    matrix: arb_mat,
    shift: arb,
) -> tuple[int, list[arb]]:
    """Certified inertia of matrix-shift*I via interval LDL congruence."""

    dimension = matrix.nrows()
    work = [
        [
            matrix[row, column] - (shift if row == column else 0)
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    pivots: list[arb] = []
    negative = 0
    for pivot_index in range(dimension):
        pivot = work[pivot_index][pivot_index]
        require(
            not pivot.contains(0),
            f"LDL pivot {pivot_index} excludes zero at shift {shift}",
        )
        pivots.append(pivot)
        if pivot.upper() < 0:
            negative += 1
        else:
            require(pivot.lower() > 0, "LDL pivot has strict sign")

        for row in range(pivot_index + 1, dimension):
            for column in range(row, dimension):
                reduced = (
                    work[row][column]
                    - work[row][pivot_index]
                    * work[column][pivot_index]
                    / pivot
                )
                work[row][column] = reduced
                work[column][row] = reduced
    return negative, pivots


def inertia_record(
    label: str,
    matrix: arb_mat,
    shift: arb,
) -> dict[str, Any]:
    negative, pivots = shifted_ldl_inertia(matrix, shift)
    return {
        "label": label,
        "shift": ball(shift),
        "negative_pivots": negative,
        "pivots": [ball(pivot) for pivot in pivots],
        "meaning": (
            "number of eigenvalues strictly below the shift, by "
            "Sylvester inertia under the certified LDL congruence"
        ),
    }


def build_certificate() -> dict[str, Any]:
    ctx.prec = PRECISION_BITS
    length = arb(6).log()
    endpoint = arb(2) ** -ENDPOINT_POWER

    require(length.lower() > 1, "log(6)>1")
    require(length.upper() < 2, "log(6)<2")
    require(arb.pi().upper() < 4, "pi<4")
    require(endpoint.lower() > 0, "positive integral endpoint")

    entries: dict[tuple[int, int], arb] = {}
    entry_audits: dict[str, Any] = {}
    for pair in DIRECT_ENTRY_PAIRS:
        value, audit = semilocal_entry(*pair, length, endpoint)
        entries[pair] = value
        entry_audits[f"{pair[0]},{pair[1]}"] = audit

    root_two = arb(2).sqrt()
    even = arb_mat(
        [
            [
                entries[(0, 0)],
                root_two * entries[(0, 1)],
                root_two * entries[(0, 2)],
            ],
            [
                root_two * entries[(0, 1)],
                entries[(1, 1)] + entries[(1, -1)],
                entries[(1, 2)] + entries[(1, -2)],
            ],
            [
                root_two * entries[(0, 2)],
                entries[(1, 2)] + entries[(1, -2)],
                entries[(2, 2)] + entries[(2, -2)],
            ],
        ]
    )
    odd = arb_mat(
        [
            [
                entries[(1, 1)] - entries[(1, -1)],
                entries[(1, 2)] - entries[(1, -2)],
            ],
            [
                entries[(1, 2)] - entries[(1, -2)],
                entries[(2, 2)] - entries[(2, -2)],
            ],
        ]
    )

    inertia = {
        "even_below_ground_lower": inertia_record(
            "even", even, GROUND_LOWER
        ),
        "even_below_ground_upper": inertia_record(
            "even", even, GROUND_UPPER
        ),
        "odd_below_ground_upper": inertia_record(
            "odd", odd, GROUND_UPPER
        ),
        "even_below_separation_level": inertia_record(
            "even", even, SEPARATION_LEVEL
        ),
        "odd_below_separation_level": inertia_record(
            "odd", odd, SEPARATION_LEVEL
        ),
    }
    require(
        inertia["even_below_ground_lower"]["negative_pivots"] == 0,
        "no even eigenvalue below ground lower endpoint",
    )
    require(
        inertia["even_below_ground_upper"]["negative_pivots"] == 1,
        "exactly one even eigenvalue below ground upper endpoint",
    )
    require(
        inertia["odd_below_ground_upper"]["negative_pivots"] == 0,
        "no odd eigenvalue below ground upper endpoint",
    )
    require(
        inertia["even_below_separation_level"]["negative_pivots"] == 1,
        "only one even eigenvalue below separation level",
    )
    require(
        inertia["odd_below_separation_level"]["negative_pivots"] == 0,
        "no odd eigenvalue below separation level",
    )

    ground = arb(
        (GROUND_LOWER + GROUND_UPPER) / 2,
        (GROUND_UPPER - GROUND_LOWER) / 2,
    )
    block_11 = even[1, 1] - ground
    block_12 = even[1, 2]
    block_22 = even[2, 2] - ground
    bottom_determinant = block_11 * block_22 - block_12**2
    require(
        not bottom_determinant.contains(0),
        "ground-vector bottom system invertible",
    )

    coordinate_1 = (
        -even[0, 1] * block_22 + block_12 * even[0, 2]
    ) / bottom_determinant
    coordinate_2 = (
        -even[0, 2] * block_11 + block_12 * even[0, 1]
    ) / bottom_determinant

    xi_0 = arb(1)
    xi_1 = coordinate_1 / root_two
    xi_2 = coordinate_2 / root_two
    coefficient_sum = xi_0 + 2 * xi_1 + 2 * xi_2
    require(
        coefficient_sum.lower() > 0,
        "Dirichlet/source overlap is strictly nonzero",
    )

    d_1 = 2 * arb.pi() / length
    d_2 = 2 * d_1
    require(d_2.upper() < COUNTING_POINT, "both free positive modes below 8")

    # P(x) is the exact numerator in x=z^2:
    #   -xi0(d1^2-x)(d2^2-x)
    #   +2 xi1 x(d2^2-x)+2 xi2 x(d1^2-x).
    leading = -coefficient_sum
    linear = (
        xi_0 * (d_1**2 + d_2**2)
        + 2 * xi_1 * d_2**2
        + 2 * xi_2 * d_1**2
    )
    constant = -xi_0 * d_1**2 * d_2**2
    discriminant = linear**2 - 4 * leading * constant
    value_at_64 = (
        leading * SQUARED_COUNTING_POINT**2
        + linear * SQUARED_COUNTING_POINT
        + constant
    )
    derivative_at_64 = 2 * leading * SQUARED_COUNTING_POINT + linear

    require(leading.upper() < 0, "active polynomial leading coefficient <0")
    require(linear.lower() > 0, "active polynomial linear coefficient >0")
    require(constant.upper() < 0, "active polynomial constant <0")
    require(discriminant.lower() > 0, "active polynomial has two real roots")
    require(value_at_64.upper() < 0, "P(64)<0")
    require(derivative_at_64.lower() > 0, "P'(64)>0")

    return {
        "certificate": "implementation_c_finite_replica_lambda2_6_N2",
        "certified": True,
        "theorem": "X_{sqrt(6),2}(8)=2",
        "method": (
            "direct nine-entry parity blocks; interval LDL inertia; "
            "two-by-two projective ground solve; characteristic-polynomial "
            "sign proof without matrix diagonalization or quadratic roots"
        ),
        "implementation": {
            "python": sys.version.split()[0],
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
            "precision_bits": PRECISION_BITS,
            "script": SCRIPT.relative_to(PHASE_ROOT).as_posix(),
            "script_sha256": sha256(SCRIPT),
            "frozen_specification_sha256": sha256(FROZEN_SPECIFICATION),
        },
        "parameters": {
            "lambda_squared": 6,
            "N": 2,
            "L_log_6": ball(length),
            "integral_endpoint": f"2^-{ENDPOINT_POWER}",
            "prime_power_rule": "p^m <= lambda^2, inclusive, no half weight",
            "retained_prime_powers": list(PRIME_POWERS),
            "counting_point_t": 8,
        },
        "direct_form_entries": entry_audits,
        "parity_blocks": {
            "basis_even": [
                "V_0",
                "(V_1+V_-1)/sqrt(2)",
                "(V_2+V_-2)/sqrt(2)",
            ],
            "basis_odd": [
                "(V_1-V_-1)/sqrt(2)",
                "(V_2-V_-2)/sqrt(2)",
            ],
            "even": matrix_ledger(even),
            "odd": matrix_ledger(odd),
        },
        "ground_state_certificate": {
            "inertia": inertia,
            "ground_eigenvalue_enclosure": ball(ground),
            "simple_isolated_inversion_even": True,
            "next_spectrum_strictly_above": ball(SEPARATION_LEVEL),
            "bottom_system_determinant": ball(bottom_determinant),
            "parity_coordinates_scaled_V0_1": {
                "V0": ball(xi_0),
                "even_mode_1": ball(coordinate_1),
                "even_mode_2": ball(coordinate_2),
            },
            "fourier_coefficients_scaled_xi0_1": {
                "xi_0": ball(xi_0),
                "xi_plus_minus_1": ball(xi_1),
                "xi_plus_minus_2": ball(xi_2),
                "coefficient_sum": ball(coefficient_sum),
            },
            "dirichlet_normalization_nonzero": True,
            "finite_pair_admissible": True,
        },
        "active_frequency_certificate": {
            "free_positive_frequencies": {
                "d1": ball(d_1),
                "d2": ball(d_2),
                "d2_strictly_below_8": True,
            },
            "polynomial_variable": "x=z^2",
            "polynomial_formula": (
                "-xi0*(d1^2-x)*(d2^2-x)"
                "+2*xi1*x*(d2^2-x)+2*xi2*x*(d1^2-x)"
            ),
            "coefficients": {
                "x_squared": ball(leading),
                "x": ball(linear),
                "constant": ball(constant),
            },
            "discriminant": ball(discriminant),
            "P_at_64": ball(value_at_64),
            "P_prime_at_64": ball(derivative_at_64),
            "proof": (
                "The leading and constant coefficients are negative, the "
                "linear coefficient is positive, and the discriminant is "
                "positive, so the two roots are positive.  The polynomial "
                "is concave; P(64)<0 and P'(64)>0 place 64 strictly to the "
                "left of its smaller root.  Hence both positive active "
                "frequencies are sqrt(root)>8."
            ),
            "both_positive_active_frequencies_strictly_above_8": True,
        },
        "counting_shift": {
            "right_continuous_convention": "N_T(t)=Tr 1_{(-infinity,t]}(T)",
            "t": 8,
            "reference_spectrum_order": "(-d2,-d1,d1,d2)",
            "active_spectrum_order": "(-gamma2,-gamma1,gamma1,gamma2)",
            "N_reference_at_t": 4,
            "N_active_at_t": 2,
            "X_reference_minus_active": 2,
        },
        "scope": {
            "finite_only": True,
            "certifies_uniform_unboundedness": False,
            "RH_used": False,
            "zeta_zeros_used": False,
            "cross_lambda_input_used": False,
            "limiting_determinant_used": False,
            "numerical_diagnostics_are_load_bearing": True,
            "claim": (
                "This is a certified finite counterexample to a universal "
                "one-step interlacing claim.  It neither proves unbounded "
                "finite shifts nor any cross-lambda statement."
            ),
        },
    }


def serialized_certificate() -> str:
    return json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    rendered = serialized_certificate()
    if arguments.check:
        require(arguments.output.exists(), f"missing {arguments.output}")
        require(
            arguments.output.read_text() == rendered,
            f"{arguments.output} is stale",
        )
        print("PASS: implementation C certificate regenerates exactly")
        return 0

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered)
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
