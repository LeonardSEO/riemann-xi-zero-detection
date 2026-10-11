# Finite certificates and computational support

These programs and frozen outputs support the specifically identified finite results in the paper. The general analytic theorems require no numerical certificate. None of these files certifies RH or an off-line zeta zero.

The source bytes and output bytes are preserved. `SOURCE_HASHES.json` records their SHA-256 values. Historical file names containing `attempt`, `candidate` or a numeric prefix are retained for integrity; their mathematical scope is given below and in the manuscript.

| Files | Certified or supporting scope | Replay |
|---|---|---|
| `xi-moments/theta_moment_certificate.py` and `.out` | Degree-10 endpoint moment inequalities and their finite root brackets. Upstream program retains AGPL-3.0. | `python3 certificates/xi-moments/theta_moment_certificate.py` |
| `finite-interlacing/` | Two separate interval programs for the finite model at lambda squared 6, N=2; every prime power 2,3,4,5 is retained. Frozen specification included. | Run either program in its `code/` directory, in a disposable copy. |
| `n1_tail_bound.py`, `n1_contour_attempt.py` and outputs | First theta head: rigorous infinite tail, zero-free dyadic rectangle edges, winding count two, and real-axis signs. | See the explicit-output commands in REPRODUCIBILITY.md. |
| `head3_local_rouche.py` and output | One simple non-real zero of the third finite theta head in its exact dyadic disk. | Explicit-output replay. |
| `heated_head3_local_rouche.py` and output | One simple non-real zero of the separately heated third finite head. | Explicit-output replay. |
| `148_mehler_*_audit.py` and outputs | A genuine modular Mehler approximant's non-real zero, and a separate elementary constant for the full-strip analytic error. | Run in a disposable certificate-directory copy. |
| `133_first_prime_core_check.py`, `133_two_prime_core_check.py`, `133_three_prime_core_check.py`, `133_prime_five_core_check.py` and outputs | Rational interval positivity cores with the full stated high-sector and Schur budgets; retained prime-power shifts through the indicated support diameter. | Standard-library programs; preserve sibling imports. |
| `155_small_core_rational_check.py`, `155_two_cell_prime_edge_check.py` | Exact rational assertions for the small gamma/pole core and the two-cell prime-edge constructions. | Standard-library programs. |
| `218_*_birth.py`, `218_adaptive_mehler_double_node.py` and outputs | Explicit theta double nodes and the individually stated non-real birth examples. | Arb replay in a disposable copy. |
| `296_theta_recursion_arb.py` and output | Positive full-theta Möbius cross increment with the complete physical tail paid. | `python3 certificates/296_theta_recursion_arb.py` |
| `348_erf_phase_sign.py` and output | The stated coherent erf-partner phase-sign example with rigorous ball errors. | Arb replay in a disposable copy. |
| `132_infinite_logconcavity_nonreal_model_check.py` | Exact rational checks for a model obstruction to unrestricted log-concavity implications. | Standard-library program. |

The separate unclosed diameter-19/10 exploratory core and its walk-filter files are not used as closed theorem inputs and are not part of this release. The third-head and Mehler outputs demonstrate limits of approximating constructions, not zeros of xi or zeta.

The numerical statements in the manuscript provide the actual contour or disk geometry, precision and analytic tail arguments. The files enable reproduction of their finite evidence. Original runtime metadata is retained; see REPRODUCIBILITY.md for the upstream metadata-comparison limitation.
