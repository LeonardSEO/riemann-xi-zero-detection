# Lean checks for the fixed-band energy addendum

These modules accompany the [24 November addendum](../paper/addendum.pdf). The project pins Lean 4.33.1 and the Mathlib revision used in the [Alpöge proof mirror](https://github.com/CocoSgt/riemann-hypothesis-lean/tree/7cbf7a0). The declarations below use Mathlib and do not invoke RH.

The full 265-page explanatory paper is [Rubinstein, version 2, 22 November 2026](https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf). Its abstract identifies the main inequality (M) as an input from the Lean development. Part II, §1, pp. 30–31 proves the one-sided interval lemma. Take θ = 1/2 + a. Its windows are [Z, Z^(3/2)] with arbitrarily large Z.

With the pinned toolchain available, run these commands in this directory.

~~~~sh
lake exe cache get Mathlib
lake build
lake env lean Axioms.lean
~~~~

| Addendum step | Lean declaration and scope |
| --- | --- |
| The pole at one | `reciprocalZeta_deriv_one` gives the derivative of the extended reciprocal at one. `reciprocalZeta_differentiableAt` treats the rest of the zero-free region. |
| Disk interpolation | `reciprocal_bound_on_disk` applies three-circle interpolation to a local logarithm of zeta whose two boundary bounds are supplied. |
| The tested contour shift | `tested_contour_shift` equates the two vertical integrals from differentiability on a strip, integrability on the two lines, and uniform decay across horizontal segments. |
| The Möbius factor in (FB2) | `mellin_mobiusWindow` takes the Mellin transform of the finite arithmetic sum for a continuous projected orbit. On Re z > 1/2 this equals the transform of its compact window kernel divided by ζ(1/2+z). |
| The numerator in (FB3) | `bandKernel_tail` bounds the explicit scalar multiplier for all band frequencies at sufficiently large Mellin height. |
| Integrating along a vertical line | `vertical_square_integrable` derives square integrability from measurable boundary data with norm decay of order −1+η for η<1/2. |
| Returning to the source | `weighted_of_inverseFourier` transfers square integrability to the weighted arithmetic variable, conditional on the almost-everywhere identification of the shifted source with the inverse Fourier transform. Its Fourier variable is in cycles per unit. |
| Every endpoint in (FB4) | `energy_le_of_weighted` and `mobiusWindow_energy_le` give the all-endpoint bound with the weighted integral as constant. The orbit parameter represents the projected translates of the convolved input. The sum uses the actual Möbius coefficients and the half-open integer window. |

**Formal status.** An end-to-end declaration from the prime-power windows to the paper's energy bound is still open. In particular, the project needs the specialization of the external interval lemma, the reciprocal-zeta vertical growth estimate, the concrete theta-convolution and projection identification in (FB2)–(FB3), and the contour shift for that source. The hypothesis of `weighted_of_inverseFourier` records exactly where the contour argument enters. The October detector (FB6) is a separate theorem in the manuscript. The axiom audit concerns the declarations listed in `Axioms.lean`; hypotheses in their types remain part of each result.
