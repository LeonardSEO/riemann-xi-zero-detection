import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Calculus.Deriv.Slope

open Filter Set
open scoped Topology

namespace ExactZeroDetection

/-- The analytic reciprocal takes its removable value at the pole. -/
noncomputable def reciprocalZeta (s : ℂ) : ℂ :=
  if s = 1 then 0 else (riemannZeta s)⁻¹

theorem reciprocalZeta_deriv_one : HasDerivAt reciprocalZeta 1 1 := by
  rw [hasDerivAt_iff_tendsto_slope, slope_fun_def_field]
  have h := riemannZeta_residue_one.inv₀ (one_ne_zero : (1 : ℂ) ≠ 0)
  rw [inv_one] at h
  apply h.congr'
  filter_upwards [self_mem_nhdsWithin] with s hs
  have hs' : s ≠ 1 := hs
  simp [reciprocalZeta, hs', div_eq_mul_inv, mul_comm]

theorem reciprocalZeta_differentiableAt {s : ℂ} (hs : riemannZeta s ≠ 0) :
    DifferentiableAt ℂ reciprocalZeta s := by
  by_cases h : s = 1
  · subst s
    exact reciprocalZeta_deriv_one.differentiableAt
  · apply ((differentiableAt_riemannZeta h).inv hs).congr_of_eventuallyEq
    filter_upwards [isOpen_ne.mem_nhds h] with w hw
    simp [reciprocalZeta, hw]

end ExactZeroDetection
