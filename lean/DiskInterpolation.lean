import Reciprocal
import Mathlib.Analysis.Complex.Hadamard
import Mathlib.Analysis.SpecialFunctions.Complex.LogDeriv

open Set Complex Complex.HadamardThreeLines

namespace ExactZeroDetection

/-- Three-circle interpolation in the translated disk coordinates used for
a branch of log zeta. l and r are the logarithms of the two radii. -/
theorem norm_le_three_circles {f : ℂ → ℂ} {l r A B : ℝ} (hlr : l < r)
    (hf : DifferentiableOn ℂ f (Metric.closedBall 0 (Real.exp r)))
    (hinner : ∀ z : ℂ, ‖z‖ = Real.exp l → ‖f z‖ ≤ A)
    (houter : ∀ z : ℂ, ‖z‖ = Real.exp r → ‖f z‖ ≤ B)
    {w : ℂ} (hw : Real.exp l ≤ ‖w‖) (hw' : ‖w‖ ≤ Real.exp r) :
    ‖f w‖ ≤ A ^ (1 - (Real.log ‖w‖ - l) / (r-l)) *
               B ^ ((Real.log ‖w‖ - l) / (r-l)) := by
  let D := Metric.closedBall (0 : ℂ) (Real.exp r)
  have hmaps : MapsTo Complex.exp (verticalClosedStrip l r) D := by
    intro s hs
    change dist (Complex.exp s) 0 ≤ _
    rw [dist_zero_right, Complex.norm_exp]
    exact Real.exp_le_exp.mpr hs.2
  have hc : closure (verticalStrip l r) = verticalClosedStrip l r := by
    rw [verticalStrip, verticalClosedStrip, closure_preimage_re, closure_Ioo hlr.ne]
  have hd : DiffContOnCl ℂ (f ∘ Complex.exp) (verticalStrip l r) := by
    constructor
    · apply hf.comp Complex.differentiable_exp.differentiableOn
      exact hmaps.mono_left (fun x hx => ⟨hx.1.le,hx.2.le⟩)
    · rw [hc]
      exact hf.continuousOn.comp Complex.continuous_exp.continuousOn hmaps
  have hbounded : BddAbove ((norm ∘ (f ∘ Complex.exp)) '' verticalClosedStrip l r) := by
    exact ((isCompact_closedBall (0 : ℂ) _).bddAbove_image hf.continuousOn.norm).mono
      (by rintro _ ⟨z,hz,rfl⟩; exact ⟨_,hmaps hz,rfl⟩)
  have hw0 : w ≠ 0 := norm_pos_iff.mp ((Real.exp_pos _).trans_le hw)
  have hlog : Complex.log w ∈ verticalClosedStrip l r := by
    change l ≤ (Complex.log w).re ∧ (Complex.log w).re ≤ r
    rw [Complex.log_re]
    exact ⟨(Real.le_log_iff_exp_le (norm_pos_iff.mpr hw0)).2 hw,
      (Real.log_le_iff_le_exp (norm_pos_iff.mpr hw0)).2 hw'⟩
  have bound := norm_le_interp_of_mem_verticalClosedStrip' hlr hlog hd hbounded
    (a := A) (b := B)
    (fun z hz => hinner _ (by rw [Complex.norm_exp, show z.re = l from hz]))
    (fun z hz => houter _ (by rw [Complex.norm_exp, show z.re = r from hz]))
  simpa [Function.comp_def, Complex.exp_log hw0, Complex.log_re] using bound

/-- Application to a local logarithm of zeta on a translated disk. Bounding
the branch on the two circles yields a reciprocal bound at the interior point. -/
theorem reciprocal_bound_on_disk {L : ℂ → ℂ} {c w : ℂ} {l r A B : ℝ}
    (hlr : l < r) (hd : DifferentiableOn ℂ L (Metric.closedBall 0 (Real.exp r)))
    (hi : ∀ z : ℂ, ‖z‖ = Real.exp l → ‖L z‖ ≤ A)
    (ho : ∀ z : ℂ, ‖z‖ = Real.exp r → ‖L z‖ ≤ B)
    (hw : Real.exp l ≤ ‖w‖) (hw' : ‖w‖ ≤ Real.exp r)
    (hlog : Complex.exp (L w) = riemannZeta (c+w)) :
    ‖reciprocalZeta (c+w)‖ ≤ Real.exp
      (A ^ (1 - (Real.log ‖w‖ - l) / (r-l)) * B ^ ((Real.log ‖w‖ - l) / (r-l))) := by
  by_cases h : c + w = 1
  · simp [reciprocalZeta, h, (Real.exp_pos _).le]
  rw [reciprocalZeta, if_neg h, ← hlog, ← Complex.exp_neg, Complex.norm_exp]
  apply Real.exp_le_exp.mpr
  exact ((neg_le_abs (L w).re).trans (Complex.abs_re_le_norm _)).trans
    (norm_le_three_circles hlr hd hi ho hw hw')

end ExactZeroDetection
