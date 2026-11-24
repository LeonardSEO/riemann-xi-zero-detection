import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

open Complex MeasureTheory Set Filter
open scoped Topology Interval

namespace ExactZeroDetection

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℂ E]

/-- Vertical contour shift for a tested Mellin transform. The horizontal
majorant is uniform across the strip, and need only vanish at infinity.
The endpoint line integrals are Bochner integrals. -/
theorem tested_contour_shift (f : ℂ → E) {a b : ℝ} (hab : a ≤ b)
    (hd : DifferentiableOn ℂ f {z | z.re ∈ Icc a b})
    (ha : Integrable (fun t : ℝ => f (a + t * I)))
    (hb : Integrable (fun t : ℝ => f (b + t * I)))
    (m : ℝ → ℝ) (hm : Tendsto m atTop (𝓝 0))
    (hbound : ∀ t : ℝ, ∀ x ∈ Icc a b, ‖f (x + t * I)‖ ≤ m |t|) :
    (∫ t : ℝ, f (a + t * I)) = ∫ t : ℝ, f (b + t * I) := by
  let h := fun t : ℝ => ∫ x in a..b, f (x + t * I)
  have hh (t : ℝ) : ‖h t‖ ≤ m |t| * |b-a| := by
    apply intervalIntegral.norm_integral_le_of_norm_le_const
    intro x hx
    rw [uIoc_of_le hab] at hx
    exact hbound t x ⟨hx.1.le,hx.2⟩
  have hlim (j : ℝ → ℝ) (hj : Tendsto (fun T => |j T|) atTop atTop) :
      Tendsto (fun T => h (j T)) atTop (𝓝 0) := by
    apply squeeze_zero_norm (fun T => hh (j T))
    simpa using (hm.comp hj).mul_const |b-a|
  have hp := hlim id (by simpa using (tendsto_abs_atTop_atTop : Tendsto abs atTop atTop))
  have hn := hlim Neg.neg (by simpa using (tendsto_abs_atTop_atTop : Tendsto abs atTop atTop))
  let v := fun (x T : ℝ) => ∫ t in -T..T, f (x + t * I)
  have hav := intervalIntegral_tendsto_integral ha tendsto_neg_atTop_atBot tendsto_id
  have hbv := intervalIntegral_tendsto_integral hb tendsto_neg_atTop_atBot tendsto_id
  have side (T : ℝ) : h (-T) - h T + I • v b T - I • v a T = 0 := by
    have H := integral_boundary_rect_eq_zero_of_differentiableOn f
      (⟨a,-T⟩ : ℂ) (⟨b,T⟩ : ℂ)
      (hd.mono (by intro z hz; simpa [uIcc_of_le hab] using hz.1))
    exact H
  have hzero := tendsto_nhds_unique
    (((hn.sub hp).add (hbv.const_smul I)).sub (hav.const_smul I))
    (tendsto_const_nhds.congr (fun T => (side T).symm))
  simp only [sub_self, zero_add, sub_eq_zero] at hzero
  exact (smul_right_injective E I_ne_zero hzero).symm

end ExactZeroDetection
