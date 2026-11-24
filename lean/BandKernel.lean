import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Tactic

namespace ExactZeroDetection

/-- The closed-form multiplier G_Q of (FB3), including its removable value. -/
noncomputable def bandKernel (Q : ℝ) (w : ℂ) : ℂ :=
  if w = 0 then Q ^ 2 / 2 else
    (1 - (1 + Q * w) * Complex.exp (-Q * w)) / w ^ 2

/-- The closed expression agrees with the actual compactly supported moment
integral. At the removable point it is Q²/2. -/
theorem bandKernel_eq_integral (Q : ℝ) (w : ℂ) :
    bandKernel Q w = ∫ t in (0 : ℝ)..Q, (t : ℂ) * Complex.exp (-w * t) := by
  by_cases hw : w = 0
  · subst w
    simp only [bandKernel, if_true, neg_zero, zero_mul, Complex.exp_zero, mul_one]
    rw [intervalIntegral.integral_ofReal, integral_id]
    push_cast
    ring
  let F := fun z : ℂ => -((1 + w*z) * Complex.exp (-w*z)) / w^2
  have hder (z : ℂ) : HasDerivAt F (z * Complex.exp (-w*z)) z := by
    have hz := (hasDerivAt_id z).const_mul w
    have he := ((hasDerivAt_id z).const_mul (-w)).cexp
    convert! (((hz.const_add 1).mul he).neg.div_const (w^2)) using 1 <;>
      first | rfl | (dsimp; field_simp; ring)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt
    (f := fun t : ℝ => F t) (fun t _ => (hder t).comp_ofReal)
    (Continuous.intervalIntegrable (by fun_prop) _ _)]
  simp [F, bandKernel, hw, neg_mul, mul_comm]
  ring

/-- Uniform high-frequency decay of the numerator multiplier. The imaginary
part of w is t+u where t is the Mellin height and |u| ≤ H is the band. -/
theorem bandKernel_tail {Q H t u : ℝ} {w : ℂ} (hQ : 0 ≤ Q)
    (hre : 0 ≤ w.re) (him : w.im = t + u) (hu : |u| ≤ H)
    (ht : 2 * (H + 1) ≤ |t|) :
    ‖bandKernel Q w‖ ≤ 2 * (2 + Q) / |t| := by
  have hH := (abs_nonneg u).trans hu
  have ht' : 0 < |t| := by linarith
  have himle := Complex.abs_im_le_norm w
  rw [him] at himle
  have hdiff := abs_sub (t + u) u
  rw [add_sub_cancel_right] at hdiff
  have htw : |t| ≤ 2 * ‖w‖ := by linarith
  have hw : 1 ≤ ‖w‖ := by linarith
  have hw' : 0 < ‖w‖ := by linarith
  have he : ‖Complex.exp (-(Q : ℂ) * w)‖ ≤ 1 := by
    rw [Complex.norm_exp]
    apply Real.exp_le_one_iff.mpr
    simp [Complex.mul_re, mul_nonneg hQ hre]
  have hn : ‖(1 : ℂ) + Q * w‖ ≤ 1 + Q * ‖w‖ := by
    simpa [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hQ]
      using norm_add_le (1 : ℂ) ((Q : ℂ) * w)
  have hnum : ‖1 - (1 + (Q : ℂ) * w) * Complex.exp (-(Q : ℂ) * w)‖ ≤
      2 + Q * ‖w‖ := by
    calc
      _ ≤ 1 + ‖(1 : ℂ) + Q * w‖ * ‖Complex.exp (-(Q : ℂ) * w)‖ := by
        simpa [norm_mul] using norm_sub_le (1 : ℂ) ((1 + (Q : ℂ) * w) * _)
      _ ≤ 1 + ‖(1 : ℂ) + Q * w‖ * 1 := by gcongr
      _ ≤ _ := by linarith
  rw [bandKernel, if_neg (norm_pos_iff.mp hw'), norm_div, norm_pow]
  calc
    _ ≤ (2 + Q * ‖w‖) / ‖w‖ ^ 2 := by gcongr
    _ ≤ (2 + Q) / ‖w‖ := by
      rw [div_le_div_iff₀ (sq_pos_of_pos hw') hw']
      nlinarith [mul_nonneg (sub_nonneg.mpr hw) hw'.le]
    _ ≤ _ := by
      rw [div_le_div_iff₀ hw' ht']
      nlinarith [mul_nonneg (by linarith : 0 ≤ 2 * ‖w‖ - |t|) (by linarith : 0 ≤ 2 + Q)]

end ExactZeroDetection
