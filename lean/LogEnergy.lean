import FixedBandEnergy
import Mathlib.MeasureTheory.Function.JacobianOneDim
import Mathlib.Analysis.Fourier.LpSpace

open MeasureTheory Set
open scoped FourierTransform

namespace ExactZeroDetection

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]
  [CompleteSpace E]

/-- If the shifted logarithmic source agrees almost everywhere on v ≥ 0
with the inverse L² Fourier transform of the boundary data, then its
weighted energy is finite. Frequencies supplied to Mathlib are in cycles
per unit, i.e. the boundary datum here is B(b+2πit). -/
theorem weighted_of_inverseFourier (S : ℝ → E) (b : ℝ)
    (boundary : Lp (α := ℝ) E 2)
    (hshift : (fun v : ℝ => Real.exp (-b*v) • S (Real.exp v))
      =ᵐ[volume.restrict (Ici 0)] (fun v : ℝ => (𝓕⁻ boundary : Lp (α := ℝ) E 2) v)) :
    IntegrableOn (fun x : ℝ => x ^ (-(2*b)) * (‖S x‖ ^ 2 / x)) (Ici 1) := by
  have hsq := (Lp.memLp (𝓕⁻ boundary)).integrable_norm_pow (p := 2) (by decide)
  have hlog : IntegrableOn (fun v : ℝ => ‖Real.exp (-b*v) • S (Real.exp v)‖ ^ 2)
      (Ici 0) := by
    apply hsq.integrableOn.congr
    filter_upwards [hshift] with v hv using congrArg (fun e : E => ‖e‖ ^ 2) hv.symm
  rw [← Real.exp_zero, ← Real.image_exp_Ici]
  apply (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ici
    (fun v _ => (Real.hasDerivAt_exp v).hasDerivWithinAt) Real.exp_injective.injOn _).2
  apply hlog.congr_fun _ measurableSet_Ici
  intro v _
  simp only [norm_smul, Real.norm_eq_abs, abs_of_pos (Real.exp_pos _), smul_eq_mul,
    mul_pow, ← Real.exp_mul, ← Real.exp_nat_mul]
  norm_num only [Nat.cast_ofNat]
  rw [show (2 : ℝ) * (-b*v) = v * -(2*b) by ring]
  field_simp

end ExactZeroDetection
