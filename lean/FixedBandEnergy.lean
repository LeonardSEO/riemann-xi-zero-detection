import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.NumberTheory.ArithmeticFunction.Moebius
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith

/-! Endpoint estimates for the fixed-band energy in (FB4). -/

open MeasureTheory Set Filter
open scoped BigOperators

namespace ExactZeroDetection

variable {E : Type*} [NormedAddCommGroup E]

/-- The energy of a projected source as in (FB4), with logarithmic measure. -/
noncomputable def energy (S : ℝ → E) (Y : ℝ) : ℝ :=
  ∫ x in Icc 1 Y, ‖S x‖ ^ 2 / x

/-- The sum in (FB1), with its lower endpoint open. The argument `orbit r`
is the band-projected translate of the convolved input by `r`. -/
noncomputable def mobiusWindow [NormedSpace ℝ E] (q X : ℝ) (orbit : ℝ → E) : E :=
  ∑ n ∈ (Finset.range (⌊X⌋₊ + 1)).filter (fun n : ℕ => X / q < (n : ℝ)),
    (((ArithmeticFunction.moebius n : ℤ) : ℝ) * Real.log (X / n) / Real.sqrt n) •
      orbit (Real.log ((n : ℝ) / X))

/-- A weighted square-integral on the full half-line bounds every endpoint
with the same constant. Only weighted integrability is required. -/
theorem energy_le_of_weighted (S : ℝ → E) {b Y : ℝ} (hb : 0 ≤ b) (hY : 1 ≤ Y)
    (hS : AEStronglyMeasurable S volume)
    (hw : IntegrableOn (fun x : ℝ => x ^ (-(2*b)) * (‖S x‖ ^ 2 / x)) (Ici 1)) :
    energy S Y ≤ Y ^ (2*b) *
      ∫ x in Ici 1, x ^ (-(2*b)) * (‖S x‖ ^ 2 / x) := by
  let w := fun x : ℝ => x ^ (-(2*b)) * (‖S x‖ ^ 2 / x)
  have hsub : Icc 1 Y ⊆ Ici (1 : ℝ) := fun x hx => hx.1
  have hy : 0 ≤ Y ^ (2*b) := Real.rpow_nonneg (by linarith) _
  have hdom (x : ℝ) (hx : x ∈ Icc 1 Y) :
      ‖S x‖ ^ 2 / x ≤ Y ^ (2*b) * w x := by
    have hxpos : 0 < x := lt_of_lt_of_le zero_lt_one hx.1
    have hpow := Real.rpow_le_rpow hxpos.le hx.2 (show 0 ≤ 2*b by linarith)
    dsimp [w]
    rw [← mul_assoc, Real.rpow_neg hxpos.le, ← div_eq_mul_inv]
    exact le_mul_of_one_le_left (by positivity) ((one_le_div (Real.rpow_pos_of_pos hxpos _)).2 hpow)
  have hnon : 0 ≤ᵐ[volume.restrict (Ici (1 : ℝ))] w := by
    filter_upwards [ae_restrict_mem measurableSet_Ici] with x hx
    have : 0 ≤ x := le_trans zero_le_one hx
    dsimp [w]
    positivity
  have hdomae := (ae_restrict_iff' (μ := volume) measurableSet_Icc).2 (Eventually.of_forall hdom)
  have hmaj : IntegrableOn (fun x => Y ^ (2*b) * w x) (Icc 1 Y) :=
    (hw.mono_set hsub).const_mul _
  have hsmall : IntegrableOn (fun x => ‖S x‖ ^ 2 / x) (Icc 1 Y) := by
    apply hmaj.mono' (((hS.norm.pow 2).aemeasurable.div measurable_id.aemeasurable).aestronglyMeasurable).restrict
    filter_upwards [hdomae, ae_restrict_mem measurableSet_Icc] with x hx hmem
    have : 0 ≤ x := le_trans zero_le_one hmem.1
    change ‖‖S x‖ ^ 2 / x‖ ≤ _
    rwa [Real.norm_of_nonneg (div_nonneg (sq_nonneg _) this)]
  calc
    energy S Y ≤ ∫ x in Icc 1 Y, Y ^ (2*b) * w x := integral_mono_ae hsmall hmaj hdomae
    _ = _ := integral_const_mul _ _
    _ ≤ _ := mul_le_mul_of_nonneg_left (setIntegral_mono_set hw hnon
      (Eventually.of_forall hsub)) hy

/-- Square integrability on a vertical Mellin line, from the addendum's
pointwise decay estimate. The strict threshold is η < 1/2. -/
theorem vertical_square_integrable (B : ℝ → E) {η C : ℝ} (hη : η < 1/2)
    (hm : AEStronglyMeasurable B volume)
    (hdecay : ∀ t, ‖B t‖ ≤ C * (1 + ‖t‖) ^ (-1 + η)) :
    Integrable (fun t => ‖B t‖ ^ 2) := by
  have hI : Integrable (fun t : ℝ => (1 + ‖t‖) ^ (-(2 - 2*η))) :=
    integrable_one_add_norm (by simpa using (show (1 : ℝ) < 2 - 2*η by linarith))
  apply (hI.const_mul (C^2)).mono' (hm.norm.pow 2)
  apply Eventually.of_forall
  intro t
  change ‖‖B t‖ ^ 2‖ ≤ _
  rw [Real.norm_of_nonneg (sq_nonneg _)]
  calc
    ‖B t‖ ^ 2 ≤ (C * (1 + ‖t‖) ^ (-1 + η)) ^ 2 :=
      pow_le_pow_left₀ (norm_nonneg _) (hdecay t) 2
    _ = C ^ 2 * (1 + ‖t‖) ^ (-(2 - 2*η)) := by
      rw [mul_pow, sq ((1 + ‖t‖) ^ _), ← Real.rpow_add (by positivity)]
      congr 2; ring

/-- The explicit endpoint estimate for the Möbius window, with an arbitrary
fixed projected orbit. Its analytic input is the weighted L² condition. -/
theorem mobiusWindow_energy_le [NormedSpace ℝ E] (orbit : ℝ → E) (q : ℝ)
    {b Y : ℝ} (hb : 0 ≤ b) (hY : 1 ≤ Y)
    (hm : AEStronglyMeasurable (fun X => mobiusWindow q X orbit) volume)
    (hw : IntegrableOn (fun X : ℝ =>
        X ^ (-(2*b)) * (‖mobiusWindow q X orbit‖ ^ 2 / X)) (Ici 1)) :
    energy (fun X => mobiusWindow q X orbit) Y ≤ Y ^ (2*b) *
      ∫ X in Ici 1, X ^ (-(2*b)) * (‖mobiusWindow q X orbit‖ ^ 2 / X) :=
  energy_le_of_weighted _ hb hY hm hw

end ExactZeroDetection
