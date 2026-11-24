import FixedBandEnergy
import Mathlib.Analysis.MellinTransform
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Tactic

open MeasureTheory Set Complex

namespace ExactZeroDetection

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℂ E]

/-- The L¹ cost of a Mellin dilation, in exactly the same half-plane as the
complex transform. This controls the exchange of sum and integral. -/
private theorem mellin_norm_dilation (f : ℝ → E) (s : ℂ) {d : ℝ} (hd : 0 < d) :
    (∫ x in Ioi (0:ℝ), ‖(x:ℂ) ^ (s-1) • f (d*x)‖) =
      d ^ (-s.re) * ∫ x in Ioi (0:ℝ), ‖(x:ℂ) ^ (s-1) • f x‖ := by
  have e (x : ℝ) (hx : x ∈ Ioi (0:ℝ)) :
      ‖(x:ℂ) ^ (s-1) • f (d*x)‖ =
        d ^ (1-s.re) * ‖((d*x:ℝ):ℂ) ^ (s-1) • f (d*x)‖ := by
    rw [norm_smul, norm_smul, norm_cpow_eq_rpow_re_of_pos hx,
      norm_cpow_eq_rpow_re_of_pos (mul_pos hd hx), sub_re, one_re, Real.mul_rpow hd.le hx.le]
    rw [← mul_assoc, ← mul_assoc, ← Real.rpow_add hd]
    simp
  rw [setIntegral_congr_fun measurableSet_Ioi e, integral_const_mul,
    integral_comp_mul_left_Ioi (fun x : ℝ => ‖(x:ℂ) ^ (s-1) • f x‖) _ hd]
  simp only [mul_zero, smul_eq_mul]
  rw [← mul_assoc, ← Real.rpow_neg_one, ← Real.rpow_add hd]
  ring_nf

/-- For a summable weighted series of dilations, the Mellin transform
factors into a Dirichlet coefficient and the transform of a single kernel. -/
theorem hasSum_mellin_dilations (f : ℝ → E) (s : ℂ)
    (c : ℕ → ℂ) (d : ℕ → ℝ) (hd : ∀ n, 0 < d n) (hf : MellinConvergent f s)
    (hcost : Summable (fun n => ‖c n‖ * (d n) ^ (-s.re))) :
    HasSum (fun n => (c n * (d n : ℂ)^(-s)) • mellin f s)
      (mellin (fun x => ∑' n, c n • f (d n * x)) s) := by
  let V := fun n (x : ℝ) => (x:ℂ)^(s-1) • (c n • f (d n*x))
  have hI (n : ℕ) : IntegrableOn (V n) (Ioi (0:ℝ)) := by
    have h := ((MellinConvergent.comp_mul_left (hd n)).2 hf).const_smul (c n)
    exact h
  have hN (n : ℕ) : (∫ x in Ioi (0:ℝ), ‖V n x‖) =
      (‖c n‖ * (d n)^(-s.re)) * ∫ x in Ioi (0:ℝ), ‖(x:ℂ)^(s-1) • f x‖ := by
    simp only [V, smul_comm _ (c n), norm_smul (c n), integral_const_mul,
      mellin_norm_dilation f s (hd n), mul_assoc]
  have h := hasSum_integral_of_summable_integral_norm hI
    (by simpa only [hN] using hcost.mul_right (∫ x in Ioi (0:ℝ), ‖(x:ℂ)^(s-1) • f x‖))
  have hval (n : ℕ) : (∫ x in Ioi (0:ℝ), V n x) =
      (c n * (d n : ℂ)^(-s)) • mellin f s := by
    change mellin (fun x => c n • f (d n*x)) s = _
    rw [mellin_const_smul (fun x => f (d n*x)), mellin_comp_mul_left _ _ (hd n), mul_smul]
  simpa only [hval, V, mellin, tsum_const_smul''] using h

private noncomputable def invIndex (n : ℕ) : ℝ := (max 1 (n:ℝ))⁻¹
private noncomputable def muCoeff (n : ℕ) : ℂ :=
  LSeries.term (fun k => ((ArithmeticFunction.moebius k : ℤ) : ℂ)) (1/2) n

/-- The Möbius dilation sum for a kernel k supported in [1,q). The zero
index has coefficient zero. For n ≥ 1 the dilation is X/n. -/
noncomputable def mobiusDilations (k : ℝ → E) (X : ℝ) : E :=
  ∑' n : ℕ, muCoeff n • k (invIndex n * X)

private theorem coeff_identity (z : ℂ) (n : ℕ) :
    muCoeff n * (invIndex n : ℂ)^z =
      LSeries.term (fun k => ((ArithmeticFunction.moebius k : ℤ) : ℂ)) (1/2+z) n := by
  by_cases hn : n = 0
  · simp [hn,muCoeff]
  have hp : 1 ≤ (n:ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr hn
  have hcn : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn
  simp only [muCoeff, LSeries.term_of_ne_zero hn, invIndex, max_eq_right hp,
    Complex.ofReal_inv, Complex.inv_cpow _ z (by simp [Real.pi_ne_zero.symm] : ((n : ℂ).arg) ≠ Real.pi),
    Complex.ofReal_natCast, Complex.cpow_add _ _ hcn, div_eq_mul_inv, mul_inv, mul_assoc]

/-- The Möbius reciprocal factor of (FB2), as an identity of actual Bochner
integrals on its original half-plane of absolute convergence. -/
theorem mellin_mobiusDilations (k : ℝ → E) {z : ℂ} (hz : 1/2 < z.re)
    (hk : MellinConvergent k (-z)) :
    mellin (mobiusDilations k) (-z) = (riemannZeta (1/2+z))⁻¹ • mellin k (-z) := by
  let μc : ℕ → ℂ := fun n => ((ArithmeticFunction.moebius n : ℤ) : ℂ)
  have hs : 1 < ((1/2:ℂ)+z).re := by simp; linarith
  have hμ : Summable (LSeries.term μc (1/2+z)) :=
    ArithmeticFunction.LSeriesSummable_moebius_iff.mpr hs
  have hd (n : ℕ) : 0 < invIndex n := inv_pos.mpr (lt_of_lt_of_le zero_lt_one (le_max_left ..))
  have cost : Summable (fun n => ‖muCoeff n‖ * invIndex n ^ (-(-z).re)) := by
    refine hμ.norm.congr fun n => ?_
    rw [← coeff_identity, norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos (hd n)]
    simp
  have h := hasSum_mellin_dilations k (-z) muCoeff invIndex hd hk cost
  simp only [neg_neg,coeff_identity] at h
  have hi := LSeries_one_mul_Lseries_moebius hs
  rw [LSeries_one_eq_riemannZeta hs] at hi
  have hsum : (∑' n, LSeries.term μc (1/2+z) n) = (riemannZeta (1/2+z))⁻¹ :=
    eq_inv_of_mul_eq_one_right hi
  exact h.unique (hsum ▸ hμ.hasSum.smul_const (mellin k (-z)))

/-- The compact kernel whose dilations cut out exactly X/q < n ≤ X. -/
noncomputable def windowSeed (q : ℝ) (orbit : ℝ → E) (x : ℝ) : E :=
  if x ∈ Ico 1 q then Real.log x • orbit (-Real.log x) else 0

private theorem dilation_term (q X : ℝ) (hq : 0 < q) (hX : 0 < X) (orbit : ℝ → E) (n : ℕ) :
    muCoeff n • windowSeed q orbit (invIndex n * X) =
      if (n : ℝ) ≤ X ∧ X/q < n then
        (((ArithmeticFunction.moebius n : ℤ) : ℝ) * Real.log (X/n) / Real.sqrt n) •
          orbit (Real.log ((n : ℝ)/X)) else 0 := by
  by_cases hn : n = 0
  · subst n
    simp [muCoeff, not_lt_of_ge (div_nonneg hX.le hq.le)]
  have h1 : 1 ≤ (n:ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr hn
  have hn' : 0 < (n:ℝ) := zero_lt_one.trans_le h1
  have hi : invIndex n * X = X/n := by simp [invIndex,max_eq_right h1,div_eq_inv_mul]
  have hmem : X/n ∈ Ico 1 q ↔ (n:ℝ) ≤ X ∧ X/q < n := by
    rw [Set.mem_Ico, one_le_div hn', div_lt_iff₀ hn', div_lt_iff₀ hq, mul_comm q]
  rw [hi,windowSeed,if_congr hmem rfl rfl]
  split_ifs with h
  · have hcoeff : muCoeff n = ((((ArithmeticFunction.moebius n : ℤ) : ℝ) / Real.sqrt n : ℝ) : ℂ) := by
      simp only [muCoeff,LSeries.term_of_ne_zero hn, Real.sqrt_eq_rpow, Complex.ofReal_div,
        Complex.ofReal_cpow hn'.le, Complex.ofReal_natCast, Complex.ofReal_intCast,
        Complex.ofReal_one, Complex.ofReal_ofNat]
    rw [hcoeff, Complex.coe_smul, smul_smul, ← Real.log_inv, inv_div]
    congr 1
    ring
  · simp

theorem mobiusDilations_windowSeed_eq (q X : ℝ) (hq : 0 < q) (hX : 0 < X) (orbit : ℝ → E) :
    mobiusDilations (windowSeed q orbit) X = mobiusWindow q X orbit := by
  unfold mobiusDilations
  simp_rw [dilation_term q X hq hX orbit]
  let S := Finset.range (⌊X⌋₊+1)
  have hmem (n : ℕ) : n ∈ S ↔ (n:ℝ) ≤ X := by
    simp only [S, Finset.mem_range, Nat.lt_succ_iff, Nat.le_floor_iff hX.le]
  rw [tsum_eq_sum (s := S) (fun n hn => if_neg (fun h => hn ((hmem n).2 h.1)))]
  unfold mobiusWindow
  rw [Finset.sum_filter]
  exact Finset.sum_congr rfl (fun n hn => by simp [(hmem n).1 hn])


/-- A continuous projected orbit gives an absolutely convergent compact
window on every Mellin line. -/
theorem windowSeed_convergent {orbit : ℝ → E} (ho : Continuous orbit) (q : ℝ) (s : ℂ) :
    MellinConvergent (windowSeed q orbit) s := by
  let f := fun x : ℝ => (x : ℂ) ^ (s-1) • (Real.log x • orbit (-Real.log x))
  have hc : ContinuousOn f (Icc 1 q) := by
    intro x hx
    have hx' : x ≠ 0 := ne_of_gt (zero_lt_one.trans_le hx.1)
    have hl := Real.continuousAt_log hx'
    exact ((Complex.continuousAt_ofReal_cpow_const x (s-1) (Or.inr hx')).smul
      (hl.smul (ho.continuousAt.comp hl.neg))).continuousWithinAt
  have hi : Integrable ((Ico 1 q).indicator f) :=
    (hc.integrableOn_Icc.mono_set Ico_subset_Icc_self).integrable_indicator measurableSet_Ico
  apply hi.integrableOn.congr_fun _ measurableSet_Ioi
  intro x _
  by_cases hx : x ∈ Ico 1 q <;> simp [f,windowSeed,hx]

/-- The Mellin identity on the original line for the same finite source
as in the endpoint-energy estimate. -/
theorem mellin_mobiusWindow (orbit : ℝ → E) {q : ℝ} (hq : 0 < q) {z : ℂ}
    (hz : 1/2 < z.re) (ho : Continuous orbit) :
    mellin (fun X => mobiusWindow q X orbit) (-z) =
      (riemannZeta (1/2+z))⁻¹ • mellin (windowSeed q orbit) (-z) := by
  rw [← mellin_mobiusDilations _ hz (windowSeed_convergent ho q (-z))]
  apply setIntegral_congr_fun measurableSet_Ioi
  intro x hx
  exact congrArg (fun v : E => (x : ℂ)^(-z-1) • v) (mobiusDilations_windowSeed_eq q x hq hx orbit).symm

end ExactZeroDetection
