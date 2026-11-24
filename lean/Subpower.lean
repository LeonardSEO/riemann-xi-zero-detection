import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Tactic

open Filter
open scoped Topology

namespace ExactZeroDetection

/-- Uniform interpolation loss for all interior points separated from the
outer circle. This is the quantitative reason the strict margin b-a matters. -/
theorem log_interpolation_loss {A C kmax η : ℝ} (hA : 1 ≤ A) (hC : 1 ≤ C)
    (hkmax : kmax < 1) (hη : 0 < η) :
    ∀ᶠ L : ℝ in atTop, ∀ k ∈ Set.Icc 0 kmax,
      A ^ (1-k) * (C*L) ^ k ≤ η*L := by
  have hK : 0 < A*C := mul_pos (zero_lt_one.trans_le hA) (zero_lt_one.trans_le hC)
  have hlim := (tendsto_rpow_neg_atTop (show 0 < 1-kmax by linarith)).eventually
    (gt_mem_nhds (div_pos hη hK))
  filter_upwards [hlim, eventually_ge_atTop (1 : ℝ)] with L hL h1
  intro k hk
  have hp : 0 < L := zero_lt_one.trans_le h1
  have ha : A ^ (1-k) ≤ A := by
    simpa using Real.rpow_le_rpow_of_exponent_le hA (show 1-k ≤ 1 by linarith [hk.1])
  have hc : C ^ k ≤ C := by
    simpa using Real.rpow_le_rpow_of_exponent_le hC (hk.2.trans hkmax.le)
  have hpow : L ^ k ≤ L ^ kmax := Real.rpow_le_rpow_of_exponent_le h1 hk.2
  have hineq : L ^ (-(1-kmax)) * (A*C) ≤ η := (le_div_iff₀ hK).mp hL.le
  rw [Real.mul_rpow (by linarith : 0 ≤ C) hp.le]
  calc
    _ ≤ A * (C * L ^ kmax) := by gcongr
    _ = (L ^ (-(1-kmax)) * (A*C)) * L := by
      have he : L ^ (-(1-kmax)) * L = L ^ kmax := by
        conv_lhs => rhs; rw [← Real.rpow_one L]
        rw [← Real.rpow_add hp]; congr 1; ring
      rw [← he]; ring
    _ ≤ _ := mul_le_mul_of_nonneg_right hineq hp.le

end ExactZeroDetection
