#!/usr/bin/env python3
"""Exact integer-interval certificate for the first-prime coercive core.

All interval endpoints are rational multiples of 10**-30.  Arithmetic
rounds outwards with integer floor/ceiling; no floating-point value is
used in the proof.  The two integrals and the analytic quadrature errors
are defined in route report 133, section 16.
"""

from fractions import Fraction
from math import factorial
import json


SCALE = 10**30


def ceildiv(a, b):
    return -((-a) // b)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi

    @staticmethod
    def rat(n, d=1):
        if isinstance(n, Fraction):
            n, d = n.numerator, n.denominator
        return I(n * SCALE // d, ceildiv(n * SCALE, d))

    @staticmethod
    def coerce(v):
        return v if isinstance(v, I) else I.rat(v)

    def __add__(self, other):
        o = I.coerce(other)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -I.coerce(other)

    def __rsub__(self, other):
        return I.coerce(other) + -self

    def __mul__(self, other):
        o = I.coerce(other)
        p = (self.lo * o.lo, self.lo * o.hi,
             self.hi * o.lo, self.hi * o.hi)
        return I(min(p) // SCALE, ceildiv(max(p), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = I.coerce(other)
        assert not (o.lo <= 0 <= o.hi), "interval divisor contains zero"
        if o.hi < 0:
            return (-self) / (-o)
        p = (Fraction(self.lo * SCALE, o.lo),
             Fraction(self.lo * SCALE, o.hi),
             Fraction(self.hi * SCALE, o.lo),
             Fraction(self.hi * SCALE, o.hi))
        low, high = min(p), max(p)
        return I(low.numerator // low.denominator,
                 ceildiv(high.numerator, high.denominator))

    def sq(self):
        if self.lo <= 0 <= self.hi:
            return I(0, ceildiv(max(self.lo**2, self.hi**2), SCALE))
        return self * self

    def widen(self, radius):
        r = I.coerce(radius)
        return I(self.lo - r.hi, self.hi + r.hi)

    def upper(self):
        return Fraction(self.hi, SCALE)

    def lower(self):
        return Fraction(self.lo, SCALE)


def horner_even(x, coefficients):
    x2 = x.sq()
    v = I.rat(coefficients[-1])
    for c in reversed(coefficients[:-1]):
        v = v * x2 + I.rat(c)
    return v


COS_COEFF = [Fraction((-1)**j, factorial(2*j)) for j in range(21)]
SIN_COEFF = [Fraction((-1)**j, factorial(2*j+1)) for j in range(21)]
HCOS_COEFF = [Fraction(1, factorial(2*j)) for j in range(13)]
HSIN_COEFF = [Fraction(1, factorial(2*j+1)) for j in range(13)]
TRIG_ERROR = Fraction(6**42, factorial(42))
HCOS_ERROR = Fraction(2, factorial(26))
HSIN_ERROR = Fraction(2, factorial(27))


def cos(x):
    assert max(abs(x.lo), abs(x.hi)) <= 6 * SCALE
    return horner_even(x, COS_COEFF).widen(TRIG_ERROR)


def sin(x):
    assert max(abs(x.lo), abs(x.hi)) <= 6 * SCALE
    return (x * horner_even(x, SIN_COEFF)).widen(TRIG_ERROR)


def cosh(x):
    assert max(abs(x.lo), abs(x.hi)) <= SCALE
    return horner_even(x, HCOS_COEFF).widen(HCOS_ERROR)


def sinh(x):
    assert max(abs(x.lo), abs(x.hi)) <= SCALE
    return (x * horner_even(x, HSIN_COEFF)).widen(HSIN_ERROR)


R = Fraction(9, 25)
H = 14
D = Fraction(27, 1000)
N = 1024
NX = 16
M = 128
PI_LOWER = Fraction(314159, 100000)
ri = I.rat(R)
cc = ri + sinh(ri)
ss = sinh(ri) - ri
half_sinh = sinh(ri / 2)
half_cosh = cosh(ri / 2)
tail_coefficient = Fraction(1, M**3) + Fraction(1, 2*M**2)


def weight(u):
    """Polynomial/rational majorant b_tilde >= a(H)-a(u), real 0<=u<=H."""
    u2 = u.sq() / 4
    delta = 49 - u2
    out = I.rat(0)
    for n in range(M):
        c = I.rat(4*n+1, 4)
        c2 = c.sq()
        out += c * delta / ((c2 + u2) * (c2 + 49))
    return out + delta * I.rat(tail_coefficient)


def coefficients(u):
    den = u.sq() + I.rat(1, 4)
    sr, cr = sin(u * ri), cos(u * ri)
    c = 2 * (half_sinh * cr / 2 + u * half_cosh * sr) / den
    d = 2 * (half_cosh * sr / 2 - u * half_sinh * cr) / den
    return c, d


def sample(u):
    b = weight(u)
    c, d = coefficients(u)
    residual = 2 * ri - c.sq() / cc - d.sq() / ss
    t = b * residual
    inner = I.rat(0)
    for j in range(NX+1):
        x = I.rat(-R + Fraction(j, NX) * D)
        real = cos(u*x) - c * cosh(x/2) / cc
        imag = -sin(u*x) + d * sinh(x/2) / ss
        z = real.sq() + imag.sq()
        coeff = 1 if j in (0, NX) else (4 if j % 2 else 2)
        inner += coeff * z
    inner *= I.rat(D / (3 * NX))
    e = 2 * b * inner
    return t, e


def as_pair(v):
    return [str(v.numerator), str(v.denominator)]


def check_analytic_constants():
    """Exact elementary premises of the Cauchy/Simpson error bounds."""
    def atan_bounds(x):
        n = 16
        head = sum((-1)**j * x**(2*j+1) / (2*j+1)
                   for j in range(n))
        # An even number of alternating terms is the lower sum.
        return head, head + x**(2*n+1) / (2*n+1)

    l5, u5 = atan_bounds(Fraction(1, 5))
    l239, u239 = atan_bounds(Fraction(1, 239))
    assert 16*l5 - 4*u239 > PI_LOWER  # exact Machin identity
    harmonic = sum(Fraction(1, n) for n in range(1, M))
    assert harmonic < 6
    assert 12 + 2*harmonic + 100*tail_coefficient < 25
    b0 = sum(Fraction(49, 1) /
             (Fraction(4*n+1, 4) *
              (Fraction(4*n+1, 4)**2 + 49)) for n in range(M))
    b0 += 49*tail_coefficient
    assert b0 < 7
    assert Fraction(36, 41) < 1  # 2R exp(2R/4) < 2R/(1-R/2)
    assert (R + Fraction(1, 20))**2 + Fraction(1, 20)**2 < Fraction(21, 50)**2
    x = Fraction(7, 10)
    exp_head = sum(x**j / factorial(j) for j in range(21))
    exp_tail = x**21 / factorial(21) / (1-x/Fraction(22, 1))
    assert exp_head + exp_tail < Fraction(21, 10)
    assert 1 / (1-Fraction(41, 200)**2/2) < Fraction(103, 100)
    assert Fraction(21, 100) / (1-Fraction(21, 100)**2/6) < Fraction(22, 100)
    assert 100*R**3/6 > 2*R
    assert 7*Fraction(533, 100)**2 < 200
    return {"status": "PASS", "analytic_constants": "exact rational checks"}


def main():
    check_analytic_constants()
    # Re psi(1/4+iH/2)-log pi exceeds 3/4.  Ca<27/5 is
    # the previously checked elementary/DLMF bound from 155.12.
    a_lower = -Fraction(27, 5)
    for n in range(M):
        c = Fraction(4*n+1, 4)
        a_lower += Fraction(49, 1) / (c * (c*c + 49))
    assert a_lower > Fraction(3, 4)

    # Log2 > 693/1000 from the positive atanh(1/3) series.
    log2_lower = sum((Fraction(2, 1) /
                      ((2*j+1) * 3**(2*j+1))) for j in range(5))
    assert log2_lower > Fraction(693, 1000)
    assert 2*R - Fraction(693, 1000) == D
    # log2<7/10 and sqrt2>7/5 imply log2/sqrt2<1/2.
    log2_upper = log2_lower + Fraction(2, 1) / (11 * 3**11) * Fraction(9, 8)
    assert log2_upper < Fraction(7, 10)

    total_t, total_e = I.rat(0), I.rat(0)
    for j in range(N+1):
        u = I.rat(Fraction(H*j, N))
        t, e = sample(u)
        coeff = 1 if j in (0, N) else (4 if j % 2 else 2)
        total_t += coeff * t
        total_e += coeff * e
    prefactor = I.rat(Fraction(H, 3*N) / PI_LOWER)
    total_t *= prefactor
    total_e *= prefactor

    # Cauchy radius 1/4, sup<=40: derivative order4<=245760.
    outer_error = Fraction(245760 * H**5, 180 * N**4) / PI_LOWER
    # Inner radius1/20, sup of b_tilde*psi_u*psi_-u <=200:
    # order4<=768000000.  Two edge intervals and u-integral retained.
    inner_error = (Fraction(768000000, 1) * D**5 /
                   (180 * NX**4)) * Fraction(2*H, 1) / PI_LOWER
    t_upper = total_t.upper() + outer_error
    e_upper = total_e.upper() + outer_error + inner_error
    assert t_upper < Fraction(13, 25)
    assert e_upper < Fraction(1, 10)
    trace_upper = Fraction(4, 3)*t_upper + Fraction(8, 3)*e_upper
    assert trace_upper < Fraction(24, 25)
    coercivity = Fraction(1, 4) * (1-trace_upper)
    assert coercivity > Fraction(1, 100)
    result = {
        "status": "PASS",
        "arithmetic": "outward integer intervals; scale=10**30",
        "R": as_pair(R), "H": H, "outer_subintervals": N,
        "inner_subintervals": NX, "digamma_head_terms": M,
        "a_H_lower": as_pair(a_lower),
        "T_upper": as_pair(t_upper), "E_upper": as_pair(e_upper),
        "outer_error": as_pair(outer_error),
        "inner_error": as_pair(inner_error),
        "trace_upper": as_pair(trace_upper),
        "coercivity_lower": as_pair(coercivity),
        "strict_claim": "q_A(f)>=||f||_2**2/100 on both-pole-null cores of diameter<=18/25",
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
