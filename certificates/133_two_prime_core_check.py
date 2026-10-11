#!/usr/bin/env python3
"""Rational-interval polynomial/LDL certificate for a two-prime core.

There is no quadrature and no floating-point arithmetic.  Every piecewise
polynomial integral is evaluated exactly through Legendre coefficients.
The mathematical high-sector reduction is in route133, section18.
"""

import argparse
from fractions import Fraction as F
from math import factorial, comb, isqrt
from pathlib import Path
import importlib.util
import json

spec = importlib.util.spec_from_file_location(
    "certified_intervals", Path(__file__).with_name("133_first_prime_core_check.py"))
interval_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(interval_module)
interval_module.SCALE = 10**150
I = interval_module.I
SCALE = interval_module.SCALE


def idiv(self, other):
    o = I.coerce(other)
    assert not (o.lo <= 0 <= o.hi)
    if o.hi < 0:
        return (-self) / (-o)
    pairs = ((self.lo, o.lo), (self.lo, o.hi),
             (self.hi, o.lo), (self.hi, o.hi))
    lows = [a*SCALE//b for a, b in pairs]
    highs = [-((-a*SCALE)//b) for a, b in pairs]
    return I(min(lows), max(highs))


I.__truediv__ = idiv


def ipow(x, n):
    r = I.rat(1)
    while n:
        if n % 2:
            r *= x
        x = x.sq()
        n //= 2
    return r


R = F(11, 20)
N = 40
P = 16
J = 64
DEG = 24
PENALTY = 100000
SCHUR_FLOOR = F(1, 2500)
KERNEL_ERROR = 7*F(11, 20)**(DEG+1)/(1-F(11, 20))
zero = I.rat(0)
one = I.rat(1)


def atanh_log(q, terms=200):
    q = I.coerce(q)
    q2 = q.sq()
    v = q
    out = zero
    for n in range(terms):
        out += 2*v/I.rat(2*n+1)
        v *= q2
    error = 2*v/I.rat(2*terms+1)/(one-q2)
    return I(out.lo, out.hi+error.hi)


def atan(q, terms=200):
    q = F(q)
    v = sum((-1)**n*q**(2*n+1)/F(2*n+1) for n in range(terms))
    e = q**(2*terms+1)/F(2*terms+1)
    return I.rat(v).widen(e)


def bernoulli(maxn):
    b = [F(1)]
    for n in range(1, maxn+1):
        b.append(-sum(F(comb(n+1, k))*b[k] for k in range(n))/F(n+1))
    return b


def constants():
    pi = 16*atan(F(1, 5))-4*atan(F(1, 239))
    log2 = atanh_log(I.rat(1, 3))
    log3 = atanh_log(I.rat(1, 2))
    logpi = atanh_log((pi-one)/(pi+one))
    m = 32
    b = bernoulli(18)
    gamma = I.rat(sum(F(1, n) for n in range(1, m+1))-F(1, 2*m))
    gamma -= 5*log2
    gamma += I.rat(sum(b[2*k]/F(2*k*m**(2*k)) for k in range(1, 9)))
    # DLMF5.11(ii): positive-real digamma remainder, first neglected term.
    gamma = gamma.widen(abs(b[18])/F(18*m**18))
    a0 = -(gamma+pi/2+3*log2+logpi)
    roots = []
    for n in (2, 3):
        r = isqrt(n*SCALE*SCALE)
        roots.append(I(r, r+1))
    return a0, log2, log3, log2/roots[0], log3/roots[1]


def exp_negative(z):
    q = I.coerce(z)/128
    assert q.lo >= 0 and q.hi <= SCALE
    out = I.rat(F((-1)**128, factorial(128)))
    for k in range(127, -1, -1):
        out = out*q+I.rat(F((-1)**k, factorial(k)))
    out = out.widen(F(1, factorial(129)))
    assert out.lo > 0
    for unused in range(7):
        out = out.sq()
    return out


def fadd(a, b):
    out = [F(0)]*max(len(a), len(b))
    for j, v in enumerate(a): out[j] += v
    for j, v in enumerate(b): out[j] += v
    return out


def fmul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for j, v in enumerate(a):
        if v:
            for k, w in enumerate(b):
                if w: out[j+k] += v*w
    return out


def legendre_polynomials(n):
    ps = [[F(1)], [F(0), F(1)]]
    for k in range(1, n-1):
        a = [F(0)]+[F(2*k+1, k+1)*v for v in ps[-1]]
        b = [-F(k, k+1)*v for v in ps[-2]]
        ps.append(fadd(a, b))
    return ps[:n]


def kernel_polynomial():
    maxn = DEG+1
    u = [F(1, factorial(n+1)) if not n % 2 else F(0)
         for n in range(maxn+1)]  # sinh(r)/r
    inv = [F(1)]
    for n in range(1, maxn+1):
        inv.append(-sum(u[j]*inv[n-j] for j in range(1, n+1)))
    e = [F(1, 2**n*factorial(n)) for n in range(maxn+1)]
    coeff = [sum(e[j]*inv[n+1-j] for j in range(n+2))/2
             for n in range(DEG+1)]
    coeff[0] -= KERNEL_ERROR
    oscillation = sum(abs(coeff[n])*(2*R)**n for n in range(1, DEG+1))
    assert coeff[0]-oscillation > F(1, 6)
    assert coeff[0]+oscillation < F(7, 20)
    return coeff


def kernel_action(poly, coeff):
    out = [F(0)]*(len(poly)+DEG+1)
    for j, aj in enumerate(coeff):
        # Integral_{-1}^1 |t-s|^j ds.
        row = [F(comb(j+1, m)*(1+(-1)**m), j+1)
               for m in range(j+2)]
        conv = [F(0)]*(len(poly)+j+1)
        for k, pk in enumerate(poly):
            if not pk: continue
            for m in range(j+1):
                conv[j-m] += pk*F(comb(j, m)*
                      (-(-1)**(k+1)+(-1)**(j-m)), k+m+1)
            if j % 2:
                conv[j+k+1] += pk*F(2*factorial(k)*factorial(j), factorial(k+j+1))
        diff = fadd(fmul(row, poly), [-v for v in conv])
        factor = aj*R**(j+1)
        for k, v in enumerate(diff): out[k] += factor*v
    while out and out[-1] == 0: out.pop()
    return out


def f_inner(a, b):
    product = fmul(a, b)
    return R*sum(F(2, n+1)*v for n, v in enumerate(product) if not n % 2)


def iadd(a, b):
    out = [zero]*max(len(a), len(b))
    for j, v in enumerate(a): out[j] = out[j]+v
    for j, v in enumerate(b): out[j] = out[j]+v
    return out


def imul(a, b):
    out = [zero]*(len(a)+len(b)-1)
    for j, v in enumerate(a):
        if v.lo == v.hi == 0: continue
        for k, w in enumerate(b):
            if w.lo != 0 or w.hi != 0: out[j+k] = out[j+k]+v*w
    return out


def compose(poly, mid, half):
    out = [zero]
    for v in reversed(poly):
        new = [zero]*(len(out)+1)
        for j, w in enumerate(out):
            new[j] = new[j]+mid*w
            new[j+1] = new[j+1]+half*w
        new[0] = new[0]+I.coerce(v)
        out = new
    out.pop()  # Horner starts with a redundant degree-zero zero.
    return out


def double_factorial(n):
    out = 1
    while n > 0:
        out *= n; n -= 2
    return out


def legendre_coefficients(poly):
    degree = len(poly)-1
    out = []
    for k in range(degree+1):
        v = zero
        for n in range(k, degree+1, 2):
            factor = F((2*k+1)*factorial(n),
                      double_factorial(n-k)*double_factorial(n+k+1))
            v += poly[n]*I.rat(factor)
        out.append(v)
    return out


def dot(a, b, weights):
    out = zero
    for j in range(min(len(a), len(b), len(weights))):
        out += a[j]*b[j]*weights[j]
    return out


def lower30(value):
    return str(F(value.lo//10**120, 10**30))


def main(output=None):
    a0, h2, h3, w2, w3 = constants()
    assert h2.lo > I.rat(R).hi
    assert h3.hi < I.rat(2*R).lo
    assert atanh_log(I.rat(3, 5)).lo > I.rat(2*R).hi # log4>2R
    assert (w2+w3).hi < I.rat(6, 5).lo
    ps = legendre_polynomials(N)
    kernel = kernel_polynomial()
    cpoly = [F(0)]*25; spoly = [F(0)]*26
    for j in range(13):
        cpoly[2*j] = (R/2)**(2*j)/factorial(2*j)
        spoly[2*j+1] = (R/2)**(2*j+1)/factorial(2*j+1)
    leak = F(PENALTY)*2*R*(
        (2*(R/2)**26/factorial(26))**2+
        (2*(R/2)**27/factorial(27))**2)
    assert leak < F(1, 10**8)

    vp = [zero]*(2*P+1)
    for j in range(J):
        lam = F(4*j+1, 2)
        decay = exp_negative(I.rat(lam*R))
        for p in range(P+1):
            vp[2*p] += decay*I.rat(2*(lam*R)**(2*p)/(lam*factorial(2*p)))
    v0 = vp[0]
    vp[0] = zero
    assert sum((v.hi for v in vp)) < 6*SCALE
    hn = [F(0)]
    for j in range(1, N+1): hn.append(hn[-1]+F(1, j))
    tailgap = I.rat(hn[N])+a0+v0-w2-w3
    assert tailgap.lo > SCALE

    images = []
    for p in ps:
        g = [I.rat(v) for v in kernel_action(p, kernel)]
        g = iadd(g, imul(vp, [I.rat(v) for v in p]))
        penalty = fadd([PENALTY*f_inner(cpoly, p)*v for v in cpoly],
                       [PENALTY*f_inner(spoly, p)*v for v in spoly])
        g = iadd(g, [I.rat(v) for v in penalty])
        images.append(g)

    # Exact breakpoints for the two prime shifts, in t=x/R coordinates.
    shifts = [h2/I.rat(R), h3/I.rat(R)]
    ends = [-one, one-shifts[1], one-shifts[0],
            -one+shifts[0], -one+shifts[1], one]
    for a, b in zip(ends, ends[1:]): assert a.hi < b.lo
    Bpieces, Wpieces, weights = [], [], []
    for piece, (a, b) in enumerate(zip(ends, ends[1:])):
        mid, half = (a+b)/2, (b-a)/2
        bb = [legendre_coefficients(compose([I.rat(v) for v in p],mid,half)) for p in ps]
        wi = []
        for j, p in enumerate(ps):
            img = images[j][:]
            for shift, weight, leftlast, rightfirst in (
                    (shifts[0],w2,1,3), (shifts[1],w3,0,4)):
                if piece <= leftlast:
                    shifted = compose([I.rat(v) for v in p], shift, one)
                    img = iadd(img, [-weight*v for v in shifted])
                if piece >= rightfirst:
                    shifted = compose([I.rat(v) for v in p], -shift, one)
                    img = iadd(img, [-weight*v for v in shifted])
            wi.append(legendre_coefficients(compose(img,mid,half)))
        degree = max(len(w) for w in wi)
        wt = [I.rat(R)*(b-a)/I.rat(2*k+1) for k in range(degree)]
        Bpieces.append(bb); Wpieces.append(wi); weights.append(wt)

    matrix = [[zero for unused in range(N)] for unused in range(N)]
    for i in range(N):
        for j in range(N):
            val = zero
            for bp,wp,wt in zip(Bpieces,Wpieces,weights):
                val += dot(bp[i],wp[j],wt)
            matrix[i][j] = val
    for i in range(N):
        for j in range(i):
            val = (matrix[i][j]+matrix[j][i])/2
            matrix[i][j] = matrix[j][i] = val

    norms = [F(2*R, 2*n+1) for n in range(N)]
    gram = [[zero for unused in range(N)] for unused in range(N)]
    for bp,wp,wt in zip(Bpieces,Wpieces,weights):
        residuals = []
        for j in range(N):
            rr = wp[j][:]
            rr += [zero] * (max(len(p) for p in bp)-len(rr))
            for n in range(N):
                factor = matrix[n][j]/I.rat(norms[n])
                for k, v in enumerate(bp[n]): rr[k] -= factor*v
            residuals.append(rr)
        for i in range(N):
            for j in range(i+1): gram[i][j] += dot(residuals[i],residuals[j],wt)
    schur = [[zero for unused in range(N)] for unused in range(N)]
    for i in range(N):
        for j in range(i+1):
            val = matrix[i][j]-gram[i][j]  # conservative high inverse <=1
            if i == j:
                val += (I.rat(hn[i])+a0+v0-I.rat(SCHUR_FLOOR))*I.rat(norms[i])
            schur[i][j] = schur[j][i] = val

    pivots = []; low = [[zero for unused in range(N)] for unused in range(N)]
    for i in range(N):
        pivot = schur[i][i]
        for k in range(i): pivot -= low[i][k].sq()*pivots[k]
        assert pivot.lo > 0, ("nonpositive certified pivot", i, lower30(pivot))
        pivots.append(pivot)
        low[i][i] = one
        for j in range(i+1,N):
            value = schur[j][i]
            for k in range(i): value -= low[j][k]*low[i][k]*pivots[k]
            low[j][i] = value/pivot

    coercivity = SCHUR_FLOOR/F(129)-leak
    assert coercivity > F(1, 10**6)
    result = {
        "status": "PASS", "arithmetic": "outward integer intervals, scale10^150; exact polynomial integrals",
        "support_diameter": "11/10", "low_legendre_modes":N,
        "gamma_heads":J,"outside_Taylor_degree":2*P,
        "smooth_kernel_degree":DEG,"smooth_kernel_error":str(KERNEL_ERROR),
        "pole_polynomial_degrees":[24,25],"pole_penalty":PENALTY,
        "tail_gap_lower":lower30(tailgap),"Schur_floor":str(SCHUR_FLOOR),
        "LDL_pivot_lower_bounds":[lower30(v) for v in pivots],
        "pool_leakage_upper":str(leak),"coercivity_lower":str(coercivity),
        "claimed_coercivity":"1/1000000", "prime_powers":"all retained; only n=2,3 have nonzero support intersections"
    }
    if output: Path(output).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--output")
    main(parser.parse_args().output)
