#!/usr/bin/env python3
"""Internal rigorous ball certificate for the actual 348 coherent erf partner.

Run with the existing runtime, without installing packages:
    PYTHONPATH=/tmp/riemann-flint python3 synthesis/certificates/348_erf_phase_sign.py

All theta indices are evaluated through Acb's complete Jacobi theta derivatives.
The two finite quadratures have explicit rigorous omitted-tail bounds below.
This certifies a phase sign of an auxiliary partner, not any zero of zeta.
"""
from pathlib import Path
import hashlib
import json
import flint
from flint import arb, acb, acb_series, ctx

ctx.prec = 192
ctx.cap = 2
PI = arb.pi()
I = acb(0, 1)
B = PI / 4
EPS = arb(1) / 256
H = B - EPS
U = arb(3)
X = arb(3)
KAPPA = arb(1)

def phi(u, analytic):
    # theta_3(z,tau) uses q=exp(pi*i*tau), with pi included in z.
    # At tau=i*exp(2u), theta_3''(0)=-8*pi^2*sum n^2 exp(-pi*n^2*exp(2u)),
    # theta_3''''(0)=32*pi^4*sum n^4 exp(-pi*n^2*exp(2u)).
    tau = I * (2 * u).exp()
    # The callback must verify holomorphy on the full complex ball.
    if not tau.imag > 0:
        return acb("nan")
    r2 = acb(0).modular_theta(tau, r=2)[2]
    r4 = acb(0).modular_theta(tau, r=4)[2]
    return ((9 * u / 2).exp() * r4 / (8 * PI * PI)
            + 3 * (5 * u / 2).exp() * r2 / (4 * PI))

ERF = acb(KAPPA * X).erf()
ERFP = 2 * KAPPA / PI.sqrt() * (-(KAPPA * X)**2).exp()
OPTIONS = dict(abs_tol=arb("1e-40"), rel_tol=arb("1e-40"),
               eval_limit=40000, depth_limit=35)

# Entire integrands: no taking a real part inside the analytic callback.
sine = acb.integral(lambda u, a: phi(u, a) * (X * u).sin(),
                    0, U, **OPTIONS)
sine_prime = acb.integral(
    lambda u, a: u * phi(u, a) * (X * u).cos(), 0, U, **OPTIONS)
cusp = acb.integral(
    lambda v, a: phi(I * v, a) *
        ((X * v).sinh() - ERF * (X * v).cosh()), 0, H, **OPTIONS)
cusp_prime = acb.integral(
    lambda v, a: phi(I * v, a) *
        (v * (X * v).cosh() - ERF * v * (X * v).sinh()
         - ERFP * (X * v).cosh()), 0, H, **OPTIONS)

# Original Xi and its exact derivative via the complete meromorphic completion.
z = acb_series([acb(X), acb(1)])
s = acb("1/2") + I * z
xi = s * (s - 1) / 2 * (-s * PI.log() / 2).exp() * (s / 2).gamma() * s.zeta()
P, PP = xi[0].real, xi[1].real

# Physical tail: for u>=U>=0,
# sum n^4 exp(-q*n^2) <= 300 exp(-q), q=pi*exp(2u)>=log 2,
# since sum n^4 (1/2)^(n-1)=300 and n^2-1>=n-1.
# Phi <= 1200*pi^2 exp(9u/2-pi*exp(2u)).
# Concavity makes the tangent-slope integral bound exact.
# u <= exp(u) bounds the first-derivative tail.
def physical_tail(j):
    power = arb(9) / 2 + j
    q = PI * (2 * U).exp()
    assert 2 * q > power
    return 1200 * PI * PI * (power * U - q).exp() / (2 * q - power)

PHYS0, PHYS1 = physical_tail(0), physical_tail(1)

# Exact Q(H),Q'(H), using the same completed modular normalization:
# Q'=i exp(iv/2) [theta_3(0,i exp(2iv))/2
#                 + exp(2iv)*theta_3''(0,i exp(2iv))/(2pi)].
# The derivative here uses all n in Z: d_w theta=theta_3''/(4pi).
w = (2 * I * H).exp()
theta_tau = I * w
theta0 = acb(0).modular_theta(theta_tau)[2]
theta2 = acb(0).modular_theta(theta_tau, r=2)[2]
QH = (I * H / 2).exp() * theta0
QPH = I * (I * H / 2).exp() * (theta0 / 2 + w * theta2 / (2 * PI))
Q_BOUND = abs(QH).upper()
QP_BOUND = abs(QPH).upper()
assert QH.real > 0

# Actual theta inequality (348.10): integral_H^b Q <= 3 EPS Q(H).
# Two exact integrations by parts eliminate Q'' from the omitted cusp.
# R=sinh(Xv)-erf(X)cosh(Xv), |R|<=2exp(Xb), |R'|<=2X exp(Xb).
# For d_X R: |d_X R|<=4exp(Xb), |d_v d_X R|<=(2+4X)exp(Xb).
# These bounds use b<1, |erf(X)|<1, |erf'(X)|<2; X=3,kappa=1.
assert B < 1
assert abs(ERF) < 1
assert abs(acb(ERFP)) < 2
CUSP0 = (X * B).exp() * (
    QP_BOUND + X * Q_BOUND + 3 * EPS * (X * X + arb(1)/4) * Q_BOUND)
CUSP1 = (X * B).exp() * (
    2 * QP_BOUND + (1 + 2 * X) * Q_BOUND
    + 3 * EPS * (2 * (X * X + arb(1)/4) + 2 * X) * Q_BOUND)

ERR0 = PHYS0 + CUSP0
ERR1 = PHYS1 + CUSP1
K = (sine + cusp).real + arb(0, ERR0.upper())
KP = (sine_prime + cusp_prime).real + arb(0, ERR1.upper())
WRONSKIAN = P * KP - PP * K
assert WRONSKIAN > arb("0.0124")
assert WRONSKIAN < arb("0.0125")

data = {
    "status": "CERTIFIED_AUXILIARY_PHASE_SIGN",
    "rh_status": "UNRESOLVED",
    "software": {"python_flint": flint.__version__, "precision_bits": ctx.prec},
    "scope": "actual coherent erf partner at tau=0, kappa=1, x=3",
    "normalization": "Xi(z)=xi(1/2+iz); E=Xi/2+iK",
    "parameters": {"x": "3", "kappa": "1", "physical_cut": "3",
                   "cusp_cut": "pi/4-1/256"},
    "ball_values": {n: str(v) for n,v in [
        ("Xi", P), ("Xi_prime", PP), ("K", K), ("K_prime", KP),
        ("Xi_K_prime_minus_Xi_prime_K", WRONSKIAN), ("Q_at_cut", QH),
        ("Q_prime_at_cut", QPH), ("physical_tail_K", PHYS0),
        ("physical_tail_K_prime", PHYS1), ("cusp_tail_K", CUSP0),
        ("cusp_tail_K_prime", CUSP1)]},
    "proven_rational_interval": ["0.0124", "0.0125"],
    "tail_dependencies": ["348A/348.10 and exact two-fold integration by parts",
                          "complete original theta series, all n>=1"],
    "not_claimed": "No off-line zeta zero, no global HB property.",
}
source = Path(__file__).read_bytes()
data["script_sha256"] = hashlib.sha256(source).hexdigest()
destination = Path(__file__).with_suffix(".json")
destination.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(data, ensure_ascii=False, indent=2))
