#!/usr/bin/env python3
"""Actual-theta cold A-birth and small-time IFT: complete Arb proof.

No truncated heat series is used.  The heated MP coefficients are
computed by their exact exponential generating function.  The proof
of each global bound and the real contraction is in route report218F.
Independent reconstruction remains a separate control status.
"""
import hashlib
import json
from pathlib import Path
from math import factorial
from flint import arb, arb_series, ctx

ctx.prec = 8192
M, N = 1400, 50
ctx.cap = M + 1
a = arb(1) / 2
A = arb("0.005095345831716476276508107776687142746042")
kappa = arb(1) / 16
x_num = 22765551187540060734255795318900976
t_num = 0
x_den, t_den = 10**33, 1
x0, t0 = arb(x_num) / x_den, arb(t_num) / t_den
radius = arb(1) / 10**18
rho = b = arb(97) / 100
rz, rk = arb(1) / 10, arb(1) / 1000

q = [arb(0) for _ in range(M + 1)]
for n in range(1, N + 1):
    pn = arb.pi() * n * n
    ee = (-pn).exp()
    l0, l1 = arb(1), a - 2 * pn
    q[0] += ee
    q[1] += ee * l1
    for m in range(1, M):
        l2 = ((2 * m + a - 2 * pn) * l1 - (m - a) * l0) / (m + 1)
        q[m + 1] += ee * l2
        l0, l1 = l1, l2
pc = arb(1)
eps_theta = arb(2).sqrt() * (-arb.pi() * (N + 1)**2 / 3).exp()
eps_theta /= 1 - (-arb.pi() * (2 * N + 3) / 3).exp()
for n in range(M + 1):
    q[n] -= a * ((-1) ** n) * pc
    if n % 2 == 0:
        q[n] += arb(0, eps_theta * (arb(2) ** n))
    pc *= (n + a) / (n + 1)

L = (1 + x0*x0).log()
beta = A * (kappa * L * L).exp()
k = (beta - 1) / (beta + 1)
hp = 4*kappa*x0*L/(1+x0*x0)
hpp = 4*kappa*((1-x0*x0)*L+2*x0*x0)/(1+x0*x0)**2
g = (1-k*k)/2
kp, kpp = g*hp, g*(hpp-k*hp*hp)

Mq = (1+rho).sqrt()/(2*(1-rho)) + 1/(2*(1-rho).sqrt())
ab, hb = b.atan(), b.atanh()
Lb2 = ab*ab + hb*hb

def cbox(X, Y, T):
    return (1-b*b)**(-arb(1)/4) * (X*ab + Y*hb + T*Lb2).exp()

def total_bound(K, X, Y, T):
    v = K/(rho*b)
    assert v < 1
    return Mq*cbox(X,Y,T)*(1+v)/(1-v)**2

v = (abs(k)+rk)/(rho*b)
assert v < 1
tail = Mq*cbox(x0+rz,rz,t0)*v**(M+1)
tail *= (2*M+3)/(1-v)+2*v/(1-v)**2

ss = arb_series([0,1])
ats = ss.atan()
base = (-(1+ss*ss).log()/4+x0*ats-t0*ats*ats).exp()
jets = {}
pp = arb_series([1])
for ell in range(4):
    coeff = base*pp
    fac = arb(1)
    sums = [arb(0),arb(0),arb(0)]
    for n in range(M+1):
        if n % 2 == 0:
            value = q[n]*k**n*((-1)**(n//2))*fac*coeff[n]
            sums[0] += value
            sums[1] += n*value/k
            sums[2] += n*(n-1)*value/(k*k)
        fac *= (n+1)/(a+n)
    for kd in range(3):
        err = factorial(ell)*factorial(kd)*tail/rz**ell/rk**kd
        jets[ell,kd] = sums[kd]+arb(0,err)
    pp *= ats

f = jets[0,0]
fx = jets[1,0]+kp*jets[0,1]
fxx = jets[2,0]+2*kp*jets[1,1]+kpp*jets[0,1]+kp*kp*jets[0,2]
physical_ft = -jets[2,0]
chi = g*jets[0,1]
ft = chi
fxt = g*(jets[1,1]+kp*jets[0,2]-k*hp*jets[0,1])

# On |z-x0|,|ell|<=1/10, A=A0 exp(ell), |k-k0|<1/50.
X = x0+rz
den = 1+x0*x0-2*x0*rz-rz*rz
arg_bound = 2*X*rz/den
log_bound = (1+X*X).log()+arg_bound
hp_bound = 4*kappa*X*log_bound/den
assert hp_bound < arb(1)/10
beta_diff = beta*((arb(11)/100).exp()-1)
k_diff = 2*beta_diff/((beta+1-beta_diff)*(beta+1))
assert k_diff < arb(1)/50
K = abs(k)+arb(1)/50
Cf = total_bound(K,x0+rz,rz,arb(0))
Rhalf = rz/2
B2 = 2*Cf/Rhalf**2
B3 = 6*Cf/Rhalf**3

det = fx*fxt-ft*fxx
assert not det.contains(0)
B = [[(fxt/det).mid(),(-ft/det).mid()],
     [(-fxx/det).mid(),(fx/det).mid()]]
J0 = [[fx,ft],[fxx,fxt]]
Bn = max((abs(B[0][0])+abs(B[0][1])).upper(),
         (abs(B[1][0])+abs(B[1][1])).upper())
center_dev = []
for row in range(2):
    dd = arb(0)
    for col in range(2):
        entry = (1 if row==col else 0) - sum((B[row][j]*J0[j][col] for j in range(2)),arb(0))
        dd += abs(entry)
    center_dev.append(dd.upper())
qcontract = max(center_dev)+4*radius*Bn*B3
resid = max(abs(B[0][0]*f+B[0][1]*fx).upper(),
            abs(B[1][0]*f+B[1][1]*fx).upper())
assert qcontract < arb(1)/2
assert resid+qcontract*radius < radius

# Node signs throughout the real contraction square.
den_node = fxx+arb(0,2*radius*B3)
Cphysical = 2*total_bound(K,x0+rz+arb(1)/100,rz+arb(1)/100,arb(0))/(arb(1)/100)**2
ft_node = physical_ft+arb(0,2*radius*Cphysical/Rhalf)
Kwide = K+arb(1)/100
Cchi = total_bound(Kwide,x0+rz,rz,arb(0))/(arb(1)/100)
chi_node = chi+arb(0,2*radius*Cchi/Rhalf)
assert not den_node.contains(0)
ratio_A = chi_node/den_node
ratio_t = ft_node/den_node
ratio_birth = ratio_A
assert ratio_A > 0 and ratio_t < 0

def out(v):
    return str(v.str(40))

report = {
    "status":"CERTIFIED_ACTUAL_THETA_DOUBLE_NODE_INTERNAL_DEVELOPMENT",
    "independent_audit":"PENDING",
    "RH_status":"UNRESOLVED",
    "A_center_exact_decimal":"0.005095345831716476276508107776687142746042", "center_log_relative_A":"0", "kappa":"1/16", "precision_bits":ctx.prec,
    "MP_degree":M,"theta_terms":N,"heat_orders":"ALL_EXACT_POLYNOMIAL_HEAT",
    "center_x":{"numerator":str(x_num),"denominator":str(x_den)},
    "center_t":{"numerator":str(t_num),"denominator":str(t_den)},
    "real_contraction_square_radius":"1/10^18 in (x,log(A/A0))",
    "beta_at_center":out(beta),"k_at_center":out(k),
    "theta_coefficient_error":out(eps_theta),"MP_tail_box":out(tail),
    "analytic_composed_bound":out(Cf),"third_jet_bound":out(B3),
    "inverse_jacobian_norm":out(Bn),"contraction_bound":out(qcontract),
    "center_newton_residual":out(resid),
    "f":out(f),"fx":out(fx),"fxx_node":out(den_node),
    "partial_logA_ratio":out(ratio_A),"partial_t_ratio":out(ratio_t),
    "path_t_derivative":"0; actual cold logA birth, IFT extends to every small t>0",
    "positive_birth_ratio":out(ratio_birth),
    "scope":"One actual composed Mehler node at finite A; no cofinal impurity or zeta counterexample.",
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
output = Path(__file__).with_suffix('.json')
output.write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
