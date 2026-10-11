#!/usr/bin/env python3
"""Candidate Acb Rouche certificate for the ACTUAL modular Mehler family.
Uses the entire E_k factor of report144, all theta coefficients with
outward-rounded Cauchy tails, and the full infinite MP-series tail.
No zeta zero is claimed. One local zero of M_alpha only.
"""
import json
from pathlib import Path
from flint import arb,acb,ctx
ctx.prec=4096
M=1500;N=50;a=arb(1)/2;k=arb(177)/200
q=[arb(0) for _ in range(M+1)]
for n in range(1,N+1):
    tau=arb.pi()*n*n;e=(-tau).exp();l0=arb(1);q[0]+=e
    l1=a-2*tau;q[1]+=e*l1
    for m in range(1,M):
        l2=((2*m+a-2*tau)*l1-(m-a)*l0)/(m+1)
        q[m+1]+=e*l2;l0,l1=l1,l2
c=arb(1)
for m in range(M+1):
    q[m]-=a*((-1)**m)*c;c*=(m+a)/(m+1)
rem=arb(2).sqrt()*(-arb.pi()*(N+1)**2/3).exp()/(1-(-arb.pi()*(2*N+3)/3).exp())
for m in range(M+1):q[m]+=arb(0,rem*(2**m))
# Only even coefficients are used: exact Jacobi parity, not numerically rounded odd values.
coeff=[q[n]*k**n*((-1)**(n//2)) if n%2==0 else arb(0) for n in range(M+1)]
def E_and_derivative(z):
    r0=acb(1);r1=z/a;dr0=acb(0);dr1=1/a
    val=coeff[0]*r0;der=acb(0)
    for n in range(1,M):
        r2=(z*r1-n*r0)/(a+n);dr2=(r1+z*dr1-n*dr0)/(a+n)
        val+=coeff[n+1]*r2;der+=coeff[n+1]*dr2
        r0,r1=r1,r2;dr0,dr1=dr1,dr2
    return val,der
# The search was performed separately. This fixed certificate uses only integers.
nums=[30561398818003719921715546903735782215,597079250754110273230973026105322206]
c=acb(arb(nums[0])/(2**120),arb(nums[1])/(2**120));rad=arb(1)/(2**60)
# Uniform Cauchy bounds: exact rational rho=b=97/100.
rho=b=arb(97)/100;t=k/(rho*b)
MQ=(1+rho).sqrt()/(2*(1-rho))+1/(2*(1-rho).sqrt())
def Cb(X,Y):return (1-b*b)**(-arb(1)/4)*(X*b.atan()+Y*b.atanh()).exp()
def tail(X,Y):
    return MQ*Cb(X,Y)*t**(M+1)*((2*M+3)/(1-t)+2*t/(1-t)**2)
X=abs(c.real);Y=abs(c.imag)
T=tail(X,Y);Td=tail(X+1,Y+1)
Cfull=MQ*Cb(X+2,Y+2)*(1+t)/(1-t)**2
val,der=E_and_derivative(c)
v_upper=abs(val)+T;d_lower=abs(der)-Td
lhs=abs(val)+T+Td*rad+Cfull*rad*rad
rhs=abs(der)*rad
assert d_lower>0
assert lhs<rhs
assert Y+rad<arb(1)/2 and Y-rad>0
out={'scope':'one actual Mehler zero; no zeta/RH claim','precision_bits':ctx.prec,
     'k_numerator':177,'k_denominator':200,'alpha_numerator':377,'alpha_denominator':23,
     'center_real_numerator':str(nums[0]),'center_imag_numerator':str(nums[1]),
     'center_denominator_exponent':120,'radius_denominator_exponent':60,
     'MP_degree':M,'theta_cutoff':N,'coefficient_circle_radius':'1/2',
     'Q_circle_radius':'97/100','MP_generating_circle_radius':'97/100',
     'theta_coefficient_tail_base':rem.str(70),'Q_circle_bound':MQ.str(70),
     'series_ratio':t.str(70),'entire_value_tail':T.str(70),'entire_derivative_tail':Td.str(70),
     'entire_cauchy_bound_radius_2':Cfull.str(70),
     'finite_value':val.str(70),'finite_derivative':der.str(70),
     'actual_value_upper':v_upper.str(70),'actual_derivative_lower':d_lower.str(70),
     'rouche_error_upper':lhs.str(70),'rouche_linear_lower':rhs.str(70),
     'rouche_ratio':(lhs/rhs).str(70),'center':c.str(70),'simple_zero_count':1}
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'center':c.str(45),'rouche_ratio':(lhs/rhs).str(25),'certificate':str(p)}))
