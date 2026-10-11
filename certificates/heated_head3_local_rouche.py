#!/usr/bin/env python3
"""Rigorous local certificate: heated head N=3, t=1/(16*pi).
A finite head with heat, not a zeta-zero certificate. Uses exact dyadic
centre/radius and Acb analytic integration plus an elementary tail bound.
"""
import math,json,argparse
from flint import arb,acb,ctx
ctx.prec=600
pi=arb.pi();t=1/(16*pi);U=3;I=acb(0,1)

def phi(k,u):
    T=pi*k*k*(2*u).exp()
    return (u/2-T).exp()*(4*T*T-6*T)

def tail(korder=0):
    # For korder=0,1,2,3, u^k <= exp(u^2), 0<=t<=1.
    # 2u^2+5u <= (pi/2) exp(2u) for u>=0.
    # Thus integrand absolute tail is <= 4pi^2*n^4 exp(-pi*n^2 exp(2u)/2).
    # e^(2u)>=e^(2U)(1+2(u-U)) integrates this bound explicitly.
    B=(2*arb(U)).exp()
    return 8*pi/B*sum((k*k*(-pi*k*k*B/2).exp() for k in (1,2,3)),arb(0))

def H(z,j=0):
    phase=(lambda w:w.cos(),lambda w:-w.sin(),lambda w:-w.cos(),lambda w:w.sin())[j%4]
    value=acb(0)
    for k in (1,2,3):
        value+=2*acb.integral(lambda u,analytic: (t*u*u).exp()*phi(k,u)*(u**j)*phase(z*u),0,U,
            rel_tol=arb(2)**(-470),abs_tol=arb(2)**(-470),depth_limit=40,eval_limit=200000)
    b=tail(j).abs_upper()
    return value+acb(arb(0,b),arb(0,b))

# Elementary all-u moment bound valid for orders 0..3 by the same inequality.
M=8*pi*sum((k*k*(-pi*k*k/2).exp() for k in (1,2,3)),arb(0))
assert M<6
c=acb(arb('67.862'),arb('.430'))
for _ in range(7):
    der=H(c,1)
    if der.contains(0):raise RuntimeError('Newton derivative overlaps zero')
    c=(c-H(c)/der).mid()
power=120
xnum=int((c.real*arb(2)**power).floor().unique_fmpz());ynum=int((c.imag*arb(2)**power).floor().unique_fmpz())
c=acb(arb(xnum)*arb(2)**(-power),arb(ynum)*arb(2)**(-power))
r=arb(2)**(-90);value=H(c);derivative=H(c,1)
err=value.abs_upper()+M*r*r/2;main=derivative.abs_lower()*r
assert c.real.rad()==0 and c.imag.rad()==0
assert c.imag-r>0 and c.imag+r<arb(1)/2
assert err<main
result={'status':'CERTIFIED_HEATED_FINITE_HEAD_NONREAL_ZERO','not_a_zeta_zero':True,'N':3,'heat_parameter':'1/(16*pi)','precision_bits':600,'centre_real_numerator':str(xnum),'centre_imag_numerator':str(ynum),'centre_denominator_power_two':120,'radius_power_two':-90,'centre':str(c),'value':str(value),'derivative':str(derivative),'moment_0_to_3_bound':str(M),'tail_bound_orders_0_to_3':str(tail()),'rouche_error':str(err),'rouche_main':str(main),'rouche_ratio':str(err/main),'zeros_with_multiplicity_in_disk':1}
p=argparse.ArgumentParser();p.add_argument('--certificate');args=p.parse_args();print(json.dumps(result,indent=2))
if args.certificate:
    with open(args.certificate,'w') as f:json.dump(result,f,indent=2);f.write('\n')
