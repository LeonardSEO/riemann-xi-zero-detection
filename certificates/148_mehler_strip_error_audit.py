#!/usr/bin/env python3
"""Outward-rounded certificate for the explicit constant in 144.32.
The analytic infinite-series/tail proof is in 144_mehler_modular_approximants.md.
This certificate only verifies its elementary final constant is <1500.
"""
import json
from pathlib import Path
from flint import arb,ctx
ctx.prec=256
pi=arb.pi();q=(-3*pi).exp()
B=[1/(1-q),pi*(1+q)/(1-q)**3,
   pi**2*(1+11*q+11*q**2+q**3)/(1-q)**5,
   pi**3*(1+57*q+302*q**2+302*q**3+57*q**4+q**5)/(1-q)**7]
CI=24*(2/pi)**5*(2*B[3]+16*B[2]+18*B[1]+arb(9)*B[0]/4)
C0=16*(2*B[2]+3*B[1])/pi**3
CT=B[2]+3*B[1]+arb(3)*B[0]/2+C0
S=(7/pi)**(arb(7)/2)*(-arb(7)/2).exp()
C=2*CI+2*CT*S+arb(3)/2
assert C<1500
assert C>1439
out={'precision_bits':ctx.prec,'B':[x.str(65) for x in B],
     'C_I':CI.str(65),'C_0_star':C0.str(65),'C_T':CT.str(65),
     'S_star':S.str(65),'C_star':C.str(65),'bound':'1439<C_star<1500',
     'scope':'elementary constant only; all-height weighted-L1 proof in144.11'}
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'C_star':C.str(30),'C_star_lt_1500':True,'certificate':str(p)}))
