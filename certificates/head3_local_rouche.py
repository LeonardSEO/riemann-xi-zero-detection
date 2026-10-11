#!/usr/bin/env python3
"""Rigorous local Rouché certificate for a nonreal zero of theta head N=3.

This is a finite-theta-head certificate, not a zeta-zero certificate.
Requires cached python-flint 0.9.0. All proof quantities use outward
rounded Arb/Acb arithmetic. The centre and radius are exact dyadic numbers.
"""
import json, argparse, math
from flint import acb, arb, ctx

ctx.prec=640
pi=arb.pi(); I=acb(0,1)

def F(n,z):
    out=acb(0); z=acb(z)
    for k in range(1,n+1):
        t=pi*k*k
        for eps in (-1,1):
            alpha=eps*I*z/2
            out += t**(-acb(1)/4-alpha)*(2*acb(t).gamma_upper(acb(9)/4+alpha)-3*acb(t).gamma_upper(acb(5)/4+alpha))
    return out

def moment_bound(n,j):
    # e^(2u)>=1+2u implies the half-line positive integrand bound;
    # the factor 2 includes both halves of the even kernel.
    return 8*pi*pi*math.factorial(j)*sum((arb(k**4)*(-pi*k*k).exp()/(2*pi*k*k-5)**(j+1) for k in range(1,n+1)),arb(0))

def midpoint_newton(n,z,h):
    for _ in range(8):
        value=F(n,z)
        dh=(F(n,z+h)-F(n,z-h))/(2*h)
        if dh.contains(0): raise RuntimeError('Derivative difference overlaps zero')
        z=(z-value/dh).mid()
    return z

def exact_floor_midpoint(x,power):
    scaled=x.mid()*arb(2)**power
    num=scaled.floor().unique_fmpz()
    if num is None:raise RuntimeError('Cannot extract exact dyadic centre')
    return int(num)

p=argparse.ArgumentParser();p.add_argument('--certificate');args=p.parse_args()
n=3; power=120; h=arb(2)**(-120); radius=arb(2)**(-85)
z0=acb(arb('67.88018965514762'),arb('0.477343841770823'))
znew=midpoint_newton(n,z0,h)
xnum=exact_floor_midpoint(znew.real,power);ynum=exact_floor_midpoint(znew.imag,power)
centre=acb(arb(xnum)*arb(2)**(-power),arb(ynum)*arb(2)**(-power))
assert centre.real.rad()==0 and centre.imag.rad()==0
m2=moment_bound(n,2);m3=moment_bound(n,3)
if not (m2>0 and m3>0):raise RuntimeError('Moment bounds not positive')
value=F(n,centre)
difference=(F(n,centre+h)-F(n,centre-h))/(2*h)
derivative_error=m3*h*h/6
derivative_lower=difference.abs_lower()-derivative_error
rouche_error=value.abs_upper()+m2*radius*radius/2
rouche_main=derivative_lower*radius
if not derivative_lower>0:raise RuntimeError('No nonzero derivative bound')
if not centre.imag-radius>0:raise RuntimeError('Disk hits real axis')
if not centre.imag+radius<arb(1)/2:raise RuntimeError('Disk leaves closed strip')
if not rouche_error<rouche_main:raise RuntimeError('Rouche inequality fails')
result={
 'status':'CERTIFIED_FINITE_HEAD_NONREAL_ZERO', 'not_a_zeta_zero':True,
 'N':n, 'precision_bits':ctx.prec, 'python_flint_version':'0.9.0',
 'fourier_convention':'integral Phi_N(u) exp(-i*z*u) du; Phi_N even',
 'centre_real_numerator':str(xnum), 'centre_imag_numerator':str(ynum),
 'centre_denominator_power_of_two':power, 'radius_power_of_two':-85,
 'derivative_step_power_of_two':-120, 'centre_display':str(centre),
 'F_centre':str(value), 'M2_upper':str(m2), 'M3_upper':str(m3),
 'difference_quotient':str(difference), 'derivative_error_upper':str(derivative_error),
 'derivative_modulus_lower':str(derivative_lower),
 'rouche_error_upper':str(rouche_error), 'rouche_main_lower':str(rouche_main),
 'rouche_ratio_upper':str(rouche_error/rouche_main),
 'number_of_zeros_with_multiplicity_in_disk':1,
 'proof':'Rouche vs F_N_prime(c)*(z-c). Value plus M2*r^2/2 is less than derivative lower bound times r. Derivative lower bound follows from symmetric finite difference with error M3*h^2/6.'}
print(json.dumps(result,indent=2))
if args.certificate:
    with open(args.certificate,'w') as out:json.dump(result,out,indent=2);out.write('\n')
