#!/usr/bin/env python3
"""Arb/Acb contour certificate attempt for the N=1 theta head.

Run with python-flint 0.9.0 (the bundled runtime may need PYTHONPATH pointing
at its installation).  The script uses the exact incomplete-gamma formula at
point endpoints and an Acb integral enclosure on each parameter rectangle.
It certifies the full rectangle [-94,94] x [-1/2,1/2] when run successfully.
"""

from flint import *
import flint
import argparse, hashlib, json, math, time
from fractions import Fraction
from pathlib import Path
ctx.prec=140  # 140-bit outward-rounded Arb/Acb arithmetic
pi=arb.pi(); I=acb(0,1)

def exact_arb(value):
    if isinstance(value, Fraction):
        return arb(value.numerator) / value.denominator
    return arb(value)

def exact_point(x, y):
    return acb(exact_arb(x), exact_arb(y))

def F_point(z):
    z=acb(z); s=acb(0)
    for eps in (-1,1):
        e=acb(eps)
        a=acb(arb(9)/4 + e*I*z/2)
        b=acb(arb(5)/4 + e*I*z/2)
        s += pi**(-e*I*z/2)*(2*acb(pi).gamma_upper(a)-3*acb(pi).gamma_upper(b))
    return pi**(-acb(1)/4)*s

def phi(u):
    u=acb(u)
    return 2*pi*acb.exp(acb(5)*u/2)*(2*pi*acb.exp(2*u)-3)*acb.exp(-pi*acb.exp(2*u))

def tail_bound(U):
    # For |Im z| <= 1/2, t = pi exp(2u), and u >= U >= 0,
    # 2 int_U^infty |phi(u) cos(zu)| du is bounded by
    # 2 pi^(-1/2) [2 Gamma(5/2,t0) + 3 Gamma(3/2,t0)].
    # Since Gamma(s,t0) <= 5 t0^(s-1) exp(-t0) for t0 >= 1 and
    # s in {3/2,5/2}, the elementary expression below is rigorous and
    # avoids evaluating gamma_upper at the very large t0.

    t0=arb.pi()*arb(2*U).exp()
    q=10*t0**(arb(3)/2)+15*t0**(arb(1)/2)
    return 2*arb.pi()**(-arb(1)/2)*q*(-t0).exp()

def F_box(x0,x1,y0,y1,U=2):
    z=acb(arb((x0+x1)/2,(x1-x0)/2),arb((y0+y1)/2,(y1-y0)/2))
    q=2*acb.integral(lambda u,a: phi(u)*(z*u).cos(),0,arb(U),rel_tol=arb(2)**(-110),abs_tol=arb(2)**(-110),depth_limit=30)
    tb=tail_bound(U)
    rad=tb.abs_upper()
    return q+acb(arb(0,rad),arb(0,rad))


def derivative_tail_bound(k, U=2):
    """Elementary bound for the omitted kth-z derivative tail."""
    # u^k <= (2k)^k exp(u/2) (with C_0=1); after t=pi exp(2u),
    # the two gamma tails are bounded by 5 t0^(s-1) exp(-t0).
    coefficient = 1 if k == 0 else (2 * k) ** k
    t0 = arb.pi() * arb(2 * U).exp()
    q = 10 * t0 ** (arb(7) / 4) + 15 * t0 ** (arb(3) / 4)
    return 2 * coefficient * arb.pi() ** (-arb(3) / 4) * q * (-t0).exp()


def derivative_point(z, k, U=2):
    """Acb enclosure of F1^(k)(z), using the finite integral and tail bound."""
    if k == 0:
        return F_point(z)
    phase = (lambda w: w.cos(), lambda w: -w.sin(),
             lambda w: -w.cos(), lambda w: w.sin())[k % 4]
    q = 2 * acb.integral(
        lambda u, analytic: phi(u) * (u ** k) * phase(z * u),
        0, arb(U), rel_tol=arb(2) ** (-110),
        abs_tol=arb(2) ** (-110), depth_limit=30)
    radius = derivative_tail_bound(k, U).abs_upper()
    return q + acb(arb(0, radius), arb(0, radius))


def absolute_derivative_bound(k, U=2):
    """Bound M_k >= sup_|Im z|<=1/2 |F1^(k)(z)|."""
    q = 2 * acb.integral(
        lambda u, analytic: phi(u) * (u ** k) * acb(u / 2).exp(),
        0, arb(U), rel_tol=arb(2) ** (-110),
        abs_tol=arb(2) ** (-110), depth_limit=30)
    return q.real.upper() + derivative_tail_bound(k, U).abs_upper()


def F_segment(x0, x1, y0, y1, order=4):
    """Acb enclosure of F1 on one horizontal or vertical segment.

    A centered Taylor polynomial at the segment midpoint is enlarged by the
    rigorous global moment bound M_order * radius^order/order!.  This avoids
    the severe dependency blow-up of direct interval quadrature at |x|~94.
    """
    za, zb = exact_point(x0, y0), exact_point(x1, y1)
    zc = (za + zb) / 2
    radius = max((exact_arb(x1) - exact_arb(x0)).abs_upper(),
                 (exact_arb(y1) - exact_arb(y0)).abs_upper()) / 2
    dz = acb(arb(0, radius), 0) if y0 == y1 else acb(0, arb(0, radius))
    assert (zc + dz).contains(za) and (zc + dz).contains(zb)
    value = acb(0)
    for k in range(order):
        value += derivative_point(zc, k) * dz ** k / math.factorial(k)
    remainder = MOMENTS[order] * radius ** order / math.factorial(order)
    return value + acb(arb(0, remainder), arb(0, remainder))

def check_seg(a,b):
    x0,y0=a;x1,y1=b
    img=F_segment(min(x0,x1),max(x0,x1),min(y0,y1),max(y0,y1))
    if img.abs_lower() <= 0: return False,'zero-image',img,None
    m=acb(img.real.mid(),img.imag.mid())
    if m==0: return False,'zero-mid',img,None
    rr=m.conjugate()*img
    if rr.real.lower() <= 0: return False,'not-halfplane',img,rr.real
    fa, fb = F_point(exact_point(x0,y0)), F_point(exact_point(x1,y1))
    if fa.abs_lower() <= 0 or fb.abs_lower() <= 0:
        return False,'zero-endpoint',img,None
    q=fb/fa
    aa=q.arg()
    if aa.lower() <= -pi.lower() or aa.upper() >= pi.lower():
        return False,'arg-branch',img,aa
    return True,q, img, aa

MOMENTS = [absolute_derivative_bound(k) for k in range(5)]

def cert(R=94, h=Fraction(1,10), maxdepth=10):
    half=Fraction(1,2)
    base=[((-R,-half),(R,-half)),((R,-half),(R,half)),
          ((R,half),(-R,half)),((-R,half),(-R,-half))]
    allq=[]; failures=[]; n=0
    def rec(a,b,d,si):
      nonlocal n
      n+=1; ok,q,img,aux=check_seg(a,b)
      if ok:
        allq.append((si,a,b,q.arg(),img)); return
      if d>=maxdepth:
        failures.append((si,a,b,q,img,aux)); return
      xm=(Fraction(a[0])+b[0])/2; ym=(Fraction(a[1])+b[1])/2
      rec(a,(xm,ym),d+1,si); rec((xm,ym),b,d+1,si)
    for si,(a,b) in enumerate(base):
      # initial subdivisions by x or y
      dist=max(abs(b[0]-a[0]),abs(b[1]-a[1])); k=1
      while Fraction(dist,k)>h: k*=2
      for j in range(k):
       t0=Fraction(j,k);t1=Fraction(j+1,k)
       aa=(a[0]+t0*(b[0]-a[0]),a[1]+t0*(b[1]-a[1]))
       bb=(a[0]+t1*(b[0]-a[0]),a[1]+t1*(b[1]-a[1]))
       rec(aa,bb,0,si)
    sm=arb(0)
    for si,a,b,qarg,img in allq: sm += qarg
    return allq,failures,sm,n
def a1_cert():
    p=arb.pi()
    return 2*p*(-p).exp()*(15*p-4*p*p-arb(15)/2)

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--certificate')
 args=parser.parse_args()
 R=94
 t=time.time(); allq,fail,sm,n=cert(R=R)
 print('time',time.time()-t,'calls',n,'accepted',len(allq),'fail',len(fail))
 print('K0=F1(0)',F_point(0))
 print('a1=phi1_prime(0+)',a1_cert())
 print('F1(R)',F_point(R))
 print('F1(0).real interval',F_point(0).real)
 print('F1(R).real interval',F_point(R).real)
 print('sum',sm,'div2pi',sm/(2*pi))
 for z in fail[:8]: print('FAIL',z[0],z[1],z[2],z[3],z[5])
 print('minabs accepted',min(float(x[4].abs_lower()) for x in allq))
 w=sm/(2*pi)
 assert not fail
 assert F_point(0).real.lower() > 0
 assert F_point(R).real.upper() < 0
 assert w.lower() > arb(3)/2 and w.upper() < arb(5)/2
 assert min(x[4].abs_lower() for x in allq) > 0
 print('ASSERTIONS: boundary-zero exclusion, winding=2, and opposite real signs certified')
 if args.certificate:
  def point_data(p):
   return [str(Fraction(v)) for v in p]
  records=[]
  for si,a,b,arg,img in allq:
   records.append({'side':si,'a':point_data(a),'b':point_data(b),
                   'argument_increment':str(arg),
                   'image_real':str(img.real),'image_imag':str(img.imag),
                   'image_abs_lower':str(img.abs_lower())})
  # Exact rational adjacency closes the complete counterclockwise contour.
  assert all(allq[j][2] == allq[(j+1)%len(allq)][1]
             for j in range(len(allq)))
  result={'function':'Fourier transform of the even first theta head Phi_1',
          'R':R,'strip_half_width':'1/2','orientation':'counterclockwise',
          'precision_bits':ctx.prec,'python_flint_version':flint.__version__,
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'segment_count':len(allq),'failure_count':len(fail),
          'winding_interval':str(w),'winding_integer':2,
          'F_at_0':str(F_point(0)),'F_at_R':str(F_point(R)),
          'moment_upper_bounds':[str(v) for v in MOMENTS],
          'minimum_image_abs_lower':str(min(x[4].abs_lower() for x in allq)),
          'conclusion':'Exactly two real simple zeros in the rectangle; combine with the separate R=94 tail certificate for the closed strip.',
          'segments':records}
  Path(args.certificate).write_text(json.dumps(result,indent=2)+'\n')
