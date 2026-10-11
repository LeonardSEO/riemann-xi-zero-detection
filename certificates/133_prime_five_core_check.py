#!/usr/bin/env python3
"""Candidate full core certificate with every shift through n=5.

No positivity theorem is claimed until main() passes.  --constants-only
checks only the new positive-weight and analytic-tail constants.  All
integrals in main() are exact piecewise-polynomial integrals, enclosed
with outward integer intervals; no floating point or quadrature occurs.
"""
import argparse
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import importlib.util
import json
import sys
import hashlib

spec = importlib.util.spec_from_file_location(
    'two_prime', Path(__file__).with_name('133_two_prime_core_check.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.interval_module.SCALE = m.SCALE = 10**400
I = m.I
zero = m.zero = I.rat(0)
one = m.one = I.rat(1)
m.R = R = F(17,20)
m.N = N = 100
m.DEG = DEG = 96
m.PENALTY = PENALTY = 100000
J = 384
OUT_DEG = 96
ALPHA = F(8,5)
PRIME_BOUND = F(9,5)
OUT_ERROR = F(1,10**22)
KERNEL_ERROR = m.KERNEL_ERROR = 17*(2*R/3)**(DEG+1)/(1-2*R/3)
old_atanh_log = m.atanh_log
m.atanh_log = lambda q: old_atanh_log(q,1000)


def lower30(value):
    return str(F(value.lo//10**370,10**30))


def upper30(value):
    return str(F(-((-value.hi)//10**370),10**30))


def exp_negative(z):
    q = I.coerce(z)/2048
    assert 0 <= q.lo <= q.hi <= m.SCALE
    out = I.rat(F(1,factorial(256)))
    for k in range(255,-1,-1):
        out = out*q+I.rat(F((-1)**k,factorial(k)))
    out = out.widen(F(1,factorial(257)))
    assert out.lo > 0
    for unused in range(11):
        out = out.sq()
    return out


def cosh(z):
    z = I.coerce(z)
    if z.hi <= 0:
        z = -z
    elif z.lo < 0:
        upper = I(0,max(-z.lo,z.hi))
        e = exp_negative(upper)
        c = (e+one/e)/2
        return I(m.SCALE,c.hi)
    e = exp_negative(z)
    return (e+one/e)/2


def kernel_polynomial():
    maxn=DEG+1
    u=[F(1,factorial(n+1)) if not n%2 else F(0) for n in range(maxn+1)]
    inv=[F(1)]
    for n in range(1,maxn+1):
        inv.append(-sum(u[j]*inv[n-j] for j in range(1,n+1)))
    exp=[F(1,2**n*factorial(n)) for n in range(maxn+1)]
    a=[sum(exp[j]*inv[n+1-j] for j in range(n+2))/2 for n in range(DEG+1)]
    a[0]-=KERNEL_ERROR
    osc=sum(abs(a[n])*(2*R)**n for n in range(1,DEG+1))
    assert a[0]-osc > F(1,16)
    assert a[0]+osc < F(7,16)
    return a


def arithmetic_constants():
    pi=16*m.atan(F(1,5))-4*m.atan(F(1,239))
    h2=m.atanh_log(I.rat(1,3))
    h3=m.atanh_log(I.rat(1,2))
    logpi=m.atanh_log((pi-one)/(pi+one))
    b=m.bernoulli(64)
    gamma=I.rat(sum(F(1,k) for k in range(1,257))-F(1,512))-8*h2
    gamma+=I.rat(sum(b[2*k]/F(2*k*256**(2*k)) for k in range(1,32)))
    gamma=gamma.widen(abs(b[64])/F(64*256**64))
    a0=-(gamma+pi/2+3*h2+logpi)
    roots=[]
    for n in (2,3):
        root=isqrt(n*m.SCALE*m.SCALE)
        roots.append(I(root,root+1))
    return a0,h2,h3,h2/roots[0],h3/roots[1]


def smooth_high_constants(kernel):
    """Legendre-Sturm bound for the full smooth integral-operator tail."""
    from math import comb
    row=[F(0)]*(len(kernel)+1)
    for n,c in enumerate(kernel):
        for k in range(0,n+2,2):
            row[k]+=2*c*R**(n+1)*F(comb(n+1,k),n+1)
    row_floor=row[0]-sum(abs(v) for v in row[2:])
    c1=sum(n*abs(kernel[n])*(2*R)**(n-1)
           for n in range(1,len(kernel)))
    c2=sum(n*(n-1)*abs(kernel[n])*(2*R)**(n-2)
           for n in range(2,len(kernel)))
    sturm_norm=2*R*c2+4*c1+2*abs(kernel[1])
    high_error=R**2*sturm_norm/F(N*(N+1))
    assert row_floor>F(3,10)
    assert c1<F(1,4) and c2<F(1,2) and sturm_norm<2
    return row_floor,c1,c2,sturm_norm,high_error


coefficients={}
def legendre_coefficients(poly):
    out=[]
    for k in range(len(poly)):
        v=zero
        for n in range(k,len(poly),2):
            key=(n,k)
            if key not in coefficients:
                coefficients[key]=I.rat(F((2*k+1)*factorial(n),
                    m.double_factorial(n-k)*m.double_factorial(n+k+1)))
            v+=poly[n]*coefficients[key]
        out.append(v)
    return out


def constants():
    sin_lower=sum((-1)**k*F(3**(2*k+1),factorial(2*k+1)) for k in range(10))
    assert sin_lower > F(7,50)
    e=F(3,2)
    exp_upper=sum(e**k/factorial(k) for k in range(21))
    exp_upper+=e**21/factorial(21)/(1-e/22)
    assert exp_upper < F(9,2)
    assert F(9,2)/(2*F(7,50))+F(1,6) < 17
    assert F(4*J-3,2)*2*R < 2048
    a0,h2,h3,w2,w3=arithmetic_constants()
    h4,w4=2*h2,h2/2
    h5=m.atanh_log(I.rat(2,3))
    root5=isqrt(5*m.SCALE*m.SCALE)
    w5=h5/I(root5,root5+1)
    assert h5.hi < I.rat(2*R).lo
    assert m.atanh_log(I.rat(5,7)).lo > I.rat(2*R).hi
    assert (2*h2).hi < I.rat(2*R).lo < (3*h2).lo
    shifts=[h/I.rat(R) for h in (h2,h3,h4,h5)]
    weights=[w2,w3,w4,w5]
    pluscuts=[one-q for q in shifts]
    minuscuts=[-one+q for q in shifts]
    prime_ends=[-one,one]+pluscuts+minuscuts
    prime_ends.sort(key=lambda z:z.lo+z.hi)
    for a,b in zip(prime_ends,prime_ends[1:]):
        assert a.hi < b.lo
    def location(c,ends):
        return next(k for k,z in enumerate(ends) if z is c)
    plusidx=[location(c,prime_ends) for c in pluscuts]
    minusidx=[location(c,prime_ends) for c in minuscuts]
    ratio_bounds=[]
    for piece,(a,b) in enumerate(zip(prime_ends,prime_ends[1:])):
        for x in (a,b):
            value=zero
            for q,w,lp,lm in zip(shifts,weights,plusidx,minusidx):
                if piece<lp:
                    value+=w*cosh(I.rat(ALPHA*R)*(x+q))/cosh(I.rat(ALPHA*R)*x)
                if piece>=lm:
                    value+=w*cosh(I.rat(ALPHA*R)*(x-q))/cosh(I.rat(ALPHA*R)*x)
            assert value.hi < I.rat(PRIME_BOUND).lo
            ratio_bounds.append(upper30(value))
    ends=prime_ends+[zero]
    for k in range(1,7):
        x=one-I.rat(F(1,2**k)/R)
        ends += [x,-x]
    ends.sort(key=lambda z:z.lo+z.hi)
    for a,b in zip(ends,ends[1:]):
        assert a.hi < b.lo
        if a.lo == -m.SCALE or b.hi == m.SCALE:
            umax=I.rat(F(4*J-3,2)*R)*(b-a)
            assert umax.hi <= 12*m.SCALE
        else:
            assert ((b-a)/(one-a)).hi < I.rat(51,100).lo
            assert ((b-a)/(one+b)).hi < I.rat(51,100).lo
    n=OUT_DEG+1
    q=F(51,100)
    interior=4*J*q**n/(1-q)
    endpoint=4*J*F(12)**n/factorial(n)/(1-F(12,n+1))
    assert interior < OUT_ERROR
    assert endpoint < OUT_ERROR
    kernel=kernel_polynomial()
    harmonic=sum(F(1,k) for k in range(1,J))
    assert 2+harmonic/2+OUT_ERROR < 6
    return a0,shifts,weights,ends,pluscuts,minuscuts,kernel,ratio_bounds


def outside_local(a,b):
    """Positive Taylor lower potential in the piece coordinate u in [-1,1]."""
    half=I.rat(R)*(b-a)/2
    out=[zero]*(OUT_DEG+1)
    for j in range(J):
        lam=F(4*j+1,2)
        left=exp_negative(I.rat(lam*R)*(one-a))
        right=exp_negative(I.rat(lam*R)*(one+b))
        powers=[one]
        v=I.rat(lam)*half
        for k in range(1,OUT_DEG+1):
            powers.append(powers[-1]*v/k)
        sums=[];value=zero
        for p in powers:
            value+=p
            sums.append(value)
        for k in range(OUT_DEG+1):
            coef=(left+right if not k%2 else left-right)/I.rat(lam)
            out[k]+=coef*powers[k]*sums[OUT_DEG-k]
    return out


def main(output=None,constants_only=False):
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a0,shifts,weights,ends,pluscuts,minuscuts,kernel,ratio_bounds=constants()
    if constants_only:
        result={'status':'CONSTANTS_PASS','full_core_status':'NOT_CHECKED',
                'prime_ratio_upper_bounds':ratio_bounds,
                'prime_bound':str(PRIME_BOUND),'kernel_error':str(KERNEL_ERROR),
                'outside_error':str(OUT_ERROR),
                'smooth_high_constants':[str(v) for v in smooth_high_constants(kernel)]}
        print(json.dumps(result,indent=2))
        return
    print('Candidate: constants and new analytic bounds passed.',file=sys.stderr,flush=True)
    ps=m.legendre_polynomials(N)
    cpoly=[F(0)]*25;spoly=[F(0)]*26
    for j in range(13):
        cpoly[2*j]=(R/2)**(2*j)/factorial(2*j)
        spoly[2*j+1]=(R/2)**(2*j+1)/factorial(2*j+1)
    leak=F(PENALTY)*2*R*((2*(R/2)**26/factorial(26))**2
                         +(2*(R/2)**27/factorial(27))**2)
    v0=zero
    for j in range(J):
        lam=F(4*j+1,2)
        v0+=2*exp_negative(I.rat(lam*R))/I.rat(lam)
    base=v0-I.rat(OUT_ERROR)
    hn=[F(0)]
    for j in range(1,N+1):
        hn.append(hn[-1]+F(1,j))
    row_floor,c1,c2,sturm_norm,high_error=smooth_high_constants(kernel)
    tailgap=I.rat(hn[N])+a0+base-I.rat(PRIME_BOUND)+I.rat(row_floor-high_error)
    assert tailgap.lo > m.SCALE
    images=[]
    for p in ps:
        img=[I.rat(v) for v in m.kernel_action(p,kernel)]
        penalty=m.fadd([PENALTY*m.f_inner(cpoly,p)*v for v in cpoly],
                       [PENALTY*m.f_inner(spoly,p)*v for v in spoly])
        images.append(m.iadd(img,[I.rat(v) for v in penalty]))
    def location(c):
        return next(j for j,e in enumerate(ends) if e is c)
    plusidx=[location(c) for c in pluscuts]
    minusidx=[location(c) for c in minuscuts]
    bpieces=[];wpieces=[];pieceweights=[]
    for piece,(a,b) in enumerate(zip(ends,ends[1:])):
        if a.hi <= 0 and b.hi <= 0:
            continue
        assert a.lo >= 0
        mid,half=(a+b)/2,(b-a)/2
        polys=[m.compose([I.rat(v) for v in p],mid,half) for p in ps]
        bb=[legendre_coefficients(p) for p in polys]
        outside=outside_local(a,b)
        outside[0]-=base
        wi=[]
        for j,p in enumerate(ps):
            img=m.compose(images[j],mid,half)
            img=m.iadd(img,m.imul(outside,polys[j]))
            for q,w,lp,lm in zip(shifts,weights,plusidx,minusidx):
                if piece<lp:
                    pp=m.compose([I.rat(v) for v in p],mid+q,half)
                    img=m.iadd(img,[-w*v for v in pp])
                if piece>=lm:
                    pp=m.compose([I.rat(v) for v in p],mid-q,half)
                    img=m.iadd(img,[-w*v for v in pp])
            wi.append(legendre_coefficients(img))
        # Reflection doubles every equal-parity product; mixed parity is exactly zero.
        wt=[I.rat(2*R)*(b-a)/I.rat(2*k+1) for k in range(max(map(len,wi)))]
        bpieces.append(bb);wpieces.append(wi);pieceweights.append(wt)
        print('Candidate: positive piece '+str(piece)+' built.',file=sys.stderr,flush=True)
    matrix=[[zero]*N for unused in range(N)]
    for i in range(N):
        for j in range(i%2,i+1,2):
            val=zero
            for bp,wp,wt in zip(bpieces,wpieces,pieceweights):
                val+=m.dot(bp[i],wp[j],wt)
            matrix[i][j]=matrix[j][i]=val
    print('Candidate: low matrix built.',file=sys.stderr,flush=True)
    norms=[F(2*R,2*n+1) for n in range(N)]
    gram=[[zero]*N for unused in range(N)]
    for bp,wp,wt in zip(bpieces,wpieces,pieceweights):
        residuals=[]
        for j in range(N):
            rr=wp[j][:]+[zero]*max(0,N-len(wp[j]))
            for n in range(j%2,N,2):
                factor=matrix[n][j]/I.rat(norms[n])
                for k,v in enumerate(bp[n]):
                    rr[k]-=factor*v
            residuals.append(rr)
        for i in range(N):
            for j in range(i%2,i+1,2):
                gram[i][j]+=m.dot(residuals[i],residuals[j],wt)
    print('Candidate: full residual Gram built.',file=sys.stderr,flush=True)
    schur=[[zero]*N for unused in range(N)]
    for i in range(N):
        for j in range(i%2,i+1,2):
            val=matrix[i][j]-gram[i][j]
            if i==j:
                val+=(I.rat(hn[i])+a0+base)*I.rat(norms[i])
            schur[i][j]=schur[j][i]=val
    if output:
        snapshot={'status':'CANDIDATE_MATRIX_NOT_A_CERTIFICATE','source_sha256':source_hash,
                  'scale':str(m.SCALE),'norms':[str(v) for v in norms],
                  'Schur_matrix':[[[str(v.lo),str(v.hi)] for v in row] for row in schur]}
        Path(output).with_suffix('.candidate_matrix.json').write_text(json.dumps(snapshot))
    pivots=[];lower=[[zero]*N for unused in range(N)]
    for i in range(N):
        pivot=schur[i][i]
        for k in range(i%2,i,2):
            pivot-=lower[i][k].sq()*pivots[k]
        if pivot.lo<=0:
            failure={'status':'NOT_CERTIFIED','source_sha256':source_hash,'pivot_index':i,
                     'pivot_lower':lower30(pivot),'pivot_upper':upper30(pivot),
                     'classification':'negative finite Schur lower bound' if pivot.hi<0
                                      else 'interval positivity unresolved',
                     'previous_pivot_lower_bounds':[lower30(v) for v in pivots]}
            if output:
                Path(output).with_suffix('.failure.json').write_text(json.dumps(failure,indent=2)+'\n')
            raise AssertionError(failure)
        pivots.append(pivot);lower[i][i]=one
        for j in range(i+2,N,2):
            value=schur[j][i]
            for k in range(i%2,i,2):
                value-=lower[j][k]*lower[i][k]*pivots[k]
            lower[j][i]=value/pivot
    inv_lower=[[zero]*N for unused in range(N)]
    for i in range(N):
        inv_lower[i][i]=one
        for j in range(i%2,i,2):
            value=zero
            for k in range(j,i,2):
                value-=lower[i][k]*inv_lower[k][j]
            inv_lower[i][j]=value
    trace=zero
    for i in range(N):
        for k in range(i,N,2):
            trace+=I.rat(norms[i])*inv_lower[k][i].sq()/pivots[k]
    assert trace.hi>0
    schur_floor=F(m.SCALE,trace.hi)
    coercivity=schur_floor/F(201)-leak
    assert coercivity>0,('pool leakage exceeds inverse-tracefloor',str(coercivity))
    claimed=F(1)
    while claimed>=coercivity:
        claimed/=10
    result={'status':'PASS','arithmetic':'outwardinteger intervals scale10^400; exact polynomial integrals',
            'support_diameter':'17/10','actual_prime_powers':[2,3,4,5],
            'prime_weight':'cosh(8*x/5)','prime_bound':str(PRIME_BOUND),
            'prime_ratio_upper_bounds':ratio_bounds,'low_legendre_modes':N,
            'gamma_heads':J,'outside_piece_degree':OUT_DEG,
            'outside_error':str(OUT_ERROR),'positive_piece_count':len(bpieces),
            'smooth_kernel_degree':DEG,'smooth_kernel_error':str(KERNEL_ERROR),
            'smooth_row_floor':str(row_floor),'smooth_derivative_bounds':[str(c1),str(c2)],
            'smooth_Sturm_norm':str(sturm_norm),'smooth_high_error':str(high_error),
            'source_sha256':source_hash,'tail_gap_lower':lower30(tailgap),
            'LDL_subtracted_floor':'0','inverse_trace_upper':str(F(trace.hi,m.SCALE)),
            'Schur_floor':str(schur_floor),
            'LDL_pivot_lower_bounds':[lower30(v) for v in pivots],
            'pool_leakage_upper':str(leak),'coercivity_lower':str(coercivity),
            'claimed_coercivity':str(claimed)}
    if output:
        Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output')
    parser.add_argument('--constants-only',action='store_true')
    args=parser.parse_args()
    main(args.output,args.constants_only)
