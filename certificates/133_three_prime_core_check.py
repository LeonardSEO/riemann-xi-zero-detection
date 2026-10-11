#!/usr/bin/env python3
"""Candidate full-Schur certificate, not a proof until its main passes.

The actual shifts are log2,log3,log4.  The script reuses the independently
audited polynomial/interval identities from133_two_prime_core_check.py.
Every integral is exact polynomial integration; no floating point occurs.
"""
import argparse
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import importlib.util
import json

spec = importlib.util.spec_from_file_location(
    'two_prime', Path(__file__).with_name('133_two_prime_core_check.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.interval_module.SCALE = m.SCALE = 10**400
I = m.I
zero = m.zero = I.rat(0)
one = m.one = I.rat(1)
m.R = R = F(7,10)
m.N = N = 100
m.P = P = 64
m.J = J = 256
m.DEG = DEG = 64
m.PENALTY = PENALTY = 100000
m.KERNEL_ERROR = KERNEL_ERROR = 7*R**(DEG+1)/(1-R)
SCHUR_FLOOR = F(1, 10**8)


def lower30(value):
    return str(F(value.lo//10**370, 10**30))


def exp_negative(z):
    q = I.coerce(z)/512
    assert 0 <= q.lo <= q.hi <= m.SCALE
    out = I.rat(F(1, factorial(256)))
    for k in range(255,-1,-1):
        out = out*q+I.rat(F((-1)**k, factorial(k)))
    out = out.widen(F(1, factorial(257)))
    assert out.lo > 0
    for unused in range(9): out = out.sq()
    return out


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
    assert a[0]-osc > F(1,8)
    assert a[0]+osc < F(3,8)
    return a


# Cache the exact monomial-to-Legendre identity; interval arithmetic is
# unchanged, only repeated rational coefficient construction is avoided.
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


def main(output=None):
    a0,h2,h3,w2,w3=m.constants()
    h4,w4=2*h2,h2/2
    assert h2.hi < I.rat(R).lo
    assert h3.lo > I.rat(R).hi
    assert h4.hi < I.rat(2*R).lo
    # No n>=5 interaction can meet this support.
    assert m.atanh_log(I.rat(2,3)).lo > I.rat(2*R).hi
    # The n=2 graph has three-node chains; its norm is sqrt2*w2=log2.
    assert (h2+w3+w4).hi < I.rat(7,4).lo
    ps=m.legendre_polynomials(N)
    kernel=kernel_polynomial()
    cpoly=[F(0)]*25;spoly=[F(0)]*26
    for j in range(13):
        cpoly[2*j]=(R/2)**(2*j)/factorial(2*j)
        spoly[2*j+1]=(R/2)**(2*j+1)/factorial(2*j+1)
    leak=F(PENALTY)*2*R*((2*(R/2)**26/factorial(26))**2
                         +(2*(R/2)**27/factorial(27))**2)
    vp=[zero]*(2*P+1)
    for j in range(J):
        lam=F(4*j+1,2)
        decay=exp_negative(I.rat(lam*R))
        for p in range(P+1):
            vp[2*p]+=decay*I.rat(2*(lam*R)**(2*p)/(lam*factorial(2*p)))
    v0=vp[0];vp[0]=zero
    assert sum(v.hi for v in vp)<3*m.SCALE
    hn=[F(0)]
    for j in range(1,N+1):hn.append(hn[-1]+F(1,j))
    tailgap=I.rat(hn[N])+a0+v0-I.rat(7,4)
    assert tailgap.lo > m.SCALE
    images=[]
    for p in ps:
        img=[I.rat(v) for v in m.kernel_action(p,kernel)]
        img=m.iadd(img,m.imul(vp,[I.rat(v) for v in p]))
        penalty=m.fadd([PENALTY*m.f_inner(cpoly,p)*v for v in cpoly],
                       [PENALTY*m.f_inner(spoly,p)*v for v in spoly])
        images.append(m.iadd(img,[I.rat(v) for v in penalty]))
    shifts=[h2/I.rat(R),h3/I.rat(R),h4/I.rat(R)]
    weights=[w2,w3,w4]
    ends=[-one,one]
    pluscuts=[];minuscuts=[]
    for shift in shifts:
        pluscuts.append(one-shift);minuscuts.append(-one+shift)
    ends+=pluscuts+minuscuts
    ends.sort(key=lambda x:x.lo+x.hi)
    for a,b in zip(ends,ends[1:]):assert a.hi<b.lo
    def location(c):
        return next(j for j,e in enumerate(ends) if e is c)
    plusidx=[location(c) for c in pluscuts]
    minusidx=[location(c) for c in minuscuts]
    bpieces=[];wpieces=[];pieceweights=[]
    for piece,(a,b) in enumerate(zip(ends,ends[1:])):
        mid,half=(a+b)/2,(b-a)/2
        bb=[legendre_coefficients(m.compose([I.rat(v) for v in p],mid,half)) for p in ps]
        wi=[]
        for j,p in enumerate(ps):
            img=images[j][:]
            for q,w,left,right in zip(shifts,weights,plusidx,minusidx):
                if piece<left:
                    pp=m.compose([I.rat(v) for v in p],q,one)
                    img=m.iadd(img,[-w*v for v in pp])
                if piece>=right:
                    pp=m.compose([I.rat(v) for v in p],-q,one)
                    img=m.iadd(img,[-w*v for v in pp])
            wi.append(legendre_coefficients(m.compose(img,mid,half)))
        wt=[I.rat(R)*(b-a)/I.rat(2*k+1) for k in range(max(map(len,wi)))]
        bpieces.append(bb);wpieces.append(wi);pieceweights.append(wt)
    matrix=[[zero]*N for unused in range(N)]
    for i in range(N):
        for j in range(N):
            val=zero
            for bp,wp,wt in zip(bpieces,wpieces,pieceweights):
                val+=m.dot(bp[i],wp[j],wt)
            matrix[i][j]=val
    for i in range(N):
        for j in range(i):
            matrix[i][j]=matrix[j][i]=(matrix[i][j]+matrix[j][i])/2
    norms=[F(2*R,2*n+1) for n in range(N)]
    gram=[[zero]*N for unused in range(N)]
    for bp,wp,wt in zip(bpieces,wpieces,pieceweights):
        residuals=[]
        for j in range(N):
            rr=wp[j][:]+[zero]*max(0,N-len(wp[j]))
            for n in range(N):
                factor=matrix[n][j]/I.rat(norms[n])
                for k,v in enumerate(bp[n]):rr[k]-=factor*v
            residuals.append(rr)
        for i in range(N):
            for j in range(i+1):
                gram[i][j]+=m.dot(residuals[i],residuals[j],wt)
    schur=[[zero]*N for unused in range(N)]
    for i in range(N):
        for j in range(i+1):
            val=matrix[i][j]-gram[i][j]
            if i==j:val+=(I.rat(hn[i])+a0+v0-I.rat(SCHUR_FLOOR))*I.rat(norms[i])
            schur[i][j]=schur[j][i]=val
    pivots=[];lower=[[zero]*N for unused in range(N)]
    for i in range(N):
        pivot=schur[i][i]
        for k in range(i):pivot-=lower[i][k].sq()*pivots[k]
        assert pivot.lo>0,('nonpositive pivot',i,lower30(pivot))
        pivots.append(pivot);lower[i][i]=one
        for j in range(i+1,N):
            value=schur[j][i]
            for k in range(i):value-=lower[j][k]*lower[i][k]*pivots[k]
            lower[j][i]=value/pivot
    coercivity=SCHUR_FLOOR/F(73)-leak
    assert coercivity>F(1,10**11)
    result={'status':'PASS','arithmetic':'outwardinteger intervals scale10^400; exact polynomial integrals',
            'support_diameter':'7/5','actual_prime_powers':[2,3,4],
            'low_legendre_modes':N,'gamma_heads':J,'outside_Taylor_degree':2*P,
            'smooth_kernel_degree':DEG,'smooth_kernel_error':str(KERNEL_ERROR),
            'tail_gap_lower':lower30(tailgap),'Schur_floor':str(SCHUR_FLOOR),
            'LDL_pivot_lower_bounds':[lower30(v) for v in pivots],
            'pool_leakage_upper':str(leak),'coercivity_lower':str(coercivity),
            'claimed_coercivity':'1/100000000000'}
    if output:Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output')
    main(parser.parse_args().output)
