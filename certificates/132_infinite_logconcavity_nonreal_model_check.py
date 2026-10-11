#!/usr/bin/env python3
# Exact rational check for the auxiliary partial-theta model in route132 section46.
# This is not a zero certificate for the Riemann xi or zeta function.
from fractions import Fraction as Q

c = (Q(-71, 32), Q(81, 128))
radius = Q(1, 16)
bound_radius = Q(12, 5)
N = 8

def add(a, b):
    return a[0] + b[0], a[1] + b[1]

def mul(a, b):
    return a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]

def scale(a, q):
    return a[0]*q, a[1]*q

powers = [(Q(1), Q(0))]
for n in range(1, N + 1):
    powers.append(mul(powers[-1], c))
coeff = lambda n: Q(1, 3**(n*(n-1)//2))
P = (Q(0), Q(0))
D = (Q(0), Q(0))
M = Q(0)
for n in range(N + 1):
    P = add(P, scale(powers[n], coeff(n)))
    if n:
        D = add(D, scale(powers[n-1], n*coeff(n)))
    if n >= 2:
        M += n*(n-1)*coeff(n)*bound_radius**(n-2)

# For n>=N+1, each successive term has ratio bounded as follows.
n = N + 1
q0 = bound_radius / 3**n
q1 = Q(n+1, n)*q0
q2 = Q(n+1, n-1)*q0
assert q0 < 1 and q1 < 1 and q2 < 1
E0 = coeff(n)*bound_radius**n / (1-q0)
E1 = n*coeff(n)*bound_radius**(n-1) / (1-q1)
E2 = n*(n-1)*coeff(n)*bound_radius**(n-2) / (1-q2)

# The entire circle is in |z|<bound_radius and strictly above the real axis.
assert c[0]**2 + c[1]**2 < (bound_radius-radius)**2
assert c[1]-radius == Q(73,128) > 0

# Certified enclosures for the full infinite series, not just P_N.
assert abs(P[0]) + abs(P[1]) + E0 < Q(1,1000)
assert D[1] - E1 > Q(3,20)
assert M + E2 < Q(4,3)

# Linear Rouche comparison L(z)=f(c)+f'(c)(z-c).
linear_lower = Q(3,20)*radius-Q(1,1000)
remainder_upper = Q(4,3)*radius**2/2
assert linear_lower > remainder_upper > 0
assert Q(1,1000) < Q(3,20)*radius

print('PASS: exact-rational full-series tail and contour bounds')
print('|f(c)|<1/1000; Im fprime(c)>3/20; sup |fsecond|<4/3')
print('Rouche: exactly one zero in |z-c|<1/16, with Im z>73/128')
print('Auxiliary partial-theta model only; no xi/zeta zero is asserted.')
