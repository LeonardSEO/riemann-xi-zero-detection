"""Exacte rationale controles voor Stelling 155.7; geen nulcertificaat.

Gebruikt uitsluitend gehele getallen en fractions.Fraction. De analytische
bewijsstappen en de externe digamma-invoer staan in route155 sectie12.
"""

from fractions import Fraction as Q
from math import factorial


def exp_lower(x: Q, degree: int) -> Q:
    return sum((x**k / factorial(k) for k in range(degree + 1)), Q(0))


def exp_upper(x: Q, degree: int) -> Q:
    ratio = x / (degree + 2)
    assert 0 < ratio < 1
    return exp_lower(x, degree) + x ** (degree + 1) / factorial(degree + 1) / (1 - ratio)


# gamma < H_30 - log(30) < 4 - 17/5 = 3/5.
assert sum((Q(1, k) for k in range(1, 31)), Q(0)) < 4
assert exp_upper(Q(17, 5), 20) < 30

# pi=4*int_0^1(1+x^2)^(-1)dx: eight geometric-series terms give pi>3.
assert 4 * sum((Q((-1) ** k, 2 * k + 1) for k in range(8)), Q(0)) > 3
# Polynomial division gives int_0^1 x^4(1-x)^4/(1+x^2)dx=22/7-pi>0.
assert Q(11, 7) < Q(393, 250)
assert exp_lower(Q(347, 500), 8) > 2
assert exp_lower(Q(573, 500), 8) > Q(22, 7)
ca_upper = Q(3, 5) + Q(393, 250) + 3 * Q(347, 500) + Q(573, 500)
assert ca_upper == Q(27, 5)

# Monotone positive digamma series: a(8) > -27/5 + these 16 terms.
series16 = sum(
    (Q(16) / (Q(4 * n + 1, 4) * (Q((4 * n + 1) ** 2, 16) + 16)) for n in range(16)),
    Q(0),
)
assert series16 - ca_upper > Q(9, 50)

# R=3/20, H=8; cosh(R/2) <= 1/(1-(R/2)^2/2).
radius, height, cosh_upper = Q(3, 20), Q(8), Q(1003, 1000)
assert Q(3200, 3191) < cosh_upper
i1 = height**5 / 20 + cosh_upper * height**3 / 24 + cosh_upper**2 * height / 64
i2 = height**7 / 252 + cosh_upper * height**5 / 360 + cosh_upper**2 * height**3 / 1728
band_bound = (2 * radius**5 * i1 / 5 + 2 * radius**7 * i2 / 7) / 3
assert band_bound < Q(19, 1000)

# a(8)(1-band_bound)-C_a*band_bound > 7/100.
assert Q(9, 50) - (Q(9, 50) + ca_upper) * Q(19, 1000) > Q(7, 100)

print("PASS: alle rationale inequalities exact; q_B >= 7/100 ||f||_2^2 op de gegeven poolvrije core.")
