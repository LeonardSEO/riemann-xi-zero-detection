"""Nieuwe exacte constantcontrole voor route155 Stelling155.10.

Geen nulcertificaat. De reeds bewezen bounds C_a<27/5 en
log(2)<347/500 zijn invoer uit155.7; deze controle herhaalt
uitsluitend de nieuwe digamma-, band- en prime-edgeconstants.
"""

from fractions import Fraction as Q

radius, height, cosh_upper = Q(3, 40), Q(12), Q(1001, 1000)
assert 1 / (1 - (radius / 2) ** 2 / 2) < cosh_upper

# Positieve digamma-difference: H=12, nu met eerste32 terms.
series32 = sum(
    (Q(36) / (Q(4 * n + 1, 4) * (Q((4 * n + 1) ** 2, 16) + 36)) for n in range(32)),
    Q(0),
)
assert series32 - Q(27, 5) > Q(3, 5)

i1 = height**5 / 20 + cosh_upper * height**3 / 24 + cosh_upper**2 * height / 64
i2 = height**7 / 252 + cosh_upper * height**5 / 360 + cosh_upper**2 * height**3 / 1728
band_bound = (2 * radius**5 * i1 / 5 + 2 * radius**7 * i2 / 7) / 3
assert band_bound < Q(1, 200)

# e^(3/20)<=1/(1-3/20)=20/17<3/2: n=3 en groter ontbreken.
assert 1 / (1 - Q(3, 20)) < Q(3, 2)
# sqrt(2)>7/5, dus log(2)/sqrt(2)<347/700<1/2.
assert Q(7, 5) ** 2 < 2
assert Q(347, 700) < Q(1, 2)

# Twee cellen: low-bandnorm²<2/200 en genuine primecompressionnorm<1/2.
assert Q(3, 5) - (Q(3, 5) + Q(27, 5)) * Q(1, 100) - Q(1, 2) == Q(1, 25)

print("PASS: exacte nieuwe constants; twee lokaal poolvrije cellen met afstand log2 hebben q_A >= 1/25 ||f||_2^2.")

# Nieuwe globale twee-cel-interpolatie: iedere n>=2, r<=1/100, H=16.
radius, height = Q(1, 100), Q(16)
cosh_small = 1 / (1 - radius**2 / 8)
# coth(log2/4)=3+2sqrt2<6, tanh<1; de derivativecost is <4.
assert Q(3, 2) ** 2 > 2
assert (Q(7, 2) + radius / 2) * cosh_small < 4

series32 = sum(
    (Q(64) / (Q(4 * n + 1, 4) * (Q((4 * n + 1) ** 2, 16) + 64)) for n in range(32)),
    Q(0),
)
assert series32 - Q(27, 5) > Q(17, 20)
band_global = 4 * radius**3 * ((height + 4) ** 3 - 4**3) / 27
assert band_global < Q(1, 500)
# max_x log(x)/sqrt(x)=2/e<3/4, want e>8/3 from five positive terms.
assert Q(1) + 1 + Q(1, 2) + Q(1, 6) + Q(1, 24) > Q(8, 3)
assert Q(17, 20) - (Q(17, 20) + Q(27, 5)) / 500 - Q(3, 4) == Q(7, 80)

print("PASS: globale pool-null twee-cellen op afstand log(n), radius<=min(1/100,1/(4n)), hebben q_A >= 7/80 ||f||_2^2.")
