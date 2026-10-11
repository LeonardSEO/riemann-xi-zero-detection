# Frozen specification: final smoothed-Stieltjes implication

Date frozen: 2026-07-24.

Status: immutable common gate.  Changing the packet, Fourier sign,
Laplace factor, active/reference pair, polar transport, deleted mode,
spectral scale, counting-shift sign, safe interval, or prime-power
endpoint restarts every branch and audit.

This phase uses the primary semilocal formulas of:

1. A. Connes and C. Consani, *Spectral triples and zeta-cycles*,
   Enseign. Math. 69 (2023), 93--148, DOI `10.4171/LEM/1049`;
2. A. Connes, C. Consani, and H. Moscovici,
   *Zeta Spectral Triples*, arXiv `2511.22755`;
3. the locally re-audited source file
   `../../tmp/pdfs/2511.22755v1-src/mc2arXiv.tex`.

Earlier project files are provenance aids, not substitutes for the
reconstruction below.  No zeta zero is supplied to an active operator or
reference operator.

## 1. Hilbert, Fourier, and Laplace conventions

The Hilbert inner product is antilinear in its first argument.  In
logarithmic coordinates,

\[
\widehat f(t)=\int_{\mathbb R}f(x)e^{-itx}\,dx,
\qquad
f(x)=\frac1{2\pi}\int_{\mathbb R}\widehat f(t)e^{itx}\,dt.
\tag{1.1}
\]

For this phase freeze

\[
Y=\frac34,\qquad
\eta=\frac14,\qquad
S=1,\qquad
I_{\eta,S}=\left[\frac34,1\right].
\tag{1.2}
\]

The bilateral transform uses

\[
\mathcal B[F](s)=\int_{\mathbb R}F(y)e^{-sy}\,dy,
\tag{1.3}
\]

and the one-sided factor-two transform uses

\[
\mathcal L_{2,Y}[F](s)
=2\int_Y^\infty e^{-2sy}F(y)\,dy.
\tag{1.4}
\]

## 2. Frozen packet

For \(n\ge1\), put

\[
a_n=2^{-n-1},
\qquad
b_{a_n}(x)=\frac{\mathbf1_{[-a_n,a_n]}(x)}{2a_n},
\tag{2.1}
\]

\[
g=*_{n\ge1}b_{a_n},
\qquad
h=\frac g{\|g\|_2}.
\tag{2.2}
\]

Then \(h\) is real, even, nonnegative, \(C_c^\infty\), supported exactly
on \([-1/2,1/2]\), and \(\|h\|_2=1\).  Freeze

\[
K_h=h*h,\qquad
L_h(w)=\int_{-1}^{1}K_h(r)e^{wr}\,dr,\qquad
M_h=\int_{\mathbb R}h(x)\,dx.
\tag{2.3}
\]

Thus

\[
w_h(t):=L_h(it)=|\widehat h(t)|^2\ge0,
\qquad
L_h(0)=M_h^2.
\tag{2.4}
\]

The normally convergent product is

\[
L_h(w)
=\|g\|_2^{-2}
\prod_{n\ge1}
\left(\frac{\sinh(2^{-n-1}w)}{2^{-n-1}w}\right)^2.
\tag{2.5}
\]

Its exact zero divisor is

\[
Z(L_h)=4\pi i\mathbb Z\setminus\{0\},
\qquad
\operatorname{ord}_{4\pi ik}L_h
=2\bigl(v_2(|k|)+1\bigr).
\tag{2.6}
\]

In particular,

\[
\Re w\ne0\Longrightarrow L_h(w)\ne0.
\tag{2.7}
\]

Both \(w_h\) and \(w_h'\) are Schwartz functions.

## 3. Complete semilocal form

For \(\lambda>1\), let

\[
\mathcal H_\lambda
=L^2([\lambda^{-1},\lambda],du/u),
\qquad
L_\lambda=2\log\lambda.
\tag{3.1}
\]

Functions are extended by zero.  The exact closed semilocal Weil form is

\[
\begin{aligned}
QW_\lambda(f,g)
={}&
\int_{\mathbb R}
\overline{\widehat f(t)}\,\widehat g(t)
\frac{2\theta'(t)}{2\pi}\,dt\\
&+\overline{\widehat f(i/2)}\,\widehat g(-i/2)
+\overline{\widehat f(-i/2)}\,\widehat g(i/2)\\
&-\sum_{\substack{p\ {\rm prime},\,m\ge1\\p^m\le\lambda^2}}
(\log p)p^{-m/2}
\left((f^**g)(p^m)+(f^**g)(p^{-m})\right),
\end{aligned}
\tag{3.2}
\]

where

\[
\theta(t)
=\operatorname{Im}\log\Gamma\left(\frac14+\frac{it}{2}\right)
-\frac t2\log\pi,
\tag{3.3}
\]

\[
2\theta'(t)
=\operatorname{Re}\frac{\Gamma'}{\Gamma}
\left(\frac14+\frac{it}{2}\right)-\log\pi.
\tag{3.4}
\]

Every prime power \(p^m\le\lambda^2\) occurs once with coefficient
\((\log p)p^{-m/2}\).  The endpoint is inclusive and has no half weight.
If \(p^m=\lambda^2\), the support overlap is a null endpoint.

The form domain is, with equivalence of graph norms,

\[
\mathcal Q_\lambda
=
\left\{
f\in\mathcal H_\lambda:
\int_{\mathbb R}
|\widehat f(t)|^2(1+\log(1+t^2))\,dt<\infty
\right\}.
\tag{3.5}
\]

The pole and finite prime-power sectors are bounded on
\(\mathcal H_\lambda\).  The form is densely defined, closed, lower
semibounded, inversion invariant, and has an associated selfadjoint
compact-resolvent operator.

## 4. Finite active and reference operators

The Fourier basis and free spectral scale are

\[
V_j(u)
=L_\lambda^{-1/2}
\exp\left(
\frac{2\pi ij\log(\lambda u)}{L_\lambda}
\right),
\qquad
d_{\lambda,j}=\frac{2\pi j}{L_\lambda}.
\tag{4.1}
\]

Put

\[
E_N=\operatorname{span}\{V_j:|j|\le N\}.
\tag{4.2}
\]

At an admissible pair, let \(\xi_{\lambda,N}\) be the inversion-even
lowest Ritz line.  The finite active selfadjoint operator
\(A_{\lambda,N}\) is the audited ZD-P quotient operator on the
\(2N\)-dimensional active quotient.

Let

\[
F_N=\operatorname{span}\{V_j:0<|j|\le N\}.
\tag{4.3}
\]

After the frozen reducing, inversion-invariant polar-transport axiom,
let

\[
B_{\lambda,N}
=U_{\lambda,N}(D_\lambda|_{F_N})U_{\lambda,N}^*.
\tag{4.4}
\]

This is project-canonical only after that axiom; it is not claimed to be
a comparison operator literally introduced in the primary source.  Its
spectrum is

\[
\sigma(B_{\lambda,N})
=\{d_{\lambda,j}:0<|j|\le N\}.
\tag{4.5}
\]

The finite projective relative determinant is

\[
\Delta_{\lambda,N}^{\rm rel}(z)
=\det[(A_{\lambda,N}-z)(B_{\lambda,N}-z)^{-1}],
\tag{4.6}
\]

and, in coefficient coordinates,

\[
\Delta_{\lambda,N}^{\rm rel}(z)
=-\frac{zL_\lambda^{-1/2}}
{\delta_N(\xi_{\lambda,N})}
\sum_{|j|\le N}
\frac{\xi_{\lambda,N,j}}{d_{\lambda,j}-z}.
\tag{4.7}
\]

No source-normalization scalar will be used in a limiting assertion.

## 5. Finite counting shift and Stieltjes transform

For a finite selfadjoint operator \(T\), define the right-continuous
count

\[
N_T(t)=\operatorname{Tr}\mathbf1_{(-\infty,t]}(T).
\tag{5.1}
\]

Multiplicities, negative eigenvalues, and zero eigenvalues are counted
literally.  Freeze

\[
\mathfrak X_{\lambda,N}(t)
=N_{B_{\lambda,N}}(t)-N_{A_{\lambda,N}}(t).
\tag{5.2}
\]

The two finite dimensions agree, so this shift vanishes below and above
both spectra.  Its value at a jump is immaterial in Lebesgue integrals.

For \(\Re s>0\), put

\[
g_{s,Y}(t)
=\frac{w_h(t)e^{-2(s-it)Y}}{s-it},
\tag{5.3}
\]

\[
\mathcal S_{\lambda,N,Y}(s)
=\operatorname{Tr}\left(
g_{s,Y}(A_{\lambda,N})
-g_{s,Y}(B_{\lambda,N})
\right).
\tag{5.4}
\]

The frozen relative sign is active minus reference.  The counting-shift
sign is reference minus active.

## 6. Infinite fixed-\(\lambda\) objects

At fixed \(\lambda\), assume explicitly that the infinite semilocal
ground eigenvalue is simple and isolated and has an inversion-even
normalized ground eigenvector \(\phi_\lambda\).  The existence of a
cofinal family \(\Lambda\) with this property is also an explicit
hypothesis, not a theorem of this phase.

The active divisor is the real analytic divisor

\[
\mathcal Z_\lambda
=Z(\widehat{\phi_\lambda})
\subset\mathbb R,
\tag{6.1}
\]

with multiplicities \(m_\lambda(\tau)\).  The active and free diagonal
smoothed kernels are

\[
Q_\lambda^{(h)}(y)
=
\sum_{\tau\in\mathcal Z_\lambda}^{\rm dist}
m_\lambda(\tau)w_h(\tau)e^{2iy\tau},
\tag{6.2}
\]

\[
Q_{\lambda,\rm ref}^{(h)}(y)
=
\sum_{j\in\mathbb Z\setminus\{0\}}
w_h(d_{\lambda,j})e^{2iyd_{\lambda,j}}.
\tag{6.3}
\]

Freeze

\[
\widetilde Q_\lambda^{(h)}
=Q_\lambda^{(h)}-Q_{\lambda,\rm ref}^{(h)},
\tag{6.4}
\]

\[
\mathcal S_{\lambda,Y}(s)
=2\int_Y^\infty
e^{-2sy}\widetilde Q_\lambda^{(h)}(y)\,dy.
\tag{6.5}
\]

The two infinite diagonal smoothed traces converge separately.  They are
not asserted to form a same-space resolvent-comparable Krein pair.

Let

\[
\delta_\lambda
=\mu_\lambda^{\rm ref}-\mu_\lambda^{\rm act}.
\tag{6.6}
\]

The centered infinite counting primitive is, for \(t>0\),

\[
\mathfrak X_\lambda(t)
=\delta_\lambda((0,t])
+\frac12\delta_\lambda(\{0\}),
\tag{6.7}
\]

and, for \(t<0\),

\[
\mathfrak X_\lambda(t)
=-\delta_\lambda((t,0))
-\frac12\delta_\lambda(\{0\}).
\tag{6.8}
\]

It satisfies

\[
\mathfrak X_\lambda'=\delta_\lambda
\quad\text{in }\mathcal D'(\mathbb R).
\tag{6.9}
\]

The notation \(\xi_\lambda^{\rm Krein}\) is reserved for a genuine
Krein spectral-shift function; no such function is asserted here.
The Riemann function, if mentioned, is \(\Xi_{\rm Riemann}\), never
\(\mathfrak X_\lambda\).

## 7. Arithmetic zero series and target

For distinct nontrivial zeta zeros \(\rho\), with multiplicity
\(m_\rho\), put

\[
w_\rho=\rho-\frac12.
\tag{7.1}
\]

The canonical arithmetic interaction is

\[
Q_h(y)
=\sum_\rho^{\rm dist}
m_\rho L_h(w_\rho)e^{2w_\rho y}.
\tag{7.2}
\]

The sum is ordinary absolute summation.  It and every fixed real
derivative converge absolutely and locally uniformly on
\(\mathbb R\).  Functional-equation and conjugation symmetries give its
canonical real, even extension.

The factor-two arithmetic Stieltjes target, including the deleted free
zero mode, is

\[
\boxed{
\begin{aligned}
\mathcal S_{h,Y}(s)
={}&
\sum_\rho^{\rm dist}
\frac{
m_\rho L_h(\rho-\frac12)
e^{-2(s-\rho+1/2)Y}
}{
s-(\rho-\frac12)
}\\
&+\frac{M_h^2e^{-2sY}}s.
\end{aligned}}
\tag{7.3}
\]

The deleted-zero correction has a plus sign in the relative transform.

## 8. Arithmetic truncation, free alias, and neutral defect

For \(n\ge2\), define

\[
\Phi_{Y,n}(s)
=
\int_{\max\{-1,\,2Y-\log n\}}^1
K_h(r)e^{-sr}\,dr,
\tag{8.1}
\]

with value zero if the lower endpoint is at least \(1\).  Put

\[
P_h=\int_{\mathbb R}h(x)e^{x/2}\,dx,
\qquad
P_h^2=L_h(1/2).
\tag{8.2}
\]

The arithmetic truncation is defined without active data:

\[
\begin{aligned}
\mathcal S_{\lambda,h,Y}^{\rm arith}(s)
={}&
\frac1{2\pi}\int_{\mathbb R}
2\theta'(t)L_h(it)
\frac{e^{-2(s-it)Y}}{s-it}\,dt\\
&+P_h^2\left[
\frac{e^{-2(s-1/2)Y}}{s-1/2}
+\frac{e^{-2(s+1/2)Y}}{s+1/2}
\right]\\
&-\sum_{\substack{n=p^m\le\lambda^2}}
\frac{\log p}{n^{s+1/2}}\Phi_{Y,n}(s)
+\frac{M_h^2e^{-2sY}}s.
\end{aligned}
\tag{8.3}
\]

Every retained prime power is included exactly once.

For \(L_\lambda\ge2Y+1\), freeze the complete free alias

\[
\mathcal A_{\lambda,Y}(s)
=\frac{L_\lambda L_h(s)}{e^{sL_\lambda}-1}.
\tag{8.4}
\]

Let

\[
Z_{\lambda,Y}^{\rm act}(s)
=\sum_{\tau\in\mathcal Z_\lambda}^{\rm dist}
m_\lambda(\tau)
\frac{w_h(\tau)e^{-2(s-i\tau)Y}}{s-i\tau},
\tag{8.5}
\]

and let

\[
Z_{\lambda,h,Y}^{\rm arith}
=\mathcal S_{\lambda,h,Y}^{\rm arith}
-\frac{M_h^2e^{-2sY}}s.
\tag{8.6}
\]

The neutral spectral--arithmetic defect is

\[
\boxed{
\mathcal D_{\lambda,Y}^{\rm spec-arith}(s)
:=
Z_{\lambda,Y}^{\rm act}(s)
-Z_{\lambda,h,Y}^{\rm arith}(s).
}
\tag{8.7}
\]

No smallness is asserted.  The exact source-level ledger is

\[
\boxed{
\mathcal S_{\lambda,Y}
=
\mathcal S_{\lambda,h,Y}^{\rm arith}
+\mathcal D_{\lambda,Y}^{\rm spec-arith}
-\mathcal A_{\lambda,Y}.
}
\tag{8.8}
\]

When \(L_\lambda\ge2Y+1\), the omitted-prime tail is

\[
R_{\lambda,h,Y}^{\rm prime}(s)
=-L_h(s)
\sum_{n>\lambda^2}
\frac{\Lambda(n)}{n^{s+1/2}}.
\tag{8.9}
\]

Thus safe-axis identification is exactly the unproved statement

\[
\sup_{s\in I_{\eta,S}}
\left|
\mathcal D_{\lambda,Y}^{\rm spec-arith}(s)
\right|\longrightarrow0,
\tag{8.10}
\]

because the free alias and omitted-prime tail are separately explicit
and tend to zero there.

## 9. Immutable logical exclusions

No branch may:

- change \(Y,\eta,S\), the packet, or any transform factor;
- add a subtraction depending on \(Q_h\), \(\Xi_{\rm Riemann}\), or zeta
  zeros;
- put the zeta divisor into an active spectrum;
- restore the deleted reference zero mode;
- call \(\mathcal D_{\lambda,Y}^{\rm spec-arith}\) a small remainder;
- assume determinant, ground-vector, or ground-divisor convergence;
- infer cross-\(\lambda\) control from a fixed-\(\lambda\) theorem;
- call \(\mathfrak X_\lambda\) a Krein shift;
- use RH, Weil positivity, or zero completeness;
- reopen a proxy, bottom-cluster, ground-state compactness,
  reflected-packet binding, Robin, colossally abundant, or
  moment-certificate route.
