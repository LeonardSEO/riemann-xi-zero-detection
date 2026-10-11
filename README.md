# Exact Zero Detection and Compact Arithmetic Realizations for the Riemann Xi Function

**Leonard van Hemert · Mathematical preprint · Version 1.0.0 · 11 October 2026**

[Read the paper](paper/paper.pdf) · [Version 1.0.0 release](https://github.com/LeonardSEO/riemann-xi-zero-detection/releases/tag/v1.0.0) · [LaTeX source](paper/paper.tex) · [Citation metadata](CITATION.cff)

This paper proves zero-detection and arithmetic-realization theorems for the Riemann xi function. It gives exact domain and multiplicity statements, quantitative costs for finite completion, and explicit lower bounds forced by an individual hypothetical off-line zero.

**RH status: unresolved.** The paper proves the stated intermediate theorems and conditional implications. The global arithmetic upper bound needed to exclude every off-line zero remains open. No simplicity of the zeta zeros is assumed.

## Principal results

| Result | Precise scope |
|---|---|
| **Exact negative index** | The exhausted negative index of the complete semilocal Weil form is $2K_+$, where $K_+$ counts distinct off-line mirror pairs at positive height. The statement includes infinitely many pairs and arbitrary zero orders. |
| **Local index with moment constraints** | The genuinely logistic-smoothed form attains the same index on every nonempty bounded open interval, after any fixed finite family of real exponential moment and jet conditions. Increasing moment depth has an explicit cost. |
| **Compact arithmetic realization** | On the actual complete closed arithmetic range, a compact injective operator with dense range has the explicitly identified xi-zero spectrum. Every nonzero eigenvalue is semisimple with a two-dimensional eigenspace. Selfadjointness is not assumed or proved. |
| **Fractional complete arithmetic** | The full family retains the gamma term, both pole terms and every prime power. Its factorial Dirichlet–Hankel block is positive; its pool-null interaction retains the exact exceptional-zero index at every finite smoothing parameter. |
| **Fixed-band individual detection** | For a compact nonzero input, every nonempty fixed Fourier band retains every zero pole and its original multiplicity. A zero of displacement $\delta>0$ and order $m$ forces cumulative source energy on the limsup scale $Y^{2\delta}(\log Y)^{2m-2}$, with an explicit positive constant. |
| **Paid finite gamma transfer** | A single integer cap, independent of the input, simultaneously controls the unweighted and logarithmic Möbius fields on the full weighted space. At the stated cubic cap scale, the absolute energy error is $o(1)$ for the allowed growing input families. The exact primitive-sector remainder is also proved. |

The full statements, hypotheses and proofs are in [the manuscript](paper/paper.tex). The introduction identifies the principal theorems and the conclusions specify the remaining global obligations.

## Read and download

The [PDF](paper/paper.pdf) is a standalone 208-page manuscript. The tagged release provides:

- `riemann-xi-zero-detection-v1.0.0-paper.pdf` — compiled manuscript;
- `riemann-xi-zero-detection-v1.0.0-latex-source.tar.gz` — standalone LaTeX source;
- `riemann-xi-zero-detection-v1.0.0-certificates.tar.gz` — finite certificate programs, frozen outputs and reproduction instructions;
- `riemann-xi-zero-detection-v1.0.0-source.tar.gz` — complete publication source package;
- `SHA256SUMS` — checksums of the four release assets.

The GitHub release is submitted through the linked Zenodo integration for archival and a DOI. The verified DOI and its citation will be added here after archival completes.

## Repository contents

```text
paper/                   Standalone English manuscript and compiled PDF
certificates/            Used finite programs and frozen evidence
scripts/                 Deterministic release packaging
CITATION.cff             Machine-readable preferred preprint citation
.zenodo.json             Zenodo preprint metadata
REPRODUCIBILITY.md        Compilation, hashes and certificate replay
PROVENANCE.md            Source versions and verification scope
SHA256SUMS               Checksums of the publication files
NOTICE                   Rights and upstream license exception
```

The general analytic theorems require no numerical certificate. The numerical material concerns finite positivity cores, approximating theta heads and finite operator models. Each item has its own scope in [the certificate index](certificates/README.md). These computations do not certify an off-line zeta zero or settle RH.

## Reproduce this version

Use the `v1.0.0` tag when comparing against the archived release. Verify the publication files with:

```sh
shasum -a 256 -c SHA256SUMS
```

Compile the source using an existing LaTeX installation, for example Tectonic:

```sh
cd paper
tectonic paper.tex
```

Rebuild the deterministic source archives and their asset checksums with:

```sh
python3 scripts/build_release_artifacts.py
```

Certificate replay uses Python's standard library or `python-flint==0.9.0`, as specified in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). The release builder verifies the archived file names and bytes. A file hash or successful compilation is not a mathematical proof check.

## Verification and research scope

The paper contains complete proofs with explicit hypotheses. The principal established proof ranges have internal reconstruction records. The latest fixed-band and joint gamma arguments were reconstructed for their used scope during manuscript integration; a separate independent audit of the entire new range is not claimed.

Proof, originality and external review are distinct. A result-specific priority assessment has not been completed. External endorsement, journal acceptance and proof-assistant formalization are not claimed. The “Proof and reproducibility scope” paragraph in the conclusions gives the mathematical control limits in one place.

The remaining arithmetic target is an eventual subpolynomial bound on the cumulative energy of the fixed visible input, retaining the signed main term and the high primitive complement. Favorable finite compressions or almost-all modulation heights do not supply that bound.

## Related papers

Version-specific foundations and computational sources are cited in the manuscript:

- [Finite Positivity Certificates and a Cross-Endpoint Obstruction for a Riemann Xi-Defect Kernel — v0.2.1](https://github.com/LeonardSEO/certified-riemann-xi-positivity/releases/tag/v0.2.1)
- [Reflected-Packet Weil Kernels — v1.0.0](https://github.com/LeonardSEO/reflected-packet-growth-criteria/releases/tag/v1.0.0)
- [Reflected-Packet Interactions in the Semilocal Weil Form — v1.0.0](https://github.com/LeonardSEO/semilocal-reflected-packet-oscillation/releases/tag/v1.0.0)
- [Smoothed Positive Kernels for Finite Zeta Spectral Triples — v1.0.0](https://github.com/LeonardSEO/smoothed-zeta-spectral-kernels/releases/tag/v1.0.0)
- [Fixed-Parameter Spectral Approximation and Smoothed Relative Traces — v1.0.0](https://github.com/LeonardSEO/semilocal-zeta-operator-theory/releases/tag/v1.0.0)

## Citation and rights

Use [CITATION.cff](CITATION.cff), which selects the mathematical preprint as the preferred citation. Cite the version-specific DOI once Zenodo archival is confirmed.

Copyright © 2026 Leonard van Hemert. Rights are described in [NOTICE](NOTICE). The upstream Xi-moment program retains its accompanying AGPL-3.0 license; public availability does not grant a new license for the manuscript.
