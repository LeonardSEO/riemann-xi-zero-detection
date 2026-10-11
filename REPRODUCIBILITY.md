# Reproducibility

## Frozen version and integrity

For exact comparisons, check out tag `v1.0.0` or download its complete source archive. `SHA256SUMS` at the repository root covers the publication files; `certificates/SOURCE_HASHES.json` identifies the frozen computational inputs. Release-asset checksums are distributed separately with the release.

```sh
shasum -a 256 -c SHA256SUMS
python3 scripts/build_release_artifacts.py
cd dist
shasum -a 256 -c SHA256SUMS
cd ..
```

The release builder uses only Python's standard library, sets fixed archive timestamps, and checks every archived name and byte sequence against its source. This establishes package integrity and reproducibility of source archives. A freshly compiled PDF can differ in creation metadata and is not promised to have the archived PDF's byte hash.

## Compile the paper

`paper/paper.tex` is standalone. With an existing Tectonic installation:

```sh
cd paper
tectonic paper.tex
```

The archived PDF was compiled successfully with the built-in document compiler and exported with Tectonic 0.17.0. It contains 208 pages. The publication pass checked author metadata, labels, citations, unresolved references and physical clipping. Small overfull/underfull typography warnings remain.

## Certificate environment

The rational core programs use only Python's standard library. The Arb/Acb programs use `python-flint==0.9.0`; their exact arithmetic precision and geometry are specified in the programs and frozen outputs. To prepare a separate environment for replay:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
```

No dependency installation is needed for reading the paper, verifying hashes or rebuilding the source archives. The certificate index lists the commands and exact scopes. Run the short scripts first; the larger rational core and contour programs have higher computation costs.

Examples that write new output at an explicitly selected path:

```sh
python3 certificates/n1_tail_bound.py --certificate /tmp/n1-tail-replay.json
python3 certificates/n1_contour_attempt.py --certificate /tmp/n1-contour-replay.json
python3 certificates/head3_local_rouche.py --certificate /tmp/head3-replay.json
python3 certificates/heated_head3_local_rouche.py --certificate /tmp/heated-head3-replay.json
```

Some retained scripts write a sibling JSON file automatically. Replay those in a disposable copy of the `certificates/` directory, preserving its subdirectories and sibling imports. Keep the frozen release outputs when making comparisons.

## Limits of a replay

Successful replay certifies the stated finite inequalities, contour count or approximant node. The analytic tails and domain arguments are proved in the manuscript. A replay does not supply a global RH estimate.

The upstream finite-interlacing outputs include original Python runtime metadata. Earlier replays reproduced their mathematical intervals and counts but an exact metadata comparison differed between Python 3.12.11 and 3.12.14. That discrepancy is retained and is not described as a byte-identical metadata check.

Publication preparation preserved previously checked numerical outputs rather than rerunning the entire certificate collection. The new general analytic results do not depend on a numerical zero search or finite-height verification.
