# TVL.py — Topological Vortex Logic

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19683376.svg)](https://doi.org/10.5281/zenodo.19683376)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

`TVL.py` is a single-file Python implementation of the classification map studied in the
associated paper: it assigns to every integer winding vector w ∈ ℤ³ on the three-torus its
stability, shell, selected mod-three charge class, traceless-projection data, and weight
class, and it verifies the root-system and module-theoretic structures the paper proves.
It is a verification and exploration tool, not a proof: every result it exercises is
established analytically in the paper, and no proof depends on the code.

The module is deliberately split into two layers, so that the mathematics can be used
without the physical reading.

**Associated paper** (a separate Zenodo record; the software is deposited on its own). Vladimer Merebashvili, *Topological Vortex Logic: Stability, Root
Systems, and Module Structure of Winding States on T³ with a Selected Z₃ Grading*,
[10.5281/zenodo.19682633](https://doi.org/10.5281/zenodo.19682633) (concept DOI, always
resolves to the latest version).

**Requirements.** Python 3.8 or later. Standard library only — no third-party
dependencies, and no installation step: download `TVL.py` and run it.

## Mathematical core

`TVL.classify(w)` returns an immutable `TVLState` carrying only geometric quantities:
the winding vector, |w|², stability, shell (face/edge/corner), shell index, charge class
q₃ = tr(w) mod 3, the traceless-projection norm |w_t|², and the weight class. Supporting
methods cover the closed-form stability argument with an explicit favourable-split
witness for unstable inputs (`closed_form_stable`), enumeration of the 26 stable states
(`all_stable`), the coordinate-C₃ permutation-module invariants
(`z3_module_invariants`), and verifiers for both root systems
(`verify_a2_root_system`, `verify_b3_root_system`).

The winding class is typed as an element of H¹(T³; ℤ). Two distinct order-three
structures appear and are kept apart throughout: **Z₃** is the charge grading, the target
of the diagonal homomorphism w ↦ tr(w) mod 3; **C₃** is the coordinate cycle
g: (w₁,w₂,w₃) ↦ (w₂,w₃,w₁) acting on the shells. Cartan matrices follow the paper's
row-coroot convention.

## Optional conjectural adapter

`TVLInterpretation.read(w)` returns a `PhysicalReading` supplying a baryon number and a
particle-sector label. These are an interpretation matched to the mathematics, not
consequences of it; they are kept in a separate object, are never attached to the
mathematical state, and are refused for vacuum and unstable vectors.

## Results exercised by the built-in suite (65 checks)

The 26-state stable vocabulary and the 6/12/8 shells; the 8/9/9 charge-class
distribution; the weight classes, including the six-element orbit {±2μᵢ}, which is
reported as a weight orbit and not as any single irreducible representation, and its
split by charge class — the q₃ = 1 non-diagonal corners project into {−2μᵢ} and the
q₃ = 2 corners into {+2μᵢ}; the A₂ root system (root count, reducedness, integrality,
reflection closure, Cartan matrix); the B₃ root system, with its Weyl group generated
from the simple-root reflections, checked to have order 48 and to coincide with the group
of signed coordinate permutations; and the coordinate-C₃ module invariants — character
traces (0, 0, 2) and invariant-subspace dimensions (2, 4, 4) — that distinguish the three
shells.

## Command line

```bash
python TVL.py --test              # self-test + both root-system verifiers
python TVL.py --all               # list the 26 stable states
python TVL.py --map               # tabulate every state with |w|^2 <= 6
python TVL.py --all --interpret   # add the conjectural reading to a listing
python TVL.py 1,0,0 1,1,-1        # classify individual winding vectors
```

`--interpret` is a modifier: combine it with `--all` or `--map`. Output is
mathematics-only unless `--interpret` is given.

## How to cite

Cite the software by its concept DOI, which always resolves to the latest version:

> Merebashvili, V. *TVL.py — Topological Vortex Logic*. Zenodo.
> https://doi.org/10.5281/zenodo.19683376

If you use the results themselves, cite the associated paper
([10.5281/zenodo.19682633](https://doi.org/10.5281/zenodo.19682633)) as well: the proofs
are there, not in the code.

## Changelog

#### v1.0.8 — 3 August 2026

Synchronised with the associated paper's v1.0.8 release
([10.5281/zenodo.21765104](https://doi.org/10.5281/zenodo.21765104)). No numerical value
changes.

- The corner-orbit charge split is added as classifier output and verified in the suite:
  the q₃ = 1 non-diagonal corners project into {−2μᵢ} and the q₃ = 2 corners into
  {+2μᵢ}, the direction established in the paper's corrected statement.
- Group naming aligned with the paper throughout: Z₃ is the charge grading, C₃ the
  coordinate cycle.
- The module-invariant table no longer names the associated paper by series label in its
  printed output; it is referred to as the associated paper.
- Documentation-layer clarifications.
- Self-test extended to 65 checks.

#### v1.0.7 — 22 July 2026

Corrects the v1.0.6 deposit, which was published without files. The changes below
were intended for v1.0.5 and first ship here.

- Companion paper no longer included; the software is deposited on its own.
- Split into a mathematical core (`TVL`, `TVLState`) and an optional adapter (`TVLInterpretation`, `PhysicalReading`); physical labels no longer appear on the state object.
- Weight labels made representation-neutral: the corner projections are the orbit {±2μᵢ}, not an irreducible representation; the sextet label is removed.
- `generation` renamed `shell_index`; "forbidden generation" and "gauge charge" wordings removed.
- Charge-class counts (8/9/9) and the adapter's sector partition stated separately.
- Per-state provenance labels removed; provenance is now the layer a quantity belongs to.
- Coefficient field corrected: three characters over ℂ or ℚ(ω); over ℚ, trivial plus a two-dimensional rational irreducible; neither over F₃.
- Every public vector helper validates its argument as exactly three integers.
- The adapter refuses vacuum and unstable vectors; `read()` returns `None`.
- Output is mathematics-only unless `--interpret` is given.
- Verification extended: a full A₂ verifier; the B₃ Weyl group generated from the simple reflections (order 48); reducedness checked over all parallel pairs.
- Cartan matrices use the row-coroot convention, matching the paper.
- The coordinate action is named C₃, distinct from the charge target Z₃; invariants report invariant-subspace dimension.
- The winding class is typed as an element of H¹(T³; ℤ).
- Self-test reports dynamically; now 62 checks.

#### v1.0.6 — 22 July 2026

Published without files; superseded by v1.0.7.

#### v1.0.5 — 22 July 2026

Published with the v1.0.4 files in error; superseded by v1.0.6.

#### v1.0.4 — 26 April 2026 · [10.5281/zenodo.19779072](https://doi.org/10.5281/zenodo.19779072)

Synchronised with the companion paper's v1.0.4 release. No change to the classification
logic.

#### v1.0.3 — 25 April 2026 · [10.5281/zenodo.19752428](https://doi.org/10.5281/zenodo.19752428)

The repository reference in the companion paper is updated. No change to the
classification logic.

#### v1.0.2 — 25 April 2026 · [10.5281/zenodo.19751573](https://doi.org/10.5281/zenodo.19751573)

Synchronised with the companion paper's v1.0.2 release. No change to the classification
logic.

#### v1.0.1 — 22 April 2026 · [10.5281/zenodo.19688646](https://doi.org/10.5281/zenodo.19688646)

Synchronised with the companion paper's v1.0.1 release. No change to the classification
logic.

#### v1.0.0 — 21 April 2026 · [10.5281/zenodo.19683377](https://doi.org/10.5281/zenodo.19683377)

Initial release: the 26-state classification map, B₃ root-system verification, and the
Z₃-module invariant table.

## License

Released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
