# CliffordLab – Action List

This checklist prioritizes making the repo production‑grade **and** directly useful to your articles on **D‑holomorphy** and **Clifford slice regularity**.

## Repo Hygiene & Structure
- [ ] Strengthen `README.md`: add quickstart and cross‑links to the two articles (D‑holomorphy, slice‑regularity).
- [ ] Add `__init__.py` to package dirs (so Python treats them as packages): examples, tests, src.
- [ ] Create a clean top‑level layout:
  - `src/cliffordlab/` (library)
  - `examples/` (small runnable scripts)
  - `notebooks/` (exploratory)
  - `papers/` (LaTeX, figures)
  - `scripts/` (CLI utilities)
  - `tests/`

- [ ] Add `pyproject.toml` (PEP 621) with `project` metadata and dependencies.
- [ ] Add `requirements.txt` with pinned minimal versions; provide `pip install -r requirements.txt`.
- [ ] Provide optional extras: `pip install cliffordlab[papers]` to include `sympy`, `numpy`, `scipy`, `matplotlib`, `numba`, `jax` (if used).
- [ ] Add `ruff` + `pytest` + `mypy` config; wire to CI (GitHub Actions).
## Dependency Audit (from imports)
- [ ] Ensure these appear in `requirements.txt`: cliffordlab, sympy, numpy.
## Core Library – Features supporting your papers
- [ ] **D‑algebra primitives**: dual‑unit `ε` with `ε²=0`, idempotents `e±`, and algebra ops.
- [ ] **D‑holomorphic module**: 
  - `dalembert_factorization(f) -> (F,G)` such that `f(u,v)=F(u)e+ + G(v)e-`
  - Cauchy‑type reconstruction on two‑annulus domains; double Laurent expansion
  - **Characteristic residue / jump law** primitives and tests
- [ ] **Interfaces & energy law**: 2×2 reflection/transmission map; numeric example & unit tests.
- [ ] **Slice‑regular (Clifford)**: 
  - generators for slice functions; stem function construction
  - Cauchy kernel & power series; Liouville‑type results where applicable
  - validation tests vs known quaternionic/Clifford identities
- [ ] **Symbolics** (`sympy`) for exact identities; **numerics** (`numpy`) for experiments.
- [ ] **Figure reproducibility**: `scripts/make_figures.sh` to reproduce all plots used in articles.
## Notebooks & Reproducibility
- [ ] Move existing notebooks under `notebooks/` and add `nbconvert` to export Python for CI runs.
- [ ] Provide a `Makefile` target: `make reproduce` to run all experiments and stash results under `results/`.
- [ ] Save deterministic seeds; persist figures as SVG/PGF for LaTeX.
## Testing
- [ ] `tests/test_dalembert.py`: numeric check that `∂_v F(u)=0`, `∂_u G(v)=0` under decomposition.
- [ ] `tests/test_jump.py`: verify discrete jump equals line residue across a characteristic.
- [ ] `tests/test_slice.py`: slice Cauchy formula and basic equivalences.
- [ ] Target ≥90% coverage on core algebra and identities.
## CLI & Examples
- [ ] `cliffordlab d-analytic --example pulse` to generate the interface reflection/transmission plot.
- [ ] `cliffordlab slice --example cauchy` to compute a series expansion and export coefficients.
- [ ] Short `examples/*.py` mirroring the CLI.
## Docs
- [ ] Sphinx or MkDocs site with API pages and “From the papers” gallery linking to figures.
- [ ] Add cross‑refs and badges (arXiv, DOI, CI, coverage).