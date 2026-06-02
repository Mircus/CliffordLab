# Clifford Twins

Experimental Python package accompanying the note:

**Clifford Algebras and Their Nonassociative Twins: Passive Volume, Active Volume, and Folded Triality**

Core slogan:

> \(\mathrm{Cl}(8,0)\) is unfolded triality.  
> \(\mathbb O\) is folded triality.  
> \(\mathrm{Cl}(3,0)_{\mathrm{twin}}\) is the active-volume route from 3D rotor geometry to the octonionic triality cell.

## What this repo does

This repo implements a small computational laboratory for the article's central construction:

\[
\mathrm{Cl}(3,0)=H\oplus HI,
\qquad
H=\mathrm{Cl}^+(3,0)\cong\mathbb H.
\]

On the same 8-dimensional real vector space \(H\oplus HI\), we compare two products.

### 1. Clifford product: passive volume

The Clifford pseudoscalar \(I=e_1e_2e_3\) is central:

\[
Ia=aI.
\]

Therefore

\[
(a+bI)(c+dI)=(ac-bd)+(ad+bc)I.
\]

This product is associative.

### 2. Octonion twin product: active volume

Replace the passive rule by the conjugating rule:

\[
\ell a=\bar a\ell.
\]

Then

\[
(a+b\ell)\star(c+d\ell)
=
(ac-\bar d\,b)+(da+b\bar c)\ell.
\]

This is the Cayley-Dickson octonion product written on the intrinsic \(\mathrm{Cl}(3,0)\) split.

## What this repo does not claim

This is **not** a new construction of the octonions.

Known prior structures include:

- Cayley-Dickson doubling;
- \(\mathrm{Cl}^+(3,0)\cong \mathbb H\);
- octonions and Clifford algebras as twisted/quasi-twisted \((\mathbb Z_2)^n\)-graded algebras;
- octonionic Clifford representations and Spin(8) triality.

The point here is computational and conceptual:

> compare passive Clifford volume and active octonionic volume on the same intrinsic \(H\oplus HI\) Clifford body.

## Install

```bash
pip install -e ".[dev]"
```

## Run tests

```bash
pytest
```

## Run examples

```bash
python examples/cl3_twin_demo.py
python examples/folded_triality_demo.py
```

## Expected demo idea

The Clifford product is associative:

\[
[(x y) z] - [x (y z)] = 0.
\]

The octonion twin product is nonassociative but alternative. A basic witness is:

\[
(\ell i)j \neq \ell(ij).
\]

This witnesses the article's physical slogan:

> nonassociativity records where the active chirality/volume gate enters a sequence of rotations.

## Repo layout

```text
clifford-twins-repo/
├── docs/
│   ├── article/
│   └── notes/
├── examples/
├── notebooks/
├── src/
│   └── clifford_twins/
└── tests/
```

## CliffordLab integration

This repo can be used as:

```text
CliffordLab/
└── experimental/
    └── clifford-twins/
```

Later it can be folded into:

```text
CliffordLab/src/cliffordlab/twins/
```

after the interface stabilizes.
