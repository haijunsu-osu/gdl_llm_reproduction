# RR/PR Dyad Synthesis Report (Ch.5)

## 1) Commands run

```powershell
Get-Content Ch5_RR_PR_Prompt.md -TotalCount 300
rg -n "\bRR\b|\bPR\b|dyad|synthesis|Equation|Eq\.|equation|Freudenstein|precision|displacement" Ch5_GDL.txt
python rr_pr_dyad_synthesis.py
```

## 2) Chapter equations/sections used from `Ch5_GDL.txt`

- Section **5.2.1**: RR geometric constraint and perpendicular-bisector relation.
  - Eq. **(5.9)**: \((W_i-W_1)\cdot(G-P_{1i})=0\)
- Section **5.3.1**: RR algebraic design equations.
  - Eq. **(5.15)** and Eq. **(5.16)**: \(W_i-W_1=[A(\phi_{1i})-I](W_1-P_{1i})\)
  - Eq. **(5.17)** (primary RR synthesis equation used in solver):
    \[
    (G-P_{1i})\cdot[A(\phi_{1i})-I](W_1-P_{1i})=0,\ i=2,\dots,n
    \]
- Section **5.4.1**: PR design equations.
  - Eq. **(5.63)** and Eq. **(5.64)**: collinearity constraints for \(W_i\) on a line.
  - Eq. **(5.65)**: reuse of RR relative displacement relation.
  - Eq. **(5.66)** (primary PR synthesis equation used in solver):
    \[
    |S,\ [A(\phi_{1i})-I](W_1-P_{1i})|=0,\ i=2,\dots,n
    \]
  - Eq. **(5.67)–(5.70)**: bilinear form with \(S=(1,m)^T\), unknowns \(m,\lambda,\mu\).

## 3) Implementation deliverable

- Script: `rr_pr_dyad_synthesis.py`
- Includes:
  - RR equation setup and 5-position finite-solution solver.
  - PR equation setup and 4-position finite-solution solver.
  - Branch/solution-set handling via root enumeration + dedup.
  - Degeneracy checks (near-translation pole singularity checks).
  - Jacobian-rank checks, residual checks, and physical validity checks.
  - Angle normalization (`[-pi, pi]`) and inverse-trig domain clipping (`safe_acos`).

## 4) Input design data used for validation

No explicit chapter numeric table was available in this file extract, so consistent synthetic validation cases were generated.

### RR test (5 positions)

| i | \(\phi_i\) [rad] | \(d_i=(d_{xi},d_{yi})\) |
|---|---:|---|
| 1 | 1.279938950 | (1.288544523, -3.039043609) |
| 2 | -1.747753990 | (2.464658508, -2.256931326) |
| 3 | -2.491178890 | (1.305719413, 1.366112897) |
| 4 | -1.260171540 | (-0.012581529, -1.813866051) |
| 5 | 0.881624880 | (2.492716361, 0.298552581) |

### PR test (4 positions)

| i | \(\phi_i\) [rad] | \(d_i=(d_{xi},d_{yi})\) |
|---|---:|---|
| 1 | 0.500000000 | (-1.512691961, -1.976868401) |
| 2 | -1.000000000 | (0.307550937, -0.215863600) |
| 3 | 1.400000000 | (1.867671484, -0.773925269) |
| 4 | -0.300000000 | (-0.733911613, -0.929764415) |

## 5) Key computed synthesis results

### RR dyad results (all physically valid real solutions found)

1. \(G=(1.200000000,\ -0.700000000)\), \(W_1=(1.690653507,\ -2.741876376)\), branch `++-+`
2. \(G=(2.146853338,\ -0.975178386)\), \(W_1=(0.958479706,\ -2.489652218)\), branch `+---`
3. \(G=(2.154370157,\ -3.121502551)\), \(W_1=(-2.566974446,\ -4.945028209)\), branch `----`
4. \(G=(3.981729532,\ -0.801550036)\), \(W_1=(0.699991946,\ -1.190718025)\), branch `+-+-`

### PR dyad results (all physically valid real solutions found)

1. \(m=0.600000000\), \(W_1=(-0.714740803,\ -1.768844482)\)

## 6) Residual/error metrics

### RR

- For each of 4 solutions:
  - `rank(J) = 4` (full rank for 4 unknowns).
  - `max |Eq(5.17)| residual` in \([1.776e-15,\ 3.109e-15]\).
  - Radius consistency (`max(r_i)-min(r_i)`) in \([1.554e-15,\ 3.553e-15]\).

### PR

- For the 1 solution:
  - `rank(J) = 3` (full rank for 3 unknowns).
  - `max |Eq(5.66)| residual = 2.220e-16`.
  - Max line-collinearity residual = `2.776e-16`.

## 7) Console output status

`python rr_pr_dyad_synthesis.py` printed:
- Input poses for RR and PR tests,
- relative poles and relative rotation angles,
- all candidate solutions with branch signatures / parameters,
- residual metrics and validity flags,
- final summary:
  - `RR pipeline: PASS`
  - `PR pipeline: PASS`

## 8) Final validation status

- RR synthesis workflow: **PASS**
- PR synthesis workflow: **PASS**
- Overall task completion: **PASS**
