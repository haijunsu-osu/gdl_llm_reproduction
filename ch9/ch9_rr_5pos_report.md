# Chapter 9 RR Dyad Synthesis (5 Positions)

## Scope
This run implements the Section 9.4.6 algebraic elimination workflow for 5-position spherical RR synthesis using `Ch9_GDL.txt` equations:
- (9.25) design constraints, `G . (A1i - I) . W1 = 0`, lines 553-555.
- (9.55) 4x3 matrix form in `(lambda, mu, nu)`, lines 1141-1150.
- (9.56) four cubic minors `R_j(x,y,z)`, lines 1154-1161.
- (9.57) dehomogenized cubic form in `y`, lines 1169-1172.
- (9.58) coefficient matrix `D` in powers `(y^3,y^2,y,1)`, lines 1182-1213.
- (9.59) sixth-degree elimination polynomial `P(x)=det(D)`, lines 1224-1227.

## Reproducible Run
Script:
- `ch9_rr_5pos_synthesis.wls`

Command:
```powershell
wolframscript -file ch9_rr_5pos_synthesis.wls
```

Determinism:
- Fixed precision (`prec=80`, `solvePrec=60`).
- No randomness.
- Deterministic sorting of roots by `(Re, Im)`.
- Deterministic solution ordering by `(x, y)`.

## Input Orientation Data
The chapter excerpt in `Ch9_GDL.txt` does not provide a complete numeric 5-orientation table in this workspace, so the script uses a deterministic synthetic 5-orientation task generated from a known spherical RR chain (fully recorded in JSON output):
- `G_true = (0.05591009679304701719, -0.04056222708515175757, 0.99761153101319187528)`
- `W1_true = (-0.72048049858468329009, 0.63165413574547576118, 0.28621828025966870428)`
- `beta_deg = {0, 28, -33, 47, -59}`
- `alpha_deg = {0, -18, 24, -31, 42}`
- `Ai = Rot(G_true,beta_i) . Rot(W1_true,alpha_i)`, `i=1..5`

## Elimination Output
- `deg(P)=6` (as required by Eq. 9.59).
- Polynomial coefficients (descending powers):
  - `x^6 = 1.42626850821444e-9`
  - `x^5 = 3.36747989621596e-9`
  - `x^4 = 1.58557888851542e-9`
  - `x^3 = -6.42840811023694e-10`
  - `x^2 = -5.61531198082361e-11`
  - `x^1 = 5.16696527840667e-13`
  - `x^0 = 2.43025751881786e-13`

All roots of `P(x)` (sorted):
1. `-1.2940836617896266 - 0.2188362540415963 i`
2. `-1.2940836617896266 + 0.2188362540415963 i`
3. `-0.06423858321342326 - 0.04201801609255440 i`
4. `-0.06423858321342326 + 0.04201801609255440 i`
5. `0.056043956043956044`
6. `0.29955850238021263`

Multiplicity summary (grouped numerically): all roots are simple (`multiplicity = 1`).

## Burmester Axes and Associated Moving Axes
Real roots produced 2 validated Burmester solutions:

1. From `x=0.05604395604395604`, `y=-0.04065934065934066`
- `G = (0.0559100967930470, -0.0405622270851518, 0.9976115310131919)`
- `W1 = (0.7204804985846833, -0.6316541357454758, -0.2862182802596687)`
- `max |G.(A1i-I).W1|, i=2..5 = 0.e-60`

2. From `x=0.29955850238021263`, `y=-0.93539813138369824`
- `G = (0.2137140034235429, -0.6673410297637949, 0.7134299368084295)`
- `W1 = (0.7182425079435887, -0.14705582516969966, 0.6800752047138336)`
- `max |G.(A1i-I).W1|, i=2..5 = 0.e-58`

Note: axis direction sign is physically equivalent (`v` and `-v` are the same line axis), so the first `W1` is sign-opposite to the generator `W1_true` and still valid.

## Deliverables
- `ch9_rr_5pos_synthesis.wls`
- `ch9_rr_5pos_results.json`
- `ch9_rr_5pos_roots.csv`
- `ch9_rr_5pos_root_groups.csv`
- `ch9_rr_5pos_solutions.csv`
- `ch9_rr_5pos_Px.txt`
- `ch9_rr_5pos_report.md`
