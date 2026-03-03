# Chapter 9 RR Dyad Synthesis (4 Positions)

## Scope and Assumptions
I followed the Chapter 9 spherical RR finite-position formulation in `Ch9_GDL.txt`.
The provided chapter text excerpt does not include a fully numeric 4-orientation table for Section 9.4.5, so the script uses a deterministic synthetic 4-orientation task generated from a known spherical RR motion, then solves it using the Chapter 9 equations.

## Equation Traceability to `Ch9_GDL.txt`
- Design constraint (9.25): `G . (A1i - I) . W1 = 0`, lines 553-555.
- Bilinear linearized form in `W1` (9.28): lines 614-616, with coefficients in (9.29), lines 619-643.
- Bilinear linearized form in `G` (9.30): lines 649-651, with coefficients in (9.31), lines 653-680.
- Four-position center-axis matrix form (9.45): lines 907-913.
- Center-axis cone condition (9.46): lines 916-919.
- Four-position circling-axis matrix form (9.53): lines 1074-1081.
- Circling-axis cone condition (9.54): lines 1086-1088.

## Reproducible Workflow
Script:
- `ch9_rr_4pos_synthesis.wls`

Command run:
```powershell
wolframscript -file ch9_rr_4pos_synthesis.wls
```

Determinism policy in script:
- Fixed arithmetic precision: `prec = 70`, `solvePrec = 50`.
- No randomness.
- Deterministic sampling grids (`Subdivide`) and lexicographic sorting of solutions.
- Fixed axis-sign canonicalization for line-axis equivalence (`v` and `-v`).

## Input Orientation Data (Task Set)
Synthetic task generation parameters (stored in `ch9_rr_4pos_results.json`):
- True fixed axis `G_true = (0.05383896085310486, -0.03828548327331901, 0.9978154078108767)`.
- True moving axis `W1_true = (-0.6096647765340221, 0.7295988309341576, 0.30982964053368336)`.
- Joint angle schedules (degrees):
  - `beta = {0, 35, -25, 60}`
  - `alpha = {0, -20, 30, -40}`
- Orientations are generated as:
  - `Ai = Rot(G_true, beta_i) . Rot(W1_true, alpha_i)`, `i = 1..4`.

Relative rotations `A1i` axis-angle summary:
- `A12`: `phi = 0.6208859301715984`, `S = (0.5091805061500293, -0.3258184009978977, 0.7966037796345013)`
- `A13`: `phi = 0.5890209749951565`, `S = (-0.42671085013606156, 0.7833169595830453, -0.4520314051100138)`
- `A14`: `phi = 1.0855606467255544`, `S = (0.6435435870010463, -0.24628429753472163, 0.7247038680852031)`

## Solution-Family Representation
The script computes both cubic cones:
- Center-axis cone `R(x,y,z)=0` from `det(M)=0` (Eq. 9.45 -> 9.46).
- Circling-axis cone `R'(lam,mu,nu)=0` from `det(M0)=0` (Eq. 9.53 -> 9.54).

Full implicit equations are exported to:
- `ch9_rr_4pos_center_cone.txt`
- `ch9_rr_4pos_circling_cone.txt`

Center-axis cone coefficients (monomial order of Eq. 9.46 style):
- `y^3 = 0.0091764451825446636839392173`
- `x y^2 = 0.0356446520046096636654723656`
- `z y^2 = -0.0075537142060945607992550500`
- `x^2 y = 0.0062744191857913590663524068`
- `x y z = -0.0076107259690745891883770856`
- `z^2 y = 0.0213258663913679326416991443`
- `x^3 = 0.0040974513088247277386058456`
- `x^2 z = -0.0247852940370611194449833637`
- `x z^2 = 0.0239804285229741919790229698`
- `z^3 = -0.0004103825690695823960357541`

Circling-axis cone coefficients:
- `mu^3 = 0.0323743673461409520195803128`
- `lam mu^2 = 0.0685058380113524507705231939`
- `nu mu^2 = 0.0252146600099852091335492054`
- `lam^2 mu = 0.0529486255168307548483029792`
- `lam nu mu = -0.0141871107232220532575065644`
- `nu^2 mu = 0.0324112744238769070960762046`
- `lam^3 = 0.0252767547842253028458657527`
- `lam^2 nu = -0.0446852046651736053913674693`
- `lam nu^2 = 0.0364497260176964626974650897`
- `nu^3 = -0.0025752599682498106284440602`

## Deterministic Sampled Dyads and Validation
The script deterministically samples real center-axis directions on the center cone, computes associated `W1` from `NullSpace(M)`, and validates Eq. (9.25) for `i = 2,3,4`.

Summary:
- Sampled/validated dyads: `30`
- Residual magnitudes from `ch9_rr_4pos_samples.csv`:
  - `maxResidual` values: `0.e-49` or `0.e-50`
  - `centerConeResidual` values: `0.e-52` or `0.e-53`
  - `circlingConeResidual` values: `0.e-50` or `0.e-51`

First five sampled dyads (full list in CSV):

| sampleId | Gx | Gy | Gz | W1x | W1y | W1z | maxResidual |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 0.01936959571648936 | 0.9998123917824682 | 0.6190060432490997 | -0.7568521847054802 | -0.20977675974148163 | 0.e-50 |
| 2 | 0.0014209953410267843 | 0.017762441762834803 | 0.9998412256127783 | 0.618904217650725 | -0.7562434829865174 | -0.21225777680585696 | 0.e-49 |
| 3 | 0.001704303578243604 | -0.02130379472804505 | -0.9997715957554998 | 0.6191189455933349 | -0.7575756124688896 | -0.20681132125631832 | 0.e-49 |
| 4 | 0.002624777358044533 | 0.016404858487778332 | 0.9998619860569841 | 0.6188124408940215 | -0.7557238400815283 | -0.2143656701228345 | 0.e-49 |
| 5 | 0.003658150608859869 | 0.015242294203582788 | 0.999877137653189 | 0.618729592314305 | -0.7552748080202377 | -0.2161796844399748 | 0.e-49 |

## Deliverables Produced
- `ch9_rr_4pos_synthesis.wls`
- `ch9_rr_4pos_results.json`
- `ch9_rr_4pos_samples.csv`
- `ch9_rr_4pos_center_cone.txt`
- `ch9_rr_4pos_circling_cone.txt`
- `ch9_rr_4pos_report.md`
