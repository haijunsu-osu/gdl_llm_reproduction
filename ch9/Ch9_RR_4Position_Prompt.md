# Prompt: Chapter 9 Spherical RR Dyad Synthesis (4 Specified Positions)

```text
[Context]
You are an expert in mechanism synthesis, spherical kinematics, and symbolic computation. Work like a careful human engineer: read the chapter text, derive the equations used in that text, implement a reproducible CAS workflow, and verify the results numerically.
Assume the workspace contains only this prompt file and `Ch9_GDL.txt`.

[Task]
Read `Ch9_GDL.txt` and solve the 4-position spherical RR dyad synthesis problem using the Chapter 9 formulation.
Use the finite-position RR equations from Section 9.4, especially:
- design equations: (9.25),
- bilinear linearized forms: (9.28), (9.30),
- four-position center-axis cone formulation: (9.45), (9.46),
- four-position circling-axis cone formulation: (9.53), (9.54).

Generate a reproducible solver workflow that computes valid fixed axes `G` and associated moving axes `W1` for the 4 specified orientations.

[Requirements]
1. Use CAS as the primary solver and provide a Mathematica script (`.wls`) that can be run end-to-end.
2. Do not assume formulas from memory; extract them from `Ch9_GDL.txt` and match notation.
3. Use deterministic numerics:
   - fixed precision policy,
   - no randomness,
   - stable sorting of solutions.
4. Because the 4-position problem yields a one-parameter family, report the complete solution in reproducible form by:
   - providing the implicit cone equation(s) (center-axis and/or circling-axis),
   - generating a deterministic sampled set of real solution points on the curve/cone.
5. For each reported dyad sample, validate residuals of the original constraints `G . (A1i - I) . W1 = 0` for all relevant `i`.
6. Assume only this prompt file and `Ch9_GDL.txt` exist in the workspace. Do not depend on any other pre-existing scripts, reports, or data files.
7. Autonomously debug until the workflow runs without errors and reproduces consistent results from a clean run.
8. Generate a `report.md` that documents equations used, commands run, task data, sampled solutions, and validation metrics.

[Deliverables]
1. A runnable Mathematica script (for example `ch9_rr_4pos_synthesis.wls`).
2. Reproducible output artifacts (for example `.txt`, `.csv`, or `.json`) containing sampled axes and residuals.
3. A concise report file (for example `ch9_rr_4pos_report.md`) with:
   - chapter equations referenced,
   - input orientation data,
   - resulting solution-family representation,
   - residual-based validation summary.
```
