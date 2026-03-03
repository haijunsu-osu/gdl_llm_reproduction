# Prompt: Chapter 9 Spherical RR Dyad Synthesis (5 Specified Positions)

```text
[Context]
You are an expert in mechanism synthesis, spherical kinematics, and symbolic/numeric algebraic elimination. Work like a careful human engineer: read the chapter text, derive the exact elimination workflow, implement reproducible CAS code, and verify all computed solutions.
Assume the workspace contains only this prompt file and `Ch9_GDL.txt`.

[Task]
Read `Ch9_GDL.txt` and solve the 5-position spherical RR dyad synthesis problem for all solutions.
Follow Section 9.4.6 and implement the algebraic elimination method using:
- design equations: (9.25),
- matrix setup for five positions: (9.55),
- cubic minors `R_j`: (9.56),
- dehomogenized polynomial form: (9.57),
- coefficient matrix in powers of `y`: (9.58),
- sixth-degree elimination polynomial `P(x)`: (9.59).

Compute all roots of `P(x)`, determine candidate Burmester axes, and recover associated moving axes for each valid solution.

[Requirements]
1. Use CAS as the primary solver and provide a Mathematica script (`.wls`) that runs end-to-end.
2. Do not assume formulas from memory; extract and follow Chapter 9 notation directly.
3. Ensure reproducibility:
   - fixed precision / solver settings,
   - deterministic ordering of roots and solutions,
   - no manual intervention during solve.
4. Report all solutions from the elimination pipeline:
   - all polynomial roots (real and complex),
   - real Burmester axes `G_i`,
   - corresponding moving axes `W1_i`,
   - multiplicities/degenerate cases if present.
5. Validate each candidate solution using residuals of the original constraints (9.25) across all five positions.
6. Assume only this prompt file and `Ch9_GDL.txt` exist in the workspace. Do not depend on any other pre-existing scripts, reports, or data files.
7. Autonomously debug until results are reproducible from a clean re-run.
8. Generate a `report.md` summarizing equations used, run commands, full solution tables, and validation outcomes.

[Deliverables]
1. A runnable Mathematica script (for example `ch9_rr_5pos_synthesis.wls`) implementing the full elimination workflow.
2. Reproducible output files (for example `.txt`, `.csv`, or `.json`) containing:
   - polynomial coefficients,
   - roots,
   - recovered axes,
   - residual checks.
3. A concise report file (for example `ch9_rr_5pos_report.md`) with:
   - chapter-equation traceability,
   - input orientation data,
   - all computed solutions,
   - final validation/pass-fail summary.
```
