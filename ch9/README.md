# Chapter 9 Spherical RR Synthesis Codespace

This workspace contains reproducible Mathematica (`.wls`) workflows for:
- 4-position spherical RR dyad synthesis (Chapter 9, Section 9.4.5 context)
- 5-position spherical RR dyad synthesis (Chapter 9, Section 9.4.6 elimination workflow)

The source chapter text is:
- `Ch9_GDL.txt`

The autonomous task prompts are:
- `Ch9_RR_4Position_Prompt.md`
- `Ch9_RR_5Position_Prompt.md`

## Zero-Supervision Reproduction (Codex Prompting)

In Codex, run exactly these prompts (one at a time):

```text
read Ch9_RR_4Position_Prompt.md and finish the task specified in the file
```

```text
read Ch9_RR_5Position_Prompt.md and finish the task specified in the file
```

No manual editing is required. The agent should read `Ch9_GDL.txt`, implement/run the solver, validate residuals, and emit artifacts plus `report.md`.

## Expected Outputs

### 4-position run
- `ch9_rr_4pos_synthesis.wls`
- `ch9_rr_4pos_results.json`
- `ch9_rr_4pos_samples.csv`
- `ch9_rr_4pos_center_cone.txt`
- `ch9_rr_4pos_circling_cone.txt`
- `ch9_rr_4pos_report.md`

### 5-position run
- `ch9_rr_5pos_synthesis.wls`
- `ch9_rr_5pos_results.json`
- `ch9_rr_5pos_roots.csv`
- `ch9_rr_5pos_root_groups.csv`
- `ch9_rr_5pos_solutions.csv`
- `ch9_rr_5pos_Px.txt`
- `ch9_rr_5pos_report.md`

## Direct Script Re-run (Optional)

If you want to rerun generated scripts directly:

```powershell
wolframscript -file ch9_rr_4pos_synthesis.wls
wolframscript -file ch9_rr_5pos_synthesis.wls
```

## Notes on Determinism

- Fixed precision and deterministic sorting are used in both scripts.
- No randomness is used.
- If full numeric task data is not explicitly present in `Ch9_GDL.txt`, the scripts use deterministic synthetic orientation sets and document this in their reports/JSON outputs.

## Expected Upper Bounds

- 4-position spherical RR synthesis:
  - The solution is a one-parameter family (center-axis/circling-axis cubic cone curves), so there is no finite upper bound on the number of dyads.
  - Practically: infinitely many valid dyads along the curve/cone.

- 5-position spherical RR synthesis:
  - The elimination step yields a sixth-degree polynomial `P(x)`, so the upper bound is 6 discrete dyad solutions (counting complex roots with multiplicity).
  - Real solution counts are typically `0`, `2`, `4`, or `6`.
