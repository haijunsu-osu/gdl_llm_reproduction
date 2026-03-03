# Preliminary Prompt: RR and PR Dyad Synthesis

```text
[Context]
You are an expert in kinematics and mechanism design. Work like a careful human engineer: read the source material, extract the governing synthesis equations, implement them in code, and verify results with explicit checks.

[Task]
Read the provided `Ch5_GDL.txt` chapter file in the workspace. From that text, identify the formulation for synthesis of:
- planar RR dyads, and
- planar PR dyads.

Then write a complete Python solution that handles both synthesis problems and reports all physically valid solutions for each dyad type.

[Requirements]
1. Do not assume formulas from memory; first read `Ch5_GDL.txt` and follow the derivations and notation in that chapter.
2. Implement clean, runnable Python code (`.py`) with clear functions for:
   - RR dyad equation setup and solver,
   - PR dyad equation setup and solver,
   - branch/solution-set handling,
   - and residual/error checking.
3. Include robust numerical handling:
   - tolerance-based rank/degeneracy checks,
   - inverse-trig/domain safety where applicable,
   - angle normalization/consistency handling,
   - and clear behavior for infeasible or singular cases.
4. Validate both synthesis pipelines with numerical examples:
   - use chapter examples if available,
   - otherwise create reasonable user-defined examples and label them clearly.
5. Autonomously debug until both RR and PR workflows run successfully and validation checks pass.
6. Generate a `report.md` file summarizing:
   - commands run,
   - equations/sections used from `Ch5_GDL.txt`,
   - input design data for RR and PR tests,
   - key computed synthesis results for each,
   - residual/error metrics,
   - and final pass/fail validation status.

[Deliverables]
1. One working Python script (or a small module) that includes:
   - RR dyad synthesis solver,
   - PR dyad synthesis solver,
   - reproducible numerical validations for both.
2. A `report.md` document covering both synthesis problems with concise technical conclusions.
3. Printed console output showing intermediate and final results for RR and PR cases.
```

