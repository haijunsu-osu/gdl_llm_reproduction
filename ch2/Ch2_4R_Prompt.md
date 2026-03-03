# Preliminary Prompt: Planar 4R Kinematic Analysis

```text
[Context]
You are an expert in kinematics and mechanism design. Work like a careful human engineer: read the source material, extract the governing equations, implement them in code, and verify results with explicit checks.

[Task]
Read the provided `Ch2_GDL.txt` chapter file in the workspace (the chapter on planar linkage analysis from Geometric Design of Linkages). From that text, identify the formulation for planar 4-bar (4R) kinematic analysis, then write a complete Python script that:
- computes the linkage configuration from given link dimensions and input angle(s),
- solves for the unknown angle(s) using the chapter’s equations and conventions,
- and reports all physically valid solution branches.

[Requirements]
1. Do not assume formulas from memory; first read the provided `.txt` and follow the derivation given there.
2. Implement clean, runnable Python code (`.py`) with clear functions for:
   - equation setup,
   - solver logic,
   - branch handling,
   - and residual/error checking.
3. Include robust numerical handling:
   - angle wrapping,
   - inverse-trig domain clipping,
   - tolerance-based validation,
   - and clear behavior for non-assembly or singular cases.
4. Validate the script with at least one numerical example:
   - use a numerical example from the chapter if available,
   - otherwise create a reasonable test case and state that it is user-defined.
5. Autonomously debug until the script runs successfully and the validation checks pass.
6. As the final requirement, generate a `report.md` file summarizing:
   - the command(s) run,
   - the input/link parameters and test angles,
   - key branch outputs,
   - and the final validation status.

[Deliverables]
1. One working Python script containing:
   - the solver implementation,
   - a reproducible numerical example,
   - and printed validation results.
2. A short summary in plain text that includes:
   - which equations/sections from the `.txt` were used,
   - what test case was run,
   - and whether validation passed.
```
