# Planar 4R Benchmark Workspace

This workspace is set up to generate a planar 4R kinematics solver and validation report from a single prompt.

## Quick Start (Single Shot)
Use this exact prompt in Codex:

```text
read Ch2_4R_Prompt.md and finish the task specified in the file
```

This is intended to run end-to-end without additional human supervision.

## What the Agent Uses
- `Ch2_4R_Prompt.md`: task instructions and deliverables.
- `Ch2_GDL.txt`: source chapter (Geometric Design of Linkages, Chapter 2 content).

## Expected Outputs
After the one-shot run, you should have:
- `Ch2_4R_solver.py`: Python implementation of planar 4R position analysis.
- `report.md`: run report with commands, parameters, branch outputs, and validation status.
- `Ch2_4R_Summary.txt`: short plain-text summary of equations used, test case, and pass/fail.

## Manual Re-Run (Optional)
If you want to re-run the generated solver manually:

```powershell
python Ch2_4R_solver.py
```

You can also pass custom parameters, for example:

```powershell
python Ch2_4R_solver.py --a 1 --b 2 --g 2.5 --h 2 --theta-deg 20 40 80
```

## Workspace Files
- `Ch2_4R_Prompt.md` - master prompt specification.
- `Ch2_GDL.txt` - reference text to extract equations from.
- `Ch2_4R_solver.py` - generated solver script.
- `report.md` - execution and validation report.
- `Ch2_4R_Summary.txt` - short textual summary.
