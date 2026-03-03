# 4R Solver Run Report

## Run Metadata
- Date/Time: 2026-03-03 11:23:03 -05:00
- Script: `Ch2_4R_solver.py`
- Command: `python Ch2_4R_solver.py`

## Test Configuration
- Link lengths: `a=1.0, b=2.0, g=2.5, h=2.0`
- Input angles (deg): `20, 40, 80, 120, 160`
- Validation tolerance: `1e-10`
- Example type: user-defined (no explicit numeric 4R example found in section 2.3)

## Summary Results
- Solver status for all tested input angles: `ok`
- Physically valid branches per input angle: `2`
- Overall validation: `PASS`

| theta (deg) | status | branches | max \|(2.50)\| | max \|loop_x\| | max \|loop_y\| | max \|length\| |
|---:|:---|---:|---:|---:|---:|---:|
| 20.000000 | ok | 2 | 3.553e-15 | 4.441e-16 | 1.110e-15 | 4.441e-15 |
| 40.000000 | ok | 2 | 1.776e-15 | 4.441e-16 | 4.441e-16 | 1.776e-15 |
| 80.000000 | ok | 2 | 1.776e-15 | 2.220e-16 | 3.331e-16 | 8.882e-16 |
| 120.000000 | ok | 2 | 3.553e-15 | 6.661e-16 | 7.772e-16 | 3.553e-15 |
| 160.000000 | ok | 2 | 3.553e-15 | 9.992e-16 | 6.661e-16 | 4.441e-15 |

## Branch Angles
| theta (deg) | branch | psi (deg) | phi (deg) | assembly | valid |
|---:|---:|---:|---:|:---|:---:|
| 20 | 1 | -125.900535 | -98.826900 | B_below_AC | True |
| 20 | 2 | +101.173100 | +34.099465 | B_above_AC | True |
| 40 | 1 | -137.876700 | -122.803384 | B_below_AC | True |
| 40 | 2 | +97.196616 | +2.123300 | B_above_AC | True |
| 80 | 1 | -152.109094 | -153.779420 | B_below_AC | True |
| 80 | 2 | +106.220580 | -52.109094 | B_above_AC | True |
| 120 | 1 | -157.419926 | -174.784301 | B_below_AC | True |
| 120 | 2 | +125.215699 | -97.419926 | B_above_AC | True |
| 160 | 1 | -155.465796 | +164.108905 | B_below_AC | True |
| 160 | 2 | +144.108905 | -135.465796 | B_above_AC | True |

## Conclusion
The generated 4R solver executed successfully and produced two valid assembly branches for each tested input angle. All equation and loop residuals are at machine-precision scale (`~1e-15`), consistent with a numerically correct implementation.
