# GDL LLM Reproduction Data (Prepared by Hai-Jun Su)

This workspace is a data-and-artifacts repository for reproducing linkage analysis and synthesis workflows from:

- **J. M. McCarthy and G. S. Soh**, *Geometric Design of Linkages*, Interdisciplinary Applied Mathematics 11, Springer, 2011.
- Book DOI: **10.1007/978-1-4419-7892-9**

The goal is to support the paper by **Haijun Su** by providing reproducible prompts, extracted chapter text, solver scripts, and generated outputs.

## What this repository contains

Each chapter folder contains a combination of:
- chapter text extracted from the digital book (`*_GDL.txt`),
- prompt files used to run autonomous coding workflows,
- implementation scripts (`.py` / `.wls`),
- reports (`report.md`),
- generated structured outputs (`.json`, `.csv`, `.txt`).

## Folder overview

### `ch2/` — Planar 4R position analysis
Key files:
- `Ch2_GDL.txt`: extracted chapter text (Chapter 2 content)
- `Ch2_4R_Prompt.md`: autonomous task prompt
- `Ch2_4R_solver.py`: planar 4R solver implementation
- `report.md`: run report and validation summary
- `Ch2_4R_Summary.txt`: compact equation/test/pass-fail summary

Status from existing artifacts:
- 4R solver runs successfully and reports machine-precision residual checks with PASS.

### `ch5/` — Planar RR/PR dyad synthesis
Key files:
- `Ch5_GDL.txt`: extracted chapter text (Chapter 5 content)
- `Ch5_RR_PR_Prompt.md`: autonomous task prompt
- `rr_pr_dyad_synthesis.py`: RR and PR synthesis implementation
- `report.md`: equations used, test data, residuals, and validation

Status from existing artifacts:
- RR and PR synthesis pipelines report PASS with residual/rank checks.

### `ch6/` — Multiloop planar linkage source text
Key files:
- `Ch6_GDL.txt`: extracted chapter text (Chapter 6 content)

This folder currently provides source chapter text for future reproduction tasks.

### `ch9/` — Spherical RR synthesis (4-position and 5-position)
Key files:
- `Ch9_GDL.txt`: extracted chapter text (Chapter 9 content)
- `Ch9_RR_4Position_Prompt.md`, `Ch9_RR_5Position_Prompt.md`: autonomous task prompts
- `ch9_rr_4pos_synthesis.wls`, `ch9_rr_5pos_synthesis.wls`: Mathematica scripts
- `ch9_rr_4pos_report.md`, `ch9_rr_5pos_report.md`: detailed run reports
- Output artifacts: `.json`, `.csv`, and polynomial/cone text files

Status from existing artifacts:
- 4-position workflow exports cone equations and validated sampled dyads.
- 5-position workflow produces a 6th-degree elimination polynomial and validated real solutions.

## Reproduction workflow

Each chapter directory includes a prompt file that can be used directly in an autonomous coding agent. Typical command-style prompt:

```text
read <PromptFile>.md and finish the task specified in the file
```

Then run the generated/maintained script in that folder, for example:

- Python workflows: `python <script>.py`
- Mathematica workflows: `wolframscript -file <script>.wls`

## Notes on extracted chapter text

The `*_GDL.txt` files are chapter extracts from the digital edition of *Geometric Design of Linkages* and are used as the equation source for reproduction.

Chapter-level DOI suffixes appear inside the extracted files (e.g., `_2`, `_5`, `_6`, `_9`) and map to the corresponding book chapters.

## Citation

If you use this dataset/reproduction workflow, cite the book and this repository.

### Book
McCarthy, J. M., & Soh, G. S. (2011). *Geometric Design of Linkages*. Interdisciplinary Applied Mathematics, 11. Springer. https://doi.org/10.1007/978-1-4419-7892-9

### Repository
Su, H. *gdl_llm_reproduction*. GitHub. https://github.com/haijunsu-osu/gdl_llm_reproduction
