# Codex output index

This folder contains the new illustrated versions and the supporting methodology/visual audit. The original blueprint files in the project root were not used as output targets and remain separate.

## Final PDFs

- `output/pdf/longitudinal_residual_rank_blueprint_concise_new.pdf` — 7-page illustrated concise blueprint, with 6 figures.
- `output/pdf/longitudinal_residual_rank_blueprint_detail_new.pdf` — 38-page illustrated detailed blueprint, with 12 figures.

## Editable sources

- `src/longitudinal_residual_rank_blueprint_concise_new.tex`
- `src/longitudinal_residual_rank_blueprint_detail_new.tex`
- `src/visual_style.tex` — shared visual language.
- `figures/` — reusable TikZ figure sources.

## Review record

- `notes/methodology_and_visual_audit.md` — methodology reconstruction, claim boundaries, visual gaps, and revision decisions.
- `render/final_concise/` and `render/final_detail/` — final rendered-page images used for visual QA.

## Final checks

- Both PDFs compile without LaTeX warnings, undefined references, or overfull boxes.
- Automated page inspection found no blank or near-empty pages.
- Key figure pages, title/contents pages, results dashboard, and final reference pages were visually inspected after the final build.
