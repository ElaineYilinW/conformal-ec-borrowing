# CC_ppt — presentation deck

A 30-slide talk on *Safe Selective Borrowing of External Controls via
Covariate-Shift-Adjusted Longitudinal Residual-Rank Screening*, built as a
**Google-Slides-compatible `.pptx`** (not Beamer). Every formula and every figure
is an **image** (LaTeX-compiled or from the blueprint pipeline) — no equations are
typed directly in the slides. Bullets are written verbose, to be read aloud.

## The file to open
**`CC_presentation.pptx`** — open it in Google Slides via **File ▸ Import slides**
(or drag it into Google Drive and open with Slides), or in PowerPoint/Keynote.

## What is in this folder
| item | what it is |
|---|---|
| `CC_presentation.pptx` | the deck (30 slides, 16:9) |
| `build_ppt.py` | the script that assembles the deck (all slide text + image placement) |
| `gen_eq.sh` | renders every formula / results table to a PNG (standalone LaTeX → `eq/*.png`) |
| `gen_figs.sh` | converts the blueprint's figure PDFs to PNG (`../fig_*.pdf` → `figs/*.png`) |
| `preview.py` | approximate visual preview of each slide (for layout checking) |
| `audit.py` | flags any text overflow or off-slide shape |
| `figs/` | the 13 figure PNGs used in the deck |
| `eq/`   | the 22 formula / table PNGs used in the deck |

## Structure of the talk
1. **Part I — Problem & background** (detailed): motivation · why longitudinal ·
   real-world ECs · identification & the four compatibility concepts · what the
   screen actually diagnoses · existing design-stage propensity screen.
2. **Part II — Method, step by step**: the basic idea · why rank beats a raw
   residual (distribution-free, and robust to *model misspecification*) · the full
   six-step method (scores · covariate-shift weighting · CV⁺ · symmetric +
   *adaptive* threshold · transported AIPW). The one-sided screen is omitted; the
   adaptive threshold is covered.
3. **Part III — Simulation**: the data-generating process (all math) · the screen
   in action · the signal-to-noise limit · results — only the good method plus the
   RCT-only / full-borrow / oracle anchors (poorly-performing variants and
   not-yet-run setups are omitted).
4. **Part IV — Discussion**: scope, limits, positioning · open problems · summary.

## To regenerate
```bash
bash gen_figs.sh      # figure PNGs from ../fig_*.pdf   (needs Ghostscript)
bash gen_eq.sh        # formula/table PNGs               (needs LaTeX + Ghostscript)
python3 build_ppt.py  # assemble CC_presentation.pptx    (needs python-pptx + Pillow)
python3 audit.py      # optional: overflow / bounds check
```
