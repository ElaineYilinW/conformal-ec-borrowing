#!/bin/bash
# Convert the blueprint's figure PDFs (in ../paper/figures) to PNGs for the deck.
cd "$(dirname "$0")" || exit 1
SRC="../paper/figures"
for f in fig_design fig_long fig_compat fig_rank fig_rankres fig_misspec \
         fig_pipeline fig_shift fig_cv fig_trimming fig_outcome fig_pval fig_snr; do
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r250 -dGraphicsAlphaBits=4 -dTextAlphaBits=4 \
     -sOutputFile="figs/$f.png" "$SRC/$f.pdf" 2>/dev/null && echo "ok $f" || echo "FAIL $f"
done
