#!/usr/bin/env bash
# Convert a revtex .tex to .docx with correctly-rendered Word equations.
# Usage: ./tex2docx.sh methods.tex   ->   methods.docx
#
# Handles two pandoc/OMML quirks:
#  1. \Bigl..\Bigr are not converted -> rewritten to \left..\right.
#  2. Bare operator scripts (^*, ^+, lone $+$, trailing ,+}) render as "¿" in
#     Word/LibreOffice -> wrapped in \text{} so they display as plain symbols.
# The source .tex is left untouched; edits are made on a temp copy only.
# PDF \includegraphics figures are rasterized to PNG (Word embeds PNG, not PDF).
set -euo pipefail

src="$1"
[ -f "$src" ] || { echo "no such file: $src" >&2; exit 1; }
base="$(basename "$src" .tex)"
dir="$(cd "$(dirname "$src")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

work="$tmp/$base.tex"
here="$(cd "$(dirname "$0")" && pwd)"

# --- equation numbering + reference rewriting (see number_refs.py) ---
python3 "$here/number_refs.py" < "$src" > "$tmp/numbered.tex"

# --- math fixes (see header) ---
sed -e 's/\\Bigl/\\left/g' -e 's/\\Bigr/\\right/g' \
    -e 's/\\bigl/\\left/g' -e 's/\\bigr/\\right/g' \
    -e 's/,+}/,\\text{+}}/g' \
    -e 's/\^{+}/^{\\text{+}}/g' \
    -e 's/\^\*/^{\\text{*}}/g' \
    -e 's/\^+/^{\\text{+}}/g' \
    -e 's/\$+\$/$\\text{+}$/g' \
    "$tmp/numbered.tex" > "$work"

# --- rasterize PDF figures referenced by \includegraphics (if any) ---
mapfile -t figs < <(grep -oE '\\includegraphics(\[[^]]*\])?\{[^}]*\.pdf\}' "$src" \
                    | sed -E 's/.*\{([^}]*)\}/\1/' | sort -u || true)
for fig in "${figs[@]:-}"; do
    [ -n "$fig" ] || continue
    abspdf="$dir/$fig"; [ -f "$abspdf" ] || abspdf="$fig"
    [ -f "$abspdf" ] || { echo "warn: figure not found: $fig" >&2; continue; }
    png="$tmp/$(basename "${fig%.pdf}").png"
    pdftoppm -png -r 300 -singlefile "$abspdf" "${png%.png}"
    sed -i "s#$fig#$png#g" "$work"
done

bib=""; [ -f "$dir/smelling.bib" ] && bib="--bibliography=$dir/smelling.bib --citeproc"
pandoc "$work" -o "$dir/$base.docx" --number-sections $bib
echo "wrote $dir/$base.docx"
