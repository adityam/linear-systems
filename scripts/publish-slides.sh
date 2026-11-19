#!/bin/sh
# Post-render hook (see _quarto-website.yml).
#
# The lecture slides are a separate Quarto project in slides/, which the book
# render skips. Render them and copy them into the book's output directory, so
# that they are published next to the notes at slides/<deck>.html.
set -e

# Full renders of the book only (as in CI), not single-chapter renders
[ "${QUARTO_PROJECT_RENDER_ALL:-}" = "1" ] || exit 0

quarto render slides

out="$QUARTO_PROJECT_OUTPUT_DIR/slides"
mkdir -p "$out"
for deck in slides/*.html; do
  cp "$deck" "$out/"
  cp -R "${deck%.html}_files" "$out/"
  # Figures used by the deck (slides/figures is a link to ../figures). The book
  # inlines its GIFs as base64, so they are not in the output; copy them here.
  grep -o '"figures/[^"]*"' "$deck" | tr -d '"' | sort -u | while read -r f; do
    mkdir -p "$out/$(dirname "$f")"
    cp "slides/$f" "$out/$f"
  done
done
