# Presentation deck

The ISEF slide deck, published as a Claude artifact at
<https://claude.ai/code/artifact/a4675615-62ab-4587-ab77-6b81ee6f9ba7>.

**`deck.html` is generated. Edit `tpl.html`, never `deck.html`.**

## Build

```sh
python3 mk3d.py      # optional: re-render the slide-07 backbone traces from the PDBs
python3 build.py     # tpl.html + *.json -> deck.html
```

```sh
python3 export.py    # deck.html -> PPIL4_deck.pdf + PPIL4_deck.pptx
```

`build.py` substitutes `__FIG_<name>__` (inline SVG, from `figs.json`) and
`__PNG_<name>__` (base64 data URIs, from `png.json`), and injects the strip-chart
data from `pts.json` at `__DATA__`. It fails loudly on any unresolved placeholder.

## Exporting

`export.py` produces both formats and needs no third-party library.

* **PDF** is printed straight from Chrome, so its text stays vector and
  selectable. Printing needs three things the on-screen deck does not do:
  every slide laid out as its own page-sized block, `print-color-adjust:
  exact` so panel and note backgrounds actually reproduce, and the entrance
  animation suppressed so a half-played frame is never captured.
* **PPTX** is one full-bleed 2560x1440 image per slide. PowerPoint cannot
  represent this deck's inline SVG and CSS layout natively, so anything else
  would be a lossy re-drawing rather than the deck. The OOXML is written
  directly - `python-pptx` is not a dependency.

Both land in `deck/` and are gitignored; regenerate rather than commit them.

## Files

| file | what it is |
|---|---|
| `tpl.html` | the deck itself — one `<section class="slide">` per page, plus all CSS and JS |
| `build.py` | placeholder substitution; writes `deck.html` |
| `mk3d.py` | reads `../CRBN-Thalidomide-SALL4_(Ryan).pdb` and `../PPIL4_alphafold_(Ryan).pdb`, projects the backbones onto their principal axes and writes depth-sorted SVG traces into `figs.json` |
| `figs.json` | inline SVG fragments (2D molecule depictions, backbone traces) |
| `png.json` | base64 JPEGs of the rendered ternary complexes |
| `pts.json` | the 27 locked-protocol docking results that drive the strip chart |
| `export.py` | `deck.html` -> PDF (vector, via Chrome print) and PPTX (image per slide, hand-written OOXML) |

## Publishing

Publish `deck.html` to the artifact URL above. Publishing without that URL creates
a *new* artifact rather than updating the existing one.

## Conventions

* Colours come from CSS custom properties only (`--accent`, `--s-ctrl`, `--s-cand`,
  `--s-cross`, `--muted`, `--rule` …) so every figure survives the light/dark switch.
  Never hard-code a hex value in a figure.
* The stage is a fixed 1280x720 box scaled to fit; content must stay above y=685
  or it collides with the footer.
