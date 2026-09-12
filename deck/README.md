# Presentation deck

The ISEF slide deck, published as a Claude artifact at
<https://claude.ai/code/artifact/a4675615-62ab-4587-ab77-6b81ee6f9ba7>.

**`deck.html` is generated. Edit `tpl.html`, never `deck.html`.**

## Build

```sh
python3 mk3d.py      # optional: re-render the slide-07 backbone traces from the PDBs
python3 build.py     # tpl.html + *.json -> deck.html
```

`build.py` substitutes `__FIG_<name>__` (inline SVG, from `figs.json`) and
`__PNG_<name>__` (base64 data URIs, from `png.json`), and injects the strip-chart
data from `pts.json` at `__DATA__`. It fails loudly on any unresolved placeholder.

## Files

| file | what it is |
|---|---|
| `tpl.html` | the deck itself — one `<section class="slide">` per page, plus all CSS and JS |
| `build.py` | placeholder substitution; writes `deck.html` |
| `mk3d.py` | reads `../CRBN-Thalidomide-SALL4_(Ryan).pdb` and `../PPIL4_alphafold_(Ryan).pdb`, projects the backbones onto their principal axes and writes depth-sorted SVG traces into `figs.json` |
| `figs.json` | inline SVG fragments (2D molecule depictions, backbone traces) |
| `png.json` | base64 JPEGs of the rendered ternary complexes |
| `pts.json` | the 27 locked-protocol docking results that drive the slide-11 strip chart |

## Publishing

Publish `deck.html` to the artifact URL above. Publishing without that URL creates
a *new* artifact rather than updating the existing one.

## Conventions

* Colours come from CSS custom properties only (`--accent`, `--s-ctrl`, `--s-cand`,
  `--s-cross`, `--muted`, `--rule` …) so every figure survives the light/dark switch.
  Never hard-code a hex value in a figure.
* The stage is a fixed 1280x720 box scaled to fit; content must stay above y=685
  or it collides with the footer.
