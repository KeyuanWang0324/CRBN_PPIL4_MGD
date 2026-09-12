# CRBN_Project_site

**This directory is a working copy. The live site is a different repository.**

| | |
|---|---|
| served from | <https://github.com/KeyuanWang0324/CRBN_Project> (branch `main`, path `/`) |
| live at | <https://keyuanwang0324.github.io/CRBN_Project/> |

The contents of this directory are the root of that repository. Committing
here does **not** publish anything — this repo is the research project, that
one is the website. They drift: at the time of writing, the live site was
three commits behind this directory.

## Publishing

```sh
git clone https://github.com/KeyuanWang0324/CRBN_Project.git /tmp/site
cp index.html glossary.js viewer.js structures.json /tmp/site/
cp -R img glb lib /tmp/site/
cd /tmp/site && git add -A && git commit && git push origin main
```

Copy files explicitly rather than syncing the whole tree: the site repo has
its own `README.md` that does not exist here, and `rsync --delete` would
remove it. GitHub Pages redeploys on push; allow a minute or two.

## Where the figures come from

The mechanism figures (`figure.mechfig`) are ported from the slide deck in
`../deck/` and translated to English. The deck is the upstream — if a figure
changes there, it has to be re-translated and re-pasted here, then published
to the site repo as above.

`.mechfig` redeclares the deck's figure palette (`--s-cand`, `--s-ctrl`,
`--warn-bg`) in all three theme states so the artwork survives the light/dark
switch without touching the site's global tokens.

## glossary.js

`annotate()` must never descend into an SVG: wrapping a term in an HTML
`<span>` there makes the text disappear, because SVG will not render it.
`SKIP_TAGS` cannot catch this on its own — SVG elements report lower-case
tagNames (`text`, `tspan`), never `SVG` — so the guard is `p.ownerSVGElement`.
