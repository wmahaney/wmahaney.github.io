# wmahaney.github.io

Personal academic website of William E. Mahaney, built with [Quarto](https://quarto.org) and
published with GitHub Pages. Live site: <https://wmahaney.github.io>

## Where things go

| Folder / file | What it is | What to do |
|---|---|---|
| `index.qmd` | Home page: intro, isogeny graphs, news, about | Edit text, links, news items |
| `research/index.qmd` | Research page: papers and talks | Add papers and talks (templates inside the file) |
| `writing/index.qmd` | **Writing** page: lists every note and post | Leave alone |
| `notes/articles/<name>/index.qmd` | Expository notes, one folder per note | Copy `notes/_template.qmd` |
| `blog/posts/<name>/index.qmd` | Blog posts, one folder per post | Copy `blog/_template.qmd` |
| `teaching/index.qmd` | Teaching page | Fill in the hidden skeleton (instructions inside the file) |
| `cv/index.qmd` | CV page | Edit text |
| `cv_william_mahaney.pdf` | The CV PDF (top-level, so the old link keeps working) | Replace with a file of the same name |
| `images/profile/` | **Upload your profile photo here** as `profile.jpg` | Upload |
| `images/notes/`, `images/blog/` | **Upload figures** for notes and posts | Upload |
| `files/papers/` | **Upload paper PDFs** | Upload |
| `files/talks/` | **Upload slides** | Upload |
| `files/teaching/` | **Upload syllabi and handouts** | Upload |
| `files/cv/` | **Upload research / teaching / personal statements** | Upload |
| `assets/` | Colors, fonts, and page styling | Rarely needed |
| `assets/js/`, `assets/data/` | Script and data for the home-page isogeny graphs | See below |
| `_scripts/isogeny_graphs.py` | Computes the isogeny graph data | Run to change the primes |
| `images/social-card.png` | Picture shown when a link to the site is shared | Replace with a 1200 x 630 image |
| `_quarto.yml` | Site settings: menu, footer, theme | Rarely needed |
| `.github/workflows/publish.yml` | Builds and publishes the site automatically | Leave alone |

Every upload folder contains a file called `_UPLOAD_HERE.md` explaining exactly what goes there.
Files starting with `_` (such as those notes and the templates) are never published.

## Publishing

Pushing to `main` builds the site with GitHub Actions and publishes it (give it a minute or two).
Pushing any other branch builds it too, without publishing, so mistakes show up before they go live.
The one-time switch for this: **Settings → Pages → Build and deployment → Source: GitHub Actions**.

## Previewing on your own computer (optional)

Install [Quarto](https://quarto.org/docs/get-started/), then run `quarto preview` in this folder.

## Filling in hidden skeletons

Some pages contain ready-made pieces that are hidden until you fill them in. They look like this:

```
<!-- TO FILL IN: ... then delete this comment line and the "END" line.
...the piece...
END -->
```

Edit the piece, then delete the first line (`<!-- TO FILL IN ...`) and the last line (`END -->`).

## The isogeny graphs on the home page

The home page draws three graphs for a choice of primes `p`: the supersingular 2-isogeny graph,
the (2,3)-isogeny graph, and the 2-isogeny graph with Borel level-3 structure. The data is
computed ahead of time by `_scripts/isogeny_graphs.py` (plain Python, no SageMath needed) and
stored in `assets/data/isogeny-graphs.json`. To offer different primes:

```
python3 _scripts/isogeny_graphs.py 101 431 1009
```

then commit the new JSON file. The script checks its own output on every prime (class number
formula, modular polynomials, the mass formula, Brandt-matrix symmetry, and the Ramanujan bound)
and stops if a check fails. The method is described at the top of the script.

## Changing the look

Colors: `assets/theme-light.scss` and `assets/theme-dark.scss` (the `--wm-*` values near the bottom).
Layout and cards: `assets/styles.css`. Menu font: `assets/head.html`.
