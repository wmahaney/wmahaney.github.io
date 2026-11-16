# wmahaney.github.io

Personal academic website of William E. Mahaney, built with [Quarto](https://quarto.org) and
published with GitHub Pages. Live site: <https://wmahaney.github.io>

## Where things go

| Folder / file | What it is | What to do |
|---|---|---|
| `index.qmd` | Home page: intro, isogeny graphs, news, about | Edit text, links, news items |
| `research/index.qmd` | Research page: papers and talks | Add papers and talks (templates inside the file) |
| `writing/index.qmd` | **Writing** page: lists every note and post together | Leave alone |
| `blog/index.qmd`, `notes/index.qmd` | The Blog and Expository Notes pages in the Writing menu | Leave alone |
| `notes/articles/<name>/index.qmd` | Expository notes, one folder per note | Copy `notes/_template.qmd` |
| `blog/posts/<name>/index.qmd` | Blog posts, one folder per post | Copy `blog/_template.qmd` |
| `teaching/index.qmd` | Teaching page | Fill in the hidden skeleton (instructions inside the file) |
| `statements/<name>/index.qmd` | **Statements**: research, teaching, AI, culture | Empty pages for your own text; see "Statements and AI pages" |
| `ai/index.qmd`, `ai/projects/<name>/index.qmd` | **AI** page and one folder per AI project | Copy `ai/projects/_template.qmd` |
| `cv/index.qmd` | CV page, written out in HTML | Edit entries (templates inside the file); keep the PDF in step |
| `cv_william_mahaney.pdf` | The CV PDF (top-level, so the old link keeps working) | Replace with a file of the same name |
| `images/profile/` | **Upload your profile photo here** as `profile.jpg` | Upload |
| `images/notes/`, `images/blog/` | **Upload figures** for notes and posts | Upload |
| `files/papers/` | **Upload paper PDFs** | Upload |
| `files/talks/` | **Upload slides** | Upload |
| `files/teaching/` | **Upload syllabi and handouts** | Upload |
| `files/statements/` | **Upload statement PDFs** | Upload |
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

## Statements and AI pages; the private branch

The four statements (`statements/<name>/index.qmd`: research, teaching, culture, AI) are empty
and marked `draft: true`. The public site leaves them out, and its menu has no Statements tab.

The site has two build profiles:

| Profile | File | Used for | Drafts | Statements menu |
|---|---|---|---|---|
| public (default) | `_quarto-public.yml` | `main` → https://wmahaney.github.io | hidden | hidden |
| private | `_quarto-private.yml` | the `private` branch; `quarto preview --profile private` | shown | shown |

**Test-deploying sensitive updates.** Work on the `private` branch. Every push builds the
site with the private profile and attaches it as **site-private** (Actions → the run →
Artifacts). Unzip it and open `index.html`. The `private` branch is never deployed.

**This repository is public.** Anything pushed to any branch, including `private` (the
`.qmd` source and the site-private download), can be read by anyone on GitHub. "Private"
here means "not on the live site", not "confidential". Keep confidential text out of this
repository; use `quarto preview --profile private` on your own computer, or a private repository.

**Publishing a statement:** delete its `draft: true`, then copy the "Statements" menu entry
from `_quarto-private.yml` into `_quarto-public.yml` (listing only the statements that are
published), and merge into `main`. The AI tab (`ai/index.qmd`) and the Statements overview
page (`statements/index.qmd`) are drafts too, and work the same way. A home-page card for
either can be added to the card row in `index.qmd` once it is public.

## The isogeny graphs on the home page

The home page draws three graphs for a choice of primes `p`: the supersingular 2-isogeny graph,
the (2,3)-isogeny graph, and the 2-isogeny graph with Borel level-3 structure. The data is
computed ahead of time by `_scripts/isogeny_graphs.py` (plain Python, no SageMath needed) and
stored as one small file per prime in `assets/data/isogeny/` (`p431.json`, ... plus `index.json`),
so the page only downloads the prime being viewed. The default list goes up to p = 4999; the level-3
graph is offered up to p = 5000. The walk has "history" (10 to 250 steps) and "speed" controls. To offer different primes:

```
python3 _scripts/isogeny_graphs.py 101 431 1009
```

then commit the new JSON files. The script checks its own output on every prime (class number
formula, modular polynomials, the mass formula, Brandt-matrix symmetry, and the Ramanujan bound)
and stops if a check fails. The method is described at the top of the script.

## Changing the look

Colors: `assets/theme-light.scss` and `assets/theme-dark.scss` (the `--wm-*` values near the bottom).
Layout and cards: `assets/styles.css`. Menu font: `assets/head.html`.
