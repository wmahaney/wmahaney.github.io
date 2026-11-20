# wmahaney.github.io

Personal academic website of William E. Mahaney, built with [Quarto](https://quarto.org) and
published with GitHub Pages. Live site: <https://wmahaney.github.io>

This file is the guide to editing the site. Each section says which file to open and what to
change. Most files also carry their own instructions as comments at the top.

- [How publishing works](#how-publishing-works)
- [Before you edit: four things to know](#before-you-edit-four-things-to-know)
- [Home page](#home-page)
- [Menu bar and profile links](#menu-bar-and-profile-links)
- [Research: papers, talks, slides](#research-papers-talks-slides)
- [Writing: notes and blog posts](#writing-notes-and-blog-posts)
- [Code pages](#code-pages)
- [CV](#cv)
- [Teaching](#teaching)
- [Statements, the AI page, and the private branch](#statements-the-ai-page-and-the-private-branch)
- [Animations and the isogeny graphs](#animations-and-the-isogeny-graphs)
- [Look: colours, fonts, layout](#look-colours-fonts-layout)
- [Where everything is](#where-everything-is)

## How publishing works

Pushing to `main` builds the site with GitHub Actions and publishes it. Give it a minute or two.
Pushing any other branch builds it too, without publishing, so a mistake shows up before it is live.
The workflow is `.github/workflows/publish.yml`. On every push it:

1. compiles the CV from `cv/cv.tex`, if that file exists (see [CV](#cv));
2. runs the tests of the random-walk code (`node --test _tests/*.test.mjs`);
3. runs `quarto render`, which writes the site to `_site/`;
4. on `main` only, publishes `_site/`.

If any step fails, nothing is published and the live site stays as it was. On a branch other than
`main`, the built site is attached to the run as **site-preview** (Actions → the run → Artifacts);
unzip it and open `index.html` to look at it before merging.

The one-time switch for this: **Settings → Pages → Build and deployment → Source: GitHub Actions**.

To preview on your own computer, install [Quarto](https://quarto.org/docs/get-started/) and run
`quarto preview` in this folder.

## Before you edit: four things to know

**1. The site shows what exists.** Posted papers, talks that were given, tools that are live,
released code. It does not announce work that is on the way: no "in preparation", no counts of
forthcoming papers, no "coming soon". When in doubt, leave it out.

**2. Hidden skeletons.** Some pages contain ready-made pieces that are hidden until you fill them
in. They look like this:

```
<!-- TO FILL IN: ... then delete this comment line and the "END" line.
...the piece...
END -->
```

Edit the piece, then delete the first line (`<!-- TO FILL IN ...`) and the last line (`END -->`).
Everything between `<!--` and `-->` is left out of the published pages.

**3. Drafts.** A page whose header says `draft: true` is not published and does not appear in
menus or listings. Delete that line to publish the page.

**4. Underscores.** A file or folder whose name starts with `_` is never published: the
`_UPLOAD_HERE.md` notes, the templates, `_partials/`, `_scripts/`, `_tests/`.

## Home page

File: `index.qmd`.

| To change | Edit |
|---|---|
| The line above your name ("Number Theory · Cryptography · Probability"), your name, your position | The block near the top of `index.qmd` that starts with `<p class="wm-eyebrow">` |
| The introduction sentence | The paragraph inside `::: {.wm-bio}` |
| The job-search line | `_partials/seeking.qmd` holds the sentence. To remove it, delete the line `{{< include /_partials/seeking.qmd >}}` in `index.qmd` |
| The buttons (CV, Email, GitHub) | The block `::: {.wm-actions .wm-profile-links}`. arXiv and ORCID are added automatically; see [Menu bar and profile links](#menu-bar-and-profile-links) |
| News | The list inside `::: {.wm-news}`. Newest first, three or four items. Each item is `- [date]{.wm-news-date} text`. Delete the whole block to hide News |
| The four cards (Research, Writing, Code, Teaching) | The `<section class="wm-cards">` block. Each card is one `<a class="wm-card" ...> ... </a>`. Icon names come from [Bootstrap Icons](https://icons.getbootstrap.com) |
| About me | The paragraphs inside `::: {.wm-about-text}` |
| Profile photo (the circle beside your name) | Upload `images/profile/profile.jpg` (see the note in that folder). Until it exists, the circle shows the initials WEM |
| The isogeny graphs | They sit in a side panel, opened by the "Isogeny graphs" tab on the right edge of the page (a button at the bottom right on a phone). The panel is the block at the end of `index.qmd`; it opens and closes with `assets/js/side-panel.js`, and its look is under "side panel" in `assets/styles.css`. On a screen at least 1600px wide the page moves over to make room; on a narrower one the panel covers it. For the graph itself see [Animations and the isogeny graphs](#animations-and-the-isogeny-graphs) |

## Menu bar and profile links

**The menu entries** are in two files with the same layout:

- `_quarto-public.yml`: the site everyone sees (used for `main`);
- `_quarto-private.yml`: the private test build, which also has the Statements menu.

Change both when you add, remove or reorder a menu entry. An entry is two lines, `href:` (the
page) and `text:` (the label). A dropdown has `text:` and a `menu:` list of entries.

**The Email and GitHub icons** at the right of the menu are under `tools:` in `_quarto.yml`. The
`text:` of each is what a screen reader reads and what shows when the pointer rests on the icon.

**arXiv and ORCID.** Open `_quarto.yml`, find `profile-links:` and fill in the ids:

```yaml
profile-links:
  arxiv-author: "surname_x_1"
  orcid: "0000-0000-0000-0000"
```

Give the id only, not the full address. A link then appears in the menu bar, on the home page
and on the CV page. While an id is empty, its link is not shown anywhere. (The code is
`assets/profile-links.lua`.)

**The footer** is `page-footer:` in `_quarto.yml`. **The description and picture shown when a link
to the site is shared** are `description:` and `image:` in the same file; the picture is
`images/social-card.png` (1200 × 630).

## Research: papers, talks, slides

File: `research/index.qmd`.

**Add a paper.** Copy the block under `PAPER TEMPLATE` (in the comment at the top of the Papers
section), paste it at the top of the list, and edit it. List all authors alphabetically by
surname. Link the posted paper (arXiv or ePrint); for a PDF hosted here, put it in
`files/papers/` and link it as `[PDF](../files/papers/your-file.pdf){.wm-chip}`. Only posted
papers go on the page.

**Fill in a summary or abstract.** Each paper has a hidden `TO FILL IN` piece for each; see
"Hidden skeletons" above.

**Add a talk.** Talks are a list in `research/talks.yml`, newest first. Copy an entry (four
lines), paste it at the top, and edit it. Only talks that have been given go in the list.

**Add slides to a talk.**

1. Upload the PDF to `files/talks/`, named `year-meeting.pdf` (for example `2025-mason-vii.pdf`).
2. In `research/talks.yml`, put the file name in that talk's `slides:` line:
   `slides: "2025-mason-vii.pdf"`.

A "Slides" link then appears next to the talk. A talk with `slides: ""` shows no link.

## Writing: notes and blog posts

- **A new expository note:** make a folder `notes/articles/<short-name>/`, copy
  `notes/_template.qmd` into it as `index.qmd`, and edit it.
- **A new blog post:** the same with `blog/posts/<short-name>/` and `blog/_template.qmd`.

The card for it appears automatically on the Writing, Expository Notes or Blog page. Figures go
in `images/notes/` or `images/blog/`. To keep a half-written piece private, add `draft: true` to
its header. Math is written as in LaTeX: `$x^2$` inline, `$$ ... $$` displayed. The listing
pages themselves (`writing/index.qmd`, `notes/index.qmd`, `blog/index.qmd`) rarely need editing.

## Code pages

One data file, `code/projects.yml`, drives both Code pages:

- `code/index.qmd`: the overview, with the list of repositories;
- `code/basic-tools.qmd`: one section per browser tool.

**Add a repository.** Open `code/projects.yml`, copy one of the commented templates at the
bottom, remove the `# ` at the start of each line, and edit it. An entry with `status: live` is
shown; one with `status: hidden` is left out of the site entirely. While no entry is live, the
overview shows its introduction only, with no empty heading.

**Switch a browser tool on**, once its page loads at its `url`:

1. in `code/projects.yml`, change its `status: hidden` to `status: live`;
2. in `_quarto-public.yml` and `_quarto-private.yml`, add an entry for it to the `Code` menu,
   next to the cropmark and travel_notes entries.

To take a tool off the site, set it back to `status: hidden` and delete its menu entry.

**Move a tool** (a renamed repository, a custom domain): change its `url:` and `source:` in
`code/projects.yml`. The address is written nowhere else.

**Change what is said about a tool:** its `about:`, `usage:` and `limits:` in the same file. The
layout of the pages is in `code/_tools.ejs` and `code/_repositories.ejs`. The introduction of the
overview is the paragraph in `code/index.qmd`.

## CV

The CV is written once, in LaTeX. Put the source at `cv/cv.tex` (with any `.cls`, `.sty` or
`.bib` it needs next to it). To update the CV: edit `cv/cv.tex` and push. The site rebuilds, and
the CV page, the "Curriculum vitae" button on the home page and every link to
`cv_william_mahaney.pdf` get the new PDF.

Until `cv/cv.tex` exists, the site serves the `cv_william_mahaney.pdf` that is in the repository;
to update the CV in the meantime, replace that file. Details, building locally and using
another LaTeX engine: `cv/README.md`.

The CV page (`cv/index.qmd`) shows the PDF and has no text of its own to keep in step.

## Teaching

File: `teaching/index.qmd`. Each course is a card (`::: {.wm-course}`). To add a semester, copy
a card under a new `### Semester Year` heading. Syllabi and handouts go in `files/teaching/`.
The page has hidden `TO FILL IN` pieces for your approach, course descriptions, materials and
the Directed Reading Program.

## Statements, the AI page, and the private branch

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
page (`statements/index.qmd`) are drafts too, and work the same way. To add an AI project, copy
`ai/projects/_template.qmd`.

## Animations and the isogeny graphs

The interactive graph in the home page's side panel and on the Animations pages is one shared piece:

| Part | File |
|---|---|
| Markup: controls, the info panel, the legend, the glossary behind the "?" buttons | `_partials/isogeny-widget.qmd` |
| Drawing and animation | `assets/js/isogeny-hero.js`, `assets/js/spectral.js` |
| The random walk's transition matrix (no drawing) | `assets/js/walk.js`, tested by `_tests/walk.test.mjs` |
| Data | `assets/data/isogeny/`, written by `_scripts/isogeny_graphs.py` |
| The Animations pages | `animations/walks.qmd`, `balls.qmd`, `spectral.qmd`, `about.qmd` (definitions), `index.qmd` |

**Change the description of an animation.** The box at the top of the picture has one block per
entry of the animation menu, in `_partials/isogeny-widget.qmd` under `THE INFO PANEL`. Edit the
block; the page shows the block of the animation that is running. On an Animations page the
page title follows the menu too.

**Change a definition behind a "?" button.** Edit its `<dd>` in the glossary at the bottom of
`_partials/isogeny-widget.qmd`. The longer statements are on `animations/about.qmd`.

**Text under the picture on an Animations page** is ordinary text in that page's `.qmd` file.
Keep text about one particular animation out of the top of the page: the top belongs to the
info panel, which changes with the menu.

**Change the primes on offer.** Run

```
python3 _scripts/isogeny_graphs.py 101 431 1009      # your list of primes
python3 _scripts/isogeny_graphs.py                    # the default list
```

and commit the files it writes: one `p<prime>-<graph>.json` per prime and graph, plus
`index.json`. It also rewrites `spectra.json` (for the spectral law: histograms of the eigenvalues of T_2, T_3, T_5 and of the sums T_2 + T_3 and T_2 + T_3 + T_5, the joint histogram of (λ_2, λ_3), correlations and second moments; `spectra_row` in the script);
add `--skip-spectra` to leave that file alone, or `--spectra-only` to remake only it. The script
needs Python 3 and `numpy`, and it checks its own output on every prime (class number formulas,
modular polynomials, mass formulas, the Brandt relation, the eigenvalue bounds) and stops if a
check fails. The graphs it computes, and the largest prime for each, are listed at the top of the
script.

The default list has two parts, both at the top of the script. `LEVEL_PRIMES` are the 21 primes
that get every graph, with level structure up to the caps. `SAMPLE_PRIMES` are 200 more primes
for the graphs without level structure only: `sample_primes()` spreads them roughly
geometrically from 101 to 50000 and cycles through p mod 12, and it is deterministic, so a
full run reproduces the same files. The full default run takes about two hours on one core
(the sample's largest primes take four to six minutes each). The run deletes all graph files
first, so to add primes without remaking everything, call `compute(p, level5=False)` and
`write_prime` from Python for the new primes, then `write_index`. Graphs above 2500 vertices
(p beyond about 30000) get a lighter layout in the browser (`HUGE` in
`assets/js/isogeny-hero.js`).

**Directed edges.** Edges are counted by kernels, so the adjacency matrix is not symmetric next
to a vertex with automorphisms other than ±1. A link in the data is `[s, t, ell, m]` (m edges
each way) or `[s, t, ell, m, back]` where the two counts differ, and the picture draws those
edges with arrowheads. `animations/about.qmd#directed` states the rule.

**Pause and "until".** Pause freezes the frame loop and the layout, and Play resumes where it stopped (`setPlaying` in `assets/js/isogeny-hero.js`). The "until" menu (`.wm-graph-until-sel`) switches the walk and the bloomings to stopping at the F_p-subgraph; `animations/about.qmd#until` states the rule.

**The random walk.** The walker's steps are drawn from a transition matrix `P` built in
`assets/js/walk.js` from the adjacency matrix `A`. On the full graphs, with or without level structure, it is
the Petersson walk, `P[i][j] = A[i][j] w_j / s_i` with `w_j = |Aut(E_j)| / 2` (with level
structure, the automorphisms that preserve it) and
`s_i = sum_k A[i][k] w_k` (`peterssonWalk`); on the other views (F_p-subgraph, quotient) it is the simple random
walk (`simpleWalk`). `walkFor` chooses between them: it uses the Petersson walk when the graph's
data file has an `aut` field, which every full graph has, with or without level structure
(`|Aut(E, alpha)|` for the vertices whose order is not the generic one, `autDefault`, 2 unless stated). `animations/about.qmd#directed` states the definition.
After changing that file, run the tests:

```
node --test _tests/*.test.mjs
```

**The static figures** on the Animations pages are made by `python3 _scripts/spine_diagrams.py`
(needs `pip install networkx`) and `python3 _scripts/ball_diagrams.py`.

### The Limit laws page

`animations/limits.qmd` has its own widget, separate from the isogeny graph:

| Part | File |
|---|---|
| Markup: the menu, the info panel (one block per picture) | `_partials/limits-widget.qmd` |
| Drawing and animation | `assets/js/limits.js` |
| Tests (the laws, the distances, the data files) | `_tests/limits.test.mjs` |
| Data | `assets/data/limits/`, written by `_scripts/limit_laws.py` |
| Statements and proofs | `animations/about.qmd#limits` |

The "ℓ → ∞" entry of the spectral law's ℓ menu draws the first picture of that page inside
the isogeny widget (`createKmSt` in `assets/js/limits.js`, called from `assets/js/spectral.js`).

**Remake the data.** `pip install cypari2 numpy`, then

```
python3 _scripts/limit_laws.py                      # everything, about an hour on one core
python3 _scripts/limit_laws.py planar diagonal      # only some files
PARI_STACK=6000000000 python3 _scripts/limit_laws.py planar   # planar's largest levels need ~6 GB
```

The lists of levels, primes and weights are at the top of the script, and it stops if one of
its checks fails (they are listed at the top of the script). The first picture, Kesten–McKay to
Sato–Tate, needs no data.

## Look: colours, fonts, layout

- Colours: `assets/theme-light.scss` and `assets/theme-dark.scss` (the `--wm-*` values near the
  bottom of each). Change a colour in both files.
- Layout, cards, buttons: `assets/styles.css`, which uses those `--wm-*` values so that light
  and dark mode both work.
- Menu font: `assets/head.html`.

## Where everything is

| Folder / file | What it is |
|---|---|
| `index.qmd` | Home page |
| `research/index.qmd`, `research/talks.yml` | Research page: papers, and the list of talks |
| `writing/`, `notes/`, `blog/` | Writing pages; notes in `notes/articles/`, posts in `blog/posts/` |
| `code/` | Code pages, driven by `code/projects.yml` |
| `animations/` | Animations pages |
| `teaching/index.qmd` | Teaching page |
| `cv/` | CV page, the LaTeX source of the CV (`cv/cv.tex`), and `cv/README.md` |
| `cv_william_mahaney.pdf` | The CV PDF the site serves when there is no `cv/cv.tex` |
| `statements/`, `ai/` | Statements and the AI page (drafts) |
| `ai-disclosure/index.qmd` | How AI was used in building the site |
| `files/papers/`, `files/talks/`, `files/teaching/`, `files/statements/` | Uploads: PDFs. Each folder has an `_UPLOAD_HERE.md` note |
| `images/` | Uploads: profile photo, figures, favicon, social card |
| `_partials/` | Shared pieces: the isogeny graph, the job-search line |
| `assets/` | Styles, scripts, graph data, the two small Lua filters |
| `_scripts/` | Python scripts that compute the graph data and figures |
| `_tests/` | Tests of the random-walk code |
| `_quarto.yml`, `_quarto-public.yml`, `_quarto-private.yml` | Site settings and the two menus |
| `.github/workflows/publish.yml` | Builds and publishes the site |
| `.vscode/settings.json` | Makes VS Code's LaTeX extension write its files to a `build/` folder |
