# wmahaney.github.io

Personal academic website of William E. Mahaney, built with [Quarto](https://quarto.org) and
published with GitHub Pages. Live site: <https://wmahaney.github.io>

## Where things go

| Folder / file | What it is | What to do |
|---|---|---|
| `index.qmd` | Home page | Edit text, links |
| `research/index.qmd` | Research page: papers and talks | Add papers (instructions inside the file) |
| `notes/articles/<name>/index.qmd` | **Expository Notes**, one folder per note | Copy `notes/_template.qmd` |
| `blog/posts/<name>/index.qmd` | Blog posts, one folder per post | Copy `blog/_template.qmd` |
| `teaching/index.qmd` | Teaching page | Edit text |
| `cv/index.qmd` | CV page | Edit text |
| `cv_william_mahaney.pdf` | The CV PDF (top-level, so the old link keeps working) | Replace with a file of the same name |
| `images/profile/` | **Upload your profile photo here** as `profile.jpg` | Upload |
| `images/notes/`, `images/blog/` | **Upload figures** for notes and posts | Upload |
| `files/papers/` | **Upload paper PDFs** | Upload |
| `files/talks/` | **Upload slides** | Upload |
| `files/teaching/` | **Upload syllabi and handouts** | Upload |
| `files/cv/` | **Upload research / teaching / personal statements** | Upload |
| `assets/` | Colors, fonts, and page styling | Rarely needed |
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

## Changing the look

Colors: `assets/theme-light.scss` and `assets/theme-dark.scss` (the `--wm-*` values near the bottom).
Layout and cards: `assets/styles.css`. Menu font: `assets/head.html`.
