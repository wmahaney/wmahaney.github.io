# The CV

The CV is written in LaTeX. This folder holds its source, and the website shows the PDF
built from it, so there is one copy of the CV.

## Where the files go

| File | What it is |
|---|---|
| `cv/cv.tex` | The main file. **Upload it here under exactly this name.** |
| `cv/*.cls`, `cv/*.sty`, `cv/*.bib` | Any class, style or bibliography file `cv.tex` needs, next to it. |
| `cv/.latexmkrc` | Optional. Only needed to compile with something other than pdfLaTeX (see below). |
| `cv/index.qmd` | The CV page of the website. It shows the PDF; leave it alone. |
| `cv_william_mahaney.pdf` (top folder) | The PDF the site serves. Built from `cv/cv.tex` on every push. |

## Updating the CV

1. Edit `cv/cv.tex`.
2. Commit and push.

GitHub then compiles the CV and rebuilds the site. The CV page, the "Curriculum vitae"
button on the home page and every existing link to `cv_william_mahaney.pdf` get the new PDF.

## What happens on push

The workflow `.github/workflows/publish.yml` does this before it builds the site:

- If `cv/cv.tex` exists, it compiles it with `latexmk` in TeX Live 2025 and writes the
  result to `cv_william_mahaney.pdf` for that build.
- If `cv/cv.tex` does not exist, the step is skipped and the `cv_william_mahaney.pdf`
  already in the repository is used. So nothing changes until you upload `cv.tex`.
- If the compilation fails, the build fails and the site is not published; the live site
  keeps the previous CV. The LaTeX log is attached to the run: Actions → the failed run →
  Artifacts → **cv-latex-log**.

The compiled PDF is not committed back to the repository. The copy of
`cv_william_mahaney.pdf` in the repository is only the fallback for when there is no
`cv/cv.tex`; once `cv.tex` is in place you can leave that file as it is.

## Building it on your own computer

```sh
cd cv
latexmk -pdf -outdir=build cv.tex      # writes cv/build/cv.pdf
```

In VS Code with the LaTeX Workshop extension, saving `cv.tex` does the same: the
repository's `.vscode/settings.json` sends the PDF and all intermediate files to a
`build/` folder next to the `.tex` file. `build/` folders are ignored by git.

## Another engine

pdfLaTeX is the default. For LuaLaTeX or XeLaTeX, create `cv/.latexmkrc` with one line:

```perl
$pdf_mode = 4;   # 4 = LuaLaTeX, 5 = XeLaTeX
```

Both the workflow and a local `latexmk` read this file.
