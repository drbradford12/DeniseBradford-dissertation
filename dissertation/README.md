# Dissertation draft (UNL format)

Quarto book that builds the advisor draft with the UNL `nuthesis` class.

Build: `quarto render --to pdf` (needs R with knitr/rmarkdown and a TeX
distribution with XeLaTeX). The PDF lands in `_book/`.

- `01-introduction.qmd` - Chapter 1 (adapted from the introduction of comprehensive_exam/second_attempt.qmd)
- `02-literature-review.qmd` - Chapter 2 (from Literature Review/literature_review_numerical_ties.qmd)
- `03-method.qmd` - Chapter 3 (numerical_ties_method.qmd, corrected to match pcp_spread_ties.R,
  plus Section 4 material from the November 9, 2026 proposal)
- `A-repository-guide.qmd`, `B-source-code.qmd` - appendices
- `references.bib` - the two .bib files merged (Literature Review version wins on shared keys)
- `partials/` - UNL title page, abstract, and front matter; `_quarto.yml` holds title-page fields

Figures (run from this folder):
- `Rscript scripts/make_ch2_figures.R figures` - beeswarm and mtcars figures for Chapter 2
- `Rscript scripts/make_method_figures.R code/pcp_spread_ties.R figures` - iris, penguins, and flea beetle figures for Chapter 3
- `figures/fig_worked_example.png` is the proposal's hand-drawn diagram (spacing_polyline.png), copied as is

Also in this archive: `comprehensive_exam/revisions/` - revision logs for the proposal; drop it into the existing `comprehensive_exam/` folder.
