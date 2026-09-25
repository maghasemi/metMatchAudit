# metMatchAudit

Work in progress. This repository stays private until the first benchmark is written up.

MetLinkR (Patt et al., J. Proteome Res. 2025, 24, 2403-2407) links metabolites across datasets by giving each row a RefMet standardized name. The authors compared its output with manual curation from the COMETS consortium. They found that it sometimes groups compounds reported at different structural resolution, groups stereoisomers such as bilirubin (E,E) and bilirubin (Z,Z), and in a few cases links compounds that are not isomers at all, such as C14:1 and C14:2 acylcarnitines. Their advice is to inspect the output by hand.

This project checks those groupings against structural evidence and reports how two rows relate to each other, instead of a plain match or no match.

## Current phase

Phase 1a asks a narrow question. When MetLinkR gives several rows the same standardized name, do those rows also share the same InChIKey? The InChIKeys come from MetLinkR's own example files, so this step needs no API calls.

## Running phase 1a

Install R packages once (see "Setup" below), then run these from the project root:

```
bash scripts/00_get_data.sh
Rscript scripts/01_build_rows.R
Rscript scripts/02_inchikey_groups.R
Rscript scripts/03_make_review_sheet.R
```

`00_get_data.sh` clones MetLinkR at a fixed commit into `data-raw/`, which git ignores. The R scripts write to `results/`. The last script creates `notes/review.csv` for manual review and will not overwrite it once it exists.

## Setup

In R, from the project root:

```r
install.packages("renv")
renv::init()
install.packages(c("readr", "dplyr", "stringr"))
renv::snapshot()
```

Commit `renv.lock` afterwards. It records the package versions used for every result in this repository.

## Layout

```
scripts/   numbered scripts, run in order
results/   script output
notes/     decisions.md (what was decided and why), findings.md, review.csv
data-raw/  downloaded data, not tracked by git
```
