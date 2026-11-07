# Results

Output of the scripts in `scripts/`. Each folder or file is listed with the script that writes it. Files ending in `_log.txt` are the terminal output of a run, saved with `tee`. Results that name rows of the private test file are not here; they stay in `data-raw/own/`, which git ignores.

| path | written by |
|---|---|
| `rows.csv`, `rows_summary.csv` | `scripts/02_manual_review/01_build_rows.R` |
| `inchikey_groups.csv`, `pattern_counts.csv`, `flagged_rows.csv` | `scripts/02_manual_review/02_inchikey_groups.R` |
| `lipid_feasibility/` | `scripts/03_lipid_feasibility/01_lipid_feasibility.R` |
| `audit_v0.1/`, `audit_v0.1_sealed.rds`, `pubchem_cache.rds` | `scripts/04_v0.1/01_audit_example.R` |
| `evaluation_v0.1/` | `scripts/04_v0.1/02_evaluate.R` |
| `v02_feasibility/` | `scripts/05_v0.2/01_feasibility.R` |
| `v02_check/` | `scripts/05_v0.2/02_check_proposals.R` |
| `v03_feasibility/` | `scripts/06_v0.3/02_feasibility.R` |
| `v03_rdkit/` | `scripts/06_v0.3/04_rdkit_feasibility.R` |
| `id_resolution/` | `scripts/07_v0.4/02_id_resolution_feasibility.R` |
| `id_resolution_retest/` | `scripts/07_v0.4/03_id_resolution_retest.R` |
| `audit_v0.4.1_before.rds` | no script; the development audit kept by hand before v0.4.2 |
| `dev_v0.4.2/` | the two v0.1 scripts rerun with v0.4.2, moved here by hand (decisions.md, 2026-10-07) |
| `comets_test/` | `scripts/08_comets/02_blind_test.R` |
| `comets_v0.4.2/` | `scripts/08_comets/03_after_v0.4.2.R` |
| `comets_recount/`, `comets_recount_log.txt` | `scripts/08_comets/04_recount_all_versions.R` |
| `refmet_inventory/` | `scripts/09_v0.5/01_refmet_inventory.R` |
| `v05_check/` | `scripts/71_v05_check.R` on branch `shelved/refmet-check` |
| `biocrates_feasibility/` | `scripts/10_v0.6/02_biocrates_feasibility.R` |
| `lmsd_feasibility/` | `scripts/10_v0.6/03_lmsd_feasibility.R` |
| `v06_check/` | `scripts/10_v0.6/04_acceptance.R` |
| `blind2/` | `scripts/11_blind_test_2/04_test.R` |
| `blind2_v041_diagnostic/` | `scripts/11_blind_test_2/05_v0.4.1_diagnostic.R` |
| `v07_feasibility/`, `v07_feasibility_log.txt` | `scripts/12_v0.7/01_feasibility.R` |
| `v07_check/`, `v07_check_log.txt` | `scripts/12_v0.7/02_acceptance.R` |
| `blind3/` | `scripts/13_blind_test_3/04_test.R` |
| `v08_feasibility/`, `v08_feasibility_log.txt` | `scripts/14_v0.8/01_feasibility.R` |
| `v09_feasibility/` | `scripts/15_v0.9/01_feasibility.R` |
| `v09_check/`, `v09_check_log.txt` | `scripts/15_v0.9/02_acceptance.R` |
