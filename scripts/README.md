# Scripts

One folder per stage, in the order the work was done. Within a folder, the number gives the order to run the scripts. Every script runs from the project root and says what it needs in its header. The plan and the result of each stage are in `notes/decisions.md` and `notes/findings.md`.

| folder | stage |
|---|---|
| `01_setup/` | Download the MetLinkR example data and build the slim SwissLipids table. |
| `02_manual_review/` | Rows, InChIKey groups and the sheet for the first manual review. |
| `03_lipid_feasibility/` | Lipid evidence from SwissLipids and rgoslin, checked before any package code. |
| `04_v0.1/` | First audit of the example data and its evaluation against the manual review. |
| `05_v0.2/` | Lipid match proposals across vendors. |
| `06_v0.3/` | Stereo comparison from full structures (RDKit). |
| `07_v0.4/` | Identifier resolution through PubChem. |
| `08_comets/` | External test on five COMETS datasets, with the MetLinkR run, the rerun after v0.4.2 and a recount of all versions. |
| `09_v0.5/` | RefMet inventory for the shelved v0.5 check. |
| `10_v0.6/` | Biocrates p180 table, LIPID MAPS feasibility and v0.6 acceptance. |
| `11_blind_test_2/` | Second blind test, with the MetLinkR run and a diagnostic with v0.4.1. |
| `12_v0.7/` | v0.7 feasibility and acceptance. |
| `13_blind_test_3/` | Third blind test, with the MetLinkR run. |
| `14_v0.8/` | v0.8 feasibility (the change was dropped). |
| `15_v0.9/` | v0.9 feasibility, acceptance and the evidence for the hand check. |
| `16_v0.10/` | v0.10 feasibility (the change was dropped). |
| `17_v0.9.1/` | v0.9.1 acceptance. |
| `18_v0.11/` | v0.11 feasibility and acceptance. |
| `19_blind_test_4/` | Fourth blind test, on three Metabolomics Workbench files. |

The three `*_run_metlinkr.R` scripts are not run from this project. MetLinkR, RaMP and Java are kept in a separate R project so that they stay out of this package's `renv.lock`. The scripts are copied here unchanged so that every step of the tests can be read in one place. Their headers give the paths they used.

## Old names

Until 2026-10-07 the scripts sat in one folder with two- or three-digit prefixes. `notes/decisions.md`, `notes/findings.md`, the label notes and the logs in `results/` were written at the time and still use the old names. This table maps them.

| old name | now |
|---|---|
| `scripts/00_get_data.sh` | `scripts/01_setup/01_get_data.sh` |
| `scripts/01_build_rows.R` | `scripts/02_manual_review/01_build_rows.R` |
| `scripts/02_inchikey_groups.R` | `scripts/02_manual_review/02_inchikey_groups.R` |
| `scripts/03_make_review_sheet.R` | `scripts/02_manual_review/03_make_review_sheet.R` |
| `scripts/10_swisslipids_slim.R` | `scripts/01_setup/02_swisslipids_slim.R` |
| `scripts/11_lipid_feasibility.R` | `scripts/03_lipid_feasibility/01_lipid_feasibility.R` |
| `scripts/20_audit_example.R` | `scripts/04_v0.1/01_audit_example.R` |
| `scripts/21_evaluate_v0.1.R` | `scripts/04_v0.1/02_evaluate.R` |
| `scripts/30_v02_feasibility.R` | `scripts/05_v0.2/01_feasibility.R` |
| `scripts/31_v02_check_proposals.R` | `scripts/05_v0.2/02_check_proposals.R` |
| `scripts/40_rdkit_reference.csv` | `scripts/06_v0.3/01_rdkit_reference.csv` |
| `scripts/40_rdkit_reference.py` | `scripts/06_v0.3/01_rdkit_reference.py` |
| `scripts/40_v03_feasibility.R` | `scripts/06_v0.3/02_feasibility.R` |
| `scripts/41_inchi_stereo_functions.R` | `scripts/06_v0.3/03_inchi_stereo_functions.R` |
| `scripts/41_v03_rdkit_feasibility.R` | `scripts/06_v0.3/04_rdkit_feasibility.R` |
| `scripts/50_id_functions.R` | `scripts/07_v0.4/01_id_functions.R` |
| `scripts/50_id_resolution_feasibility.R` | `scripts/07_v0.4/02_id_resolution_feasibility.R` |
| `scripts/51_id_resolution_retest.R` | `scripts/07_v0.4/03_id_resolution_retest.R` |
| `scripts/60_comets_test.R` | `scripts/08_comets/02_blind_test.R` |
| `scripts/61_comets_after_v042.R` | `scripts/08_comets/03_after_v0.4.2.R` |
| `scripts/70_refmet_inventory.R` | `scripts/09_v0.5/01_refmet_inventory.R` |
| `scripts/80_biocrates_table.R` | `scripts/10_v0.6/01_biocrates_table.R` |
| `scripts/81_biocrates_feasibility.R` | `scripts/10_v0.6/02_biocrates_feasibility.R` |
| `scripts/82_lmsd_feasibility.R` | `scripts/10_v0.6/03_lmsd_feasibility.R` |
| `scripts/83_v06_check.R` | `scripts/10_v0.6/04_acceptance.R` |
| `scripts/91_blind2_prepare.R` | `scripts/11_blind_test_2/01_prepare.R` |
| `scripts/92_blind2_sample.R` | `scripts/11_blind_test_2/03_sample.R` |
| `scripts/93_blind2_test.R` | `scripts/11_blind_test_2/04_test.R` |
| `scripts/94_blind2_v041_diagnostic.R` | `scripts/11_blind_test_2/05_v0.4.1_diagnostic.R` |
| `scripts/95_v07_feasibility.R` | `scripts/12_v0.7/01_feasibility.R` |
| `scripts/96_v07_check.R` | `scripts/12_v0.7/02_acceptance.R` |
| `scripts/97_comets_recount.R` | `scripts/08_comets/04_recount_all_versions.R` |
| `scripts/101_blind3_prepare.R` | `scripts/13_blind_test_3/01_prepare.R` |
| `scripts/102_blind3_sample.R` | `scripts/13_blind_test_3/03_sample.R` |
| `scripts/103_blind3_test.R` | `scripts/13_blind_test_3/04_test.R` |
| `scripts/104_v08_feasibility.R` | `scripts/14_v0.8/01_feasibility.R` |
| `scripts/71_v05_check.R` | not on main; kept on branch `shelved/refmet-check` |
| `scripts/105_v08_check.R` | never written; planned for v0.8, which was dropped at feasibility |
| (kept in the separate MetLinkR project) | `scripts/08_comets/01_run_metlinkr.R`, `scripts/11_blind_test_2/02_run_metlinkr.R`, `scripts/13_blind_test_3/02_run_metlinkr.R` |
