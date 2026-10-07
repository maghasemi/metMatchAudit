
## 2026-10-07: v0.5 feasibility plan, name to RefMet, written before any code

### Question

739 of 1429 COMETS rows have no structure key, because three files carry names only. Can a row's own name give it a key that is accurate enough to audit with?

### Source

RefMet from the Metabolomics Workbench. The REST service returns the whole table at `https://www.metabolomicsworkbench.org/rest/refmet/all/xx/all/txt`, with name, regno, pubchem_cid, inchi_key, exactmass, formula, sys_name, main_class, sub_class and synonyms. The table is downloaded once, its md5, download date and row count are logged here, and lookups run offline against that file, the same way SwissLipids is used. The `match` service is not used: it changes on the server without a version, so results could not be repeated.

### Independence

MetLinkR also uses RefMet. If MetLinkR built its standardized names from the same RefMet lookup, a key from RefMet would agree with MetLinkR's grouping by construction and would audit nothing. Before the feasibility run, read the MetLinkR code at commit a6aa145 and write down how it uses RefMet. The `Standardized name` column is never used as input to the lookup.

### Lookup, fixed before the run

1. Normalize the row's name with `normalize_name()`.
2. Look it up against RefMet names and synonyms, normalized the same way. Exact matches only.
3. Collect the InChIKeys of all hits. If they share one block1, the row gets that key, and the stereo block only if all hits agree on it (otherwise the key is kept at block1 and marked stereo unknown). If block1 differs between hits, the row gets no key and is counted as ambiguous.
4. Rows skipped by `resolve_identifiers()` are skipped here too: group names, alternative names, Species and molecular lipid names, internal standards.

### Data

- Accuracy set: every row in the development data (LEOCC and Broad) and in COMETS that already has a key, from the file or from identifier resolution. COMETS expert labels are not read.
- Coverage set: the 739 COMETS rows without a key. Only the count of rows that get a key is reported, no comparison.

### Criteria

- Accuracy: among accuracy-set rows that get a RefMet key, block1 agrees with the existing key in at least 95%. Same bar as identifier resolution.
- Coverage: at least 30% of the 739 keyless rows get a key. Below that, the route is not worth building and v0.5 goes another way.
- Both must hold. Every block1 disagreement is listed and read.

If either fails, the criteria stay as written and the result is logged as a failed feasibility test.
