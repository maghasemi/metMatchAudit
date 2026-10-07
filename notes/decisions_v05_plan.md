
## 2026-10-07: v0.5 plan, rows checked against their RefMet entry, written before any code

### Why the first v0.5 plan was dropped

The plan to give keyless rows a key from their RefMet name was dropped before it was logged. MetLinkR's `Standardized name` is a RefMet name (`workflow.R` at commit a6aa145 queries RefMet name conversion first, then RaMP synonyms, then RefMet again). A key taken from RefMet would agree with MetLinkR's grouping by construction. The downloaded table also has no synonyms, so exact name lookup would have found very little.

### What v0.5 checks instead

MetLinkR claims, for each row, "this row is RefMet entry X". v0.5 tests that claim against the row's own evidence and reports conflicts as row flags. It does not change any pair result. A conflict on one row says nothing yet about the other row of a pair, so the pair rules stay as they are.

### Source

`data-raw/refmet/refmet_all.txt`, from `https://www.metabolomicsworkbench.org/rest/refmet/all/xx/all/txt`, downloaded 2026-10-07. md5 `63079d6af6aff5fec9ddc9d3f952afa2`, 208170 entries, 35618 with a valid InChIKey, no name used twice. A different file stops the run until the change is logged. Group names missing from the table are reported and left without a check.

### The check, fixed before the run

For each row whose group name is a RefMet entry, the first case that applies decides:

1. **Entry is a lipid at species level** (its name parses with rgoslin at species level, such as `CAR 5:0` or `MG 18:2`). Its InChIKey, when present, stands for one member of a sum composition and is not used.
   - If the row has a lipid profile, compare Species. Headgroup naming differences (HexCer against GlcCer) are not a conflict, as in R04h.
   - Else, if the row has a key, compare molecular formulas: the formula of the row's parent compound against the entry's formula.
2. **Entry is a lipid at molecular level or finer and the row has a lipid profile:** compare Species.
3. **Entry has an InChIKey and the row has a key:** compare block1 of the two parent compounds.
4. Otherwise there is nothing to compare.

Parent compounds come from PubChem, through the client and cache built for identifier resolution: the row's key and the entry's PubChem CID are each mapped to their parent compound (`cids_type=parent`), which removes salts and counterions. If several parents come back and they disagree on block1, or the lookup fails, the row is reported as unchecked, not flagged.

When the comparison finds a difference:

- If the row's name and the entry's name are the same after `normalize_name()`, the row gets **C07 `refmet_identifier_conflict`**. The names agree and the structures do not, so one of the identifiers is wrong or the two describe different forms (ring and open chain).
- Otherwise the row gets **C06 `refmet_assignment_conflict`**: the row's own evidence says it is not the compound MetLinkR assigned it to.

Group names, alternative names and internal standards are skipped, as in identifier resolution.

### Known limits, stated now

- block1 ignores stereo, so a row assigned to the wrong stereoisomer is not caught.
- Rows with no key and no lipid profile (most of the Biocrates, LC/MS and Nightingale files in COMETS) are not checked.
- Ring and open-chain forms with different names are flagged C06 even when they are the same sugar.

### Acceptance, on the development data

- Every C06 flag is reviewed by hand and labeled real or not real. Real means the row is not the entry's compound, or for a species-level entry, not within its sum composition. At least 80% must be real.
- These six rows, which the manual review placed in C or E groups, must get C06: `andro steroid monosulfate C19H28O6S (1)*`, `(15:0)-anacardic acid`, `docosapentaenoate (n3 DPA; 22:5n3)`, `N-acetyl-1-methylhistidine*`, `octadecadienedioate (C18:2-DC)*`, `pregnenediol sulfate (C21H34O5S)*`.
- Pair results and the development evaluation stay identical to v0.4.2, since no pair rule changes. All tests pass.
- COMETS is run and reported as after adaptation, with the number of C06 and C07 flags and the rows behind them. No expert label is used to tune anything.

If an acceptance check fails, it is logged as failed. The check is not tuned on the same data to make it pass.
