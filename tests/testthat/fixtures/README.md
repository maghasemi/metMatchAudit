# Test fixtures

`swisslipids_mini.tsv` has the same header as the SwissLipids download (`lipids.tsv`) and 21 made-up entries. Entry IDs start at SLM:000000901 so they cannot be mistaken for real ones.

The InChIKeys attached to the isomeric entries are real keys from the MetLinkR example files, so tests can use real row names:

| key | row in the example data |
|---|---|
| BLZVZPYMHLXLHG-RQOIEFAZSA-N | oleoyl-linoleoyl-glycerol (18:1/18:2) [1], LEOCC |
| GREDRAMJRDQWEJ-LTEAFHAISA-N | oleoyl-linoleoyl-glycerol (18:1/18:2) [2], LEOCC (stored here with protonation letter O, to test tier 2) |
| RZRNAYUHWVFMIP-KTKRTIGZSA-N | 1-oleoylglycerol (18:1), LEOCC |
| LKQLRGMMMAHREN-YJFXYUILSA-N | stearoyl sphingomyelin (d18:1/18:0), LEOCC, and SM 36:1;O2, Broad |
| VJLLLMIZEJJZTE-NNTBDIJYSA-N | glycosyl-N-palmitoyl-sphingosine (d18:1/16:0), LEOCC |

The hierarchy copies situations seen in the feasibility run: an entry with no Level and no parent that needs the connectivity fallback (1-oleoylglycerol), a connectivity block whose entries point to two different Species (no fallback allowed), a glycosylceramide under a HexCer Species, and a loop in the parent links.
