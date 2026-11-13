# Protocol addendum 1 (2026-11-13), still before any experimental comparison

Reason: after freezing PROTOCOL.md, the min-fill treewidth of the 12 A interaction graphs was
computed (structure only, no experiment, no model output): 4PTI 3, 2LZT 4, 3RN3 7, 1UBQ 4,
1PGA 5, 1STN 8, 1A2P 9, 2RN2 7, 4HHB 17. For 4HHB (174 sites) exact inference is therefore
feasible (2^18 tables). PROTOCOL section 5 expected 4HHB to be "beyond exact inference"; that
expectation was wrong.

Changes:
1. The LARGE protein of section 5 becomes 1OEL (7 chains, all used, 1064 sites, greedy min-fill
   treewidth upper bound 31, i.e. tables of 2^32 doubles = 32 GB, beyond the 4 GB budget).
   1OEL is not in the experimental comparison set.
2. 4HHB stays in the comparison set and becomes a MID-SIZE system with an exact (float64
   variable-elimination) reference at pH 3, 5, 7, 9, 11 (not part of the trusted path) against which the
   certificates are compared.
3. For the 1064-site system the radius hierarchy {0,4,6,8} is run at pH 7.0 for all sites;
   the pKa bracket procedure (section 4) is run for all sites at R=8 only.
4. Exact validation uses 4PTI (20 sites) plus four 18-20 site subsets; as specified.
Nothing else changes.
