# Decision: Weeks 9-12 extensions (task W01)

**Status: DRAFT, prepared by Claude. A human must choose one option below and sign it.**

## Where things stand (from tracker/esr1.log, 2026-10-20)

Already rebuilt as active code in `scripts/` on the 788-chain run, each logged DONE:

- ESR2: X01 (entries), X02 (RMSD matrix), X03 (clusters), X04 (PCA), X05 (basin sizes).
- ChEMBL: X13, X14 (earlier), X15 (potency by cluster), X16 (tests).

Also done since (2026-10-21): pockets X07, X08 (P2Rank 2.5.1 on all 788 chains), X09 (pocket table) and X12 (pocket PCA).

Not done yet:

- Pockets: X10 and X11. They need the Q10 ligand labels.
- Statistics: T01-T03. These wait on L03 and X11, which wait on C08 (choice of K) and Q10.

## Options

1. **Accept the archived outputs.** The Weeks 9-12 results stay as delivered. The README says that X07-X16 and T01-T03 are archived results and not reproduced by active code.
2. **Implement as active code.** Finish X07-X12 and T01-T03 in `scripts/`, on the same keep-all data. This needs the Q10 labels and a decision on K for the parts that use ligand class or chemistry clusters.

## Recommendation (not a decision)

Option 2, implement. Most of the extensions already run as active code. P2Rank runs here. The parts still open are blocked only by the human tasks Q10 and C08.

## Decision

Chosen option: ______  Name / date: ______
