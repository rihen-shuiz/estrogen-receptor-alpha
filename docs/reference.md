# Reference structure

Chosen reference: **7NFB, chain B, model 1** (label chain B, author chain B). Chain key `7NFB|B|1`.

Selection rule, from the original project's `choose_reference.py`, applied to `results/lbd_chains.csv` and `results/entries.csv`:

1. ESR1 chain, eligible under the 180 LBD CA and 60 core CA screen, and entry QC `pass`.
2. LBD coverage at least 95 %, core helices at 100 %, H12 at 100 % of expected positions.
3. Sort by resolution (best first), then mapped LBD CA (most first), then declared differences (fewest first), then chain key.

Result: 40 ESR1 chains pass step 2. 7NFB|B|1 ranks first by resolution (1.33 A). This matches the choice noted for the original project's clean version (task card S01).

The original project's older historical reference is 3Q95|B|1. The original notes that new matrices will be recomputed with the chosen chain instead.

Evidence for 7NFB|B|1:

- Resolution 1.33 A, R-free 0.178, mean B 26.75 A^2, all three QC flags `pass`.
- LBD: 244 of 246 expected positions have CA. The two missing residues are 549 and 550, at the C-terminus.
- Core helices: 84 of 84 CA. H12: 10 of 10 CA.
- Declared differences from the canonical sequence: 28, all annotated as engineered mutations. None are unexplained or conflicting.

Limitations:

- This is an engineered variant, not wild-type ESR1. The primary citation describes it as a PROSS-designed, extensively mutated ERalpha variant with enhanced thermal stability. The 28 substitutions include positions inside the LBD (for example 315, 381, 417, 530, 537). Any conformational result from this reference must be read with that in mind.
- The reference fixes one crystal form and one chain. It is a coordinate frame, not a claim about the native receptor's conformation.
- Resolution alone does not establish stability, pharmacological action, or biological superiority.
- Entry ligand codes are CL, EDO, GEN, GOL, NA. This list is for the whole entry. Contact-level ligand assignment for chain B was not computed here.
