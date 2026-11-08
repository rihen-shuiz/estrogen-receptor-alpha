Not done / substituted: (1) no peer-reviewed source on the residue-level effect of **L536S** (the construct mutation of 243 of our 788 chains) could be opened - searched repeatedly, found only USPTO patents (not opened: not peer-reviewed, large scanned PDFs) and PDB entry metadata; the sources used for it are the RCSB 2QXS entry page and the Bruning et al. 2010 abstract (Europe PMC record), plus MD of L536Q (Pavlin 2018) as the nearest mechanistic evidence, all flagged as such in docs/sv/literature.md; (2) 8 sources could not be opened (Nat Commun version of the 2025 ternary-switch paper and its PMC copy, the same preprint's PDF, Shiau 1998 Cell, Gangloff 2001 JBC, the PMC3268042 review, the PubMed page for PMID 20924370, one Europe PMC XML endpoint) - each is listed with the reason and with what was used instead; the ternary-switch claims therefore come from the bioRxiv preprint, not the journal version; (3) WebFetch returns opened pages through a reader model that caps verbatim quotes at 125 characters, so claims are short quotes plus close paraphrase of the opened text, not long verbatim passages; (4) rule 1.9: claims framed by ligand action class were skipped, and where an opened source labelled a *conformation* "agonist/antagonist" the claim was reworded to the structural state (H12 docked against H3/H11, H12 in the coactivator groove) - no ligand action label is carried into the table; (5) no figure, no script and no computed number: the only use of the venv was counting construct-mutation values in the SV00 chain table for the "bears on" column.

# SV08 - Literature: is only H12 different between ERalpha LBD structures?

**Answer to the step question ("check literature if only H12 is significant"): no.** Opened sources report
reproducible, ligand- or construct-associated differences in helix 11, helix 7, helix 3, the H11-H12 loop and
the His524/Leu525 pocket rotamers, and unmodelled density in two loops outside H12 (331-340 and 462-463).

## Numbers

- 16 WebSearch queries, 28 WebFetch calls, **14 distinct sources opened**, **25 claims** recorded, **8 sources
  could not be opened** (reason and substitute given for each).
- Opened sources, by what was opened: 11 publisher/PMC full texts (Lee & Barron 2017; Nettles 2008; Li 2011;
  Tanenbaum 1998; Chen 2016; Djinovic-Carugo & Carugo 2015; Fanning 2016; Hosfield 2022; Mayne 2021;
  Pavlin 2018; Sendanayake 2025), 1 preprint full text (McDougal 2025), 1 abstract record (Bruning 2010 via
  Europe PMC), 2 PDB entry pages (2QXS, 2YAT). One further full text was opened and yielded no in-scope claim
  (Dai 2009).
- Claims per question: Q1 (other regions) 8, Q2 (missing H12 density) 7, Q3 (Y537S / L536S / D538G) 6,
  Q4 (flexibility near 333 and 467) 4.
- Largest literature numbers carried in: 73 hERalpha LBD structures grouped into three His524 classes
  (Lee & Barron 2017); 185 H12-docked ERalpha structures surveyed for Leu536, buried in >98 % of Y537S monomers
  (n = 170) vs ~18 % of wild-type monomers (13 structures, 28 monomers) (Mayne 2021); PDB-wide, 69 % of entries
  contain missing strings (Djinovic-Carugo & Carugo 2015).

## Verdict per question

1. **Other regions do differ.** Helix 11 shifts with ligand chemistry down to one methyl group, and the authors
   argue from density that this is not an artifact (Nettles 2008, doi 10.1038/nchembio.76); His524 takes three
   conformational classes across 73 structures, the "open" one for ligands longer than 13 A (Lee & Barron 2017,
   doi 10.1371/journal.pone.0169607); a bulky tagged ligand shortens H3 and H8 by one turn and disorders H7
   (Li 2011, doi 10.1021/jm200192y); Y537S remodels the H11-H12 loop and the helix-3 hydrogen-bond pattern
   (Nettles 2008). Associations named: ligand chemistry/size and construct mutation; **no opened source
   attributes a non-H12 difference to the coactivator peptide**.
2. **Missing H12 density** is reported in prose per structure and by omission from the deposited model
   (Fanning 2016, doi 10.7554/eLife.12792); mobility of a modelled H12 is reported as per-position main-chain
   B-factors normalised to the structure average, and crystal contacts at or near 538 affect H12 in one monomer
   but not the other, so the uncontaminated monomer is analysed (Hosfield 2022, doi 10.7554/eLife.72512).
   Absent H12 is also read as disorder of the undocked helix, with no apo wild-type structure existing at all
   (Mayne 2021, doi 10.1158/1541-7786.MCR-20-1017; McDougal 2025, doi 10.1101/2025.02.09.637324). PDB-wide,
   missing strings sit in loops with high flanking B-factors and are usually far from packing contacts, except
   at chain termini (Djinovic-Carugo & Carugo 2015, doi 10.1080/21690707.2015.1095697).
3. **Y537S** makes a new Ser537 hydrogen bond (D351 in the apo mutant, E380 in another series), remodels the
   H11-H12 loop, buries Leu536 in nearly all deposited Y537S monomers, and was introduced as a crystallisation
   aid (Fanning 2016; Hosfield 2022; Nettles 2008; Mayne 2021). **D538G** removes the Y537-N348 bond, moves the
   start of the H12 helix from 537 to 539, displaces the H11-H12 loop toward H11 and raises its flexibility in
   MD (Fanning 2016), relaxing the loop "spring force" (Mayne 2021). **L536S** is established only as an
   engineered conformation-driving construct mutation (RCSB 2QXS; Bruning 2010, doi 10.1038/nchembio.451); its
   mechanism is a gap (line 1).
4. **Flexibility near 333 and 467 is documented**: residues 331-340 (omega loop between H2 and H3) and 462-463
   are unmodelled with very weak density in PDB 2YAT (Li 2011); the 331-340 loop is distorted by H3's lever
   motion in MD/HDX (McDougal 2025); the H8-H9 loop belongs to the dimer interface (Tanenbaum 1998,
   doi 10.1073/pnas.95.11.5998).

## What this bears on in our results (full mapping in the table)

- **SV01a / SV01b / SV01c**: the completeness outcome is a REMARK-465-style missing string; the mobility outcome
  (normalised B over 528-547) is the same measure the ERalpha literature uses for H12 mobility; the named crystal
  contacts at 538 / K531 / D538-R436 are the same mechanism as our lattice terms, and the "missing strings are
  far from packing contacts except at termini" result supports the SV01b caveat; 443 Y537S + 243 L536S + 12 D538G
  of 788 chains means the construct-mutation block carries a well-documented H12 effect; there is no apo
  wild-type baseline anywhere in the PDB.
- **SV02**: the reported shift of the H12 helix start from 537 to 539 with D538G matches our cluster consensus
  (cl7 538-545, cl4 537-544) and the data-derived window 538-544; psi at 537 is the residue the literature tracks.
- **SV03**: H11 and H3 are independently reported as variable regions, supporting H10/11 length, H11 kink and
  H11-H3 angle as genuine (non-H12) discriminators; H7 and H8 lengths vary only with specific bulky ligands,
  which matches their not separating our clusters.
- **SV07**: both genuine non-H12 spike blocks (331-340, 464-468) coincide with loops that a deposited ERalpha
  structure could not model, and the two associations we found (ligand chemotype at 331-339, dimer partner at
  334-340/467-468) match the literature mechanisms (H3 lever; H8-H9 loop in the dimer interface). The peptide
  association at 466-468 has no literature support in the opened sources.
- **SV11**: the AF-2 groove set (E380, D538-M543) contains the residues of the Y537S S537-E380 bond; the
  dimer-interface set spans exactly the elements named as the interface in 1998; pocket-residue rotamer
  variability (His524) is finer than our site-level test.

## Files

- docs/sv/literature.md (25.0 kB) - the question | claim | source (opened) | bears-on table, the short answer per
  question, the list of sources that could not be opened, and the gaps.
- results/sv/SV08/SUMMARY.md (8.0 kB, this file)
- results/sv/SV08/env.txt (0.9 kB)

Nothing enters SV09 from here unless its source was opened; every row of the table names what was opened.
