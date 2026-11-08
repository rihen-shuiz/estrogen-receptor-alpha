Not done / substituted: (1) the primary sources of the ERalpha pocket and E2 volumes could not be opened - Brzozowski et al. 1997 (Nature, paywalled; searched for an open repository copy, only metadata found) and Shiau et al. 1998 (Cell, HTTP 403 in SV08, not retried) - so the 450 A3 pocket and ~250 A3 ligand values in this file are quoted second-hand from Tanenbaum 1998 and Katzenellenbogen 2011, both opened, and flagged as such; no opened source gives a method (program, probe radius, grid) for any ERalpha volume. (2) Three Q1 claims rest on abstract records only (Eelen 2005 Mol Pharmacol, Eelen 2008 Chem Biol, Vanhooke 2004 Biochemistry: not in PMC or HTTP 403), one Q2 claim on a PDB entry page plus the abstract it carries (Pike 2001 Structure, via RCSB 1HJ1) and one on an entry page alone (Liang 2021, via RCSB 7MSA); 12 sources could not be opened and each is listed with the reason and the substitute. (3) WebFetch returns opened pages through a reader model that caps verbatim quotes at 125 characters, so claims are short quotes plus close paraphrase of the opened text. (4) Rule 1.9: unlike SV08 this step uses action-class claims, as the step says; wording keeps H12 order/position separate from degradation and cellular activity, and every source's evidence type (crystal, HDX-MS, limited proteolysis, MD, cell assay) is named in the claim. (5) No figure, no script and no computed number: the SV12-SV15 numbers in the "bears on" column are read from those steps' SUMMARY.md files. (6) Two bioRxiv preprints are used where the journal version could not be opened (Diennet 2023, no journal version found; McDougal 2025, Nat Commun version behind idp.nature.com as in SV08).

# SV08b - Literature for the class questions

**Answers to the three step questions: (1) no H12 position difference between VDR superagonists and agonists; (2) yes, SERD H12 destabilisation relative to SERMs is reported, mostly in solution, and its link to degradation is asserted rather than shown; (3) the published ERalpha volumes are 450 A3 (pocket) and ~250 A3 (E2), both without a stated method.**

## Numbers

- 20 WebSearch queries (1 extended), 31 WebFetch calls, **19 distinct sources opened**, **20 claims** recorded,
  **12 sources could not be opened** (reason and substitute given for each).
- Opened sources, by what was opened: 10 publisher/PMC/journal full texts, 2 preprint full texts, 4 abstract
  records (3 via Europe PMC, 1 PMC page carrying the abstract only), 3 PDB entry pages (3O1E, 7MSA, 1HJ1).
- Claims per question: Q1 (VDR superagonists vs agonists) 7, Q2 (SERD vs SERM H12, degradation) 8,
  Q3 (ERalpha pocket and E2 volumes) 5.
- Largest literature numbers carried in: Calpha RMSD 0.08-0.38 A between VDR superagonist and agonist complexes
  (Tocchini-Valentini 2001; Huet 2011); trypsin half-lives apo 2.2 / fulvestrant 2.3 / OP-1074 3.5 / 4-OHT 6.0 min
  and ERalpha remaining 28 % (fulvestrant) vs 49 % (OP-1074) of vehicle (Fanning 2018 Nat Commun); ERalpha pocket
  450 A3 with ligands 250 A3 (Katzenellenbogen 2011, quoting Brzozowski 1997 and Shiau 1998); VDR pocket occupancy
  61-63 % with cavity 729-776 A3 and ligands 461-474 A3 (Huet 2011).

## Verdict per question

1. **The professor's VDR point is correct as a statement about H12 position, and only about position.** Every
   opened VDR source reports the same for superagonists: "the position of helix H12 ... strictly maintained"
   across 1,25(OH)2D3, MC1288 and KH1060 (Calpha RMSD 0.08 and 0.14 A; Tocchini-Valentini 2001,
   doi 10.1073/pnas.091018698), "The conformation of the LBD is identical in each complex" for the superpotent 2MD
   (Vanhooke 2004, doi 10.1021/bi036056y), "The position and conformation of the activation helix H12 of VDR was
   strictly maintained" for the 14-epi TX522 (Eelen 2005, doi 10.1124/mol.104.008730), H12 "strictly conserved" and
   "Both complexes adopt the canonical, active conformation" for the gemini superagonists (Huet 2011,
   doi 10.1039/C1MD00059D), and a review stating that every VDR-LBD/ligand structure reported up to then adopted
   the active conformation (Yamamoto 2019, doi 10.1248/cpb.c19-00131). Superagonism is explained by the lifetime of
   the active complex, by stronger coactivator interaction, or by extra ligand contacts with H12 residues Val444 /
   Phe448 / Leu440 called "increased stabilization of helix H12" (Eelen 2008, doi 10.1016/j.chembiol.2008.08.008;
   Huet 2011). **Slide caution:** do not write "no H12 difference" without the word *position* - two of these
   papers do claim increased H12 stabilisation, inferred from contacts in static structures, with no H12 B-factors
   and no MD behind it. In VDR, H12 differences are reported for antagonists and partial agonists, not between
   agonist subtypes.
2. **SERD vs SERM H12 destabilisation is reported, and it is not a missing-density phenomenon.** HDX-MS: fulvestrant
   raises deuterium uptake near H12 while E2 protects it, and bazedoxifene stabilises H12 less than 4-OHT
   (Fanning 2018, doi 10.7554/eLife.37161). Limited proteolysis orders the ligands apo 2.2 < fulvestrant 2.3 <
   OP-1074 3.5 < 4-OHT 6.0 min, with H12 nonetheless modelled in the AF-2 cleft in every crystal and the effect
   visible as B-factors (Fanning 2018, doi 10.1038/s41467-018-04413-3). In long-chain SERD structures the H11-H12
   loop is disordered and the H12 helix starts at 536 rather than 537, while SERM structures keep that loop "more
   ordered" (Diennet 2023, doi 10.1101/2023.10.31.564979). Against this, a panel study finds SERM and SERD
   structures "very similar" (Min 2021, doi 10.1073/pnas.2101657118). **The degradation link is asserted, not
   shown:** the same papers measure both but state that how H12 destabilisation leads to degradation "remains
   unresolved", report "poor correlation with endogenous ERalpha degradation", or note that receptor loss "may not
   be required" for antagonism. Complete loss of H12 is documented only for ERbeta with ICI164,384 (RCSB 1HJ1:
   255 deposited, 204 modelled residues; Pike 2001, doi 10.1016/s0969-2126(01)00568-8) - a different receptor and
   species.
3. **The published ERalpha volumes are few, second-hand and method-free.** 450 A3 for the cavity left by removing
   the hormone (Tanenbaum 1998, doi 10.1073/pnas.95.11.5998), 450 A3 pocket with 250 A3 ligands for the E2 and DES
   complexes (Katzenellenbogen 2011, doi 10.1021/jm200801h), a 40 % pocket enlargement by a bulky substituent
   (Nettles 2007, doi 10.1038/sj.embor.7400963); a 73-structure survey, the ERbeta structure paper and the 2025
   apo-state paper give no volumes at all. No opened ERalpha source names a program, probe radius or grid. The only
   opened sources reporting pocket, ligand and occupancy together are VDR papers (660 A3 pocket with 375-392 A3
   ligands; 729-776 A3 cavity with 461-474 A3 ligands, occupancy "61%, 63%, 63%"), again with no program named.

## What this bears on in our results (full mapping in the table)

- **SV13 (C2a / C2b, the VDR superagonist point):** both comparisons were not testable in our data (0 partial-agonist
  entries; 3 agonist entries with an SV15 dpX). The literature answer for VDR is that H12 position does not separate
  superagonists from agonists, with potency differences of 20-40x and Calpha RMSD under 0.4 A - so a null here would
  agree with the field, and the slide can cite the literature instead of our (ungated) data.
- **SV13 (C1, SERD vs SERM):** not testable (SERD 4 entries vs the gate of 5). The literature expects the difference
  in H12 *mobility* and in H11-H12 loop order, not in whether 538-547 is modelled; our complete-H12 counts
  (agonist 53/55, SERM 19/53, SERD 5/9) are confounded with cluster, construct and resolution, and the opened
  sources show the same ligand classes crystallised with H12 docked in the AF-2 cleft.
- **SV13 Part C and SV12:** the SERD / pure-antiestrogen structures in the literature use the L536S + Cys->Ser
  construct that 243 of our 788 chains carry; the SERD representative entry 7MSA has a second, different ligand in
  chains B and D (RCSB entry page).
- **SV14 (references):** our E2 ligand volume 223.7 A3 is ~11 % below the published ~250 A3; our wild-type E2 pocket
  956.2 A3 is about twice the published 450 A3, while the lower mode of our bimodal E2 pocket (505-579 A3, 8 of 17
  chains) is close to it; our median filling 33.8 % is below the ~56 % implied by the published pair and the 61-63 %
  reported for VDR. None of these can be reconciled by method, because no opened source states one. The 40 % pocket
  enlargement by a bulky ligand is the literature counterpart of our ligand-pocket volume correlation
  (rho 0.44 ligands / 0.47 entries).
- **SV15:** the VDR sources show agonists differing 20-40x in potency with no H12 position difference, which is the
  pattern C2b would have tested against dpX (DES +0.63, J3Z -0.78).
- **SV01a / SV01b / SV02:** SERD-associated disorder at 528-535 and an H12 helix starting at 536 rather than 537
  (Diennet 2023) fall inside our mobility window 528-547 and next to our H12 window 538-547; one opened structure has
  H12 unmodelled in 2 of 4 chains of the same entry (Min 2021), the within-entry discordance our chain-level
  analysis assumes.

## Files

- docs/sv/literature_class.md (28.1 kB) - the question | claim | source (opened) | bears-on table (20 claims), the
  short answer per question, the list of sources that could not be opened, and the gaps.
- results/sv/SV08b/SUMMARY.md (this file)
- results/sv/SV08b/env.txt (1.2 kB)

Nothing enters SV09 from here unless its source was opened; every row of the table names what was opened.
