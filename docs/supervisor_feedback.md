# Supervisor feedback

## Meeting 2026-10-28 (team notes, verbatim)

```
[10/28/26 1:08 PM] ayau🫧 (taylor's version): we can check distribution of completeness, and if completeness correlates with any experimental conditions/parameters, or is it biology?
we can try geometrical approach since it's more sensitive for things like the length of the helix (could affect ligand)
consensus alignment of secondary structures
take helices without the problem of completeness, find a stable core
color outliers in clustering (but no pattern is also possible)
check spikes on figure in slide 5
check literature if only H12 is significant? so far seems so
[10/28/26 1:13 PM] ayau🫧 (taylor's version): also, for friday we can use that ppt but we have to extract extra data and leave only the main findings etc
[10/28/26 1:15 PM] ㅤ`Ye: I also remember the point about somehow being able to compare incomplete residues, for example, taking the middle atoms (based on start and end of the helix) and applying methods on these data. Huge limitations exist but we are playing around
```

## Meeting notes received 2026-11-04

```
Correlation between ligand size and pocket size
APO is the baseline and do relative according to it, volume normalized to 100%.
Go back to standard natural structure that is binding to estrogen receptor and calculate ligand volume and normalize to 100%.
And then you can calculate the percentage of filling.

Different affinities since of different assays (in vitro, ...) and create different IC50, ... There are a lot of different systems and you need to think how to compare based on different systems.

6 different action class types and based on their DOIs, show in the plot. Somehow try to show that.

Prof is saying about superagonists in vdr but scientists didn't find differences in H12 between simple agonists.

Can antagonist kill the receptor by heavily distabilizing them, the H12. Good to show the picture of H12 and the ligand and compare between several class types. This is very important. Is difference on structure linked to functional something.

One idea: get the SMILES, InChiKey from api and using chembl api to get info.

Prof asked to predict affinities based on folding but he isn't sure whether we can finish.
```

## Meeting notes received 2026-11-11

```
Ye: Outliers in "how far the lid sits from H11", specifically from G4. Put ligand action class types to see correlation. Examine each outliers.
Summons mapping, Self organising maps.
UMAP of PCA.
Possibly hide H12 and do clusters again.
If you add a lot of guys in crystalization, it should be stable. So, receptor should be more dynamic inside water.
Professor says that the real treasure is hidden inside G4, two guys from G7, and gray outliers.
Basically, we need find an usual tree from forest.
Instead of RMSD, do analysis Phi/Psi, geometrical descriptors, Center of mass COM.
ISE GEOMETRICAL APPROACHES.
Dana: Biological variability (st dev of the structure ) or something significant ?
If too zoom in, we will see smth
Future work : you define the range of the helix, extension decrease of the helix tells about stability, define com center of mass virtual point for every helix , pairwise distanse, com for every protein , you see how this point is moving, combine movement (which parts are moving) which feature (rank) connected to which part (no superimposing needed)
Dana: Look at outliers, action type (even if partially) sample of structures. Agonist antagonist prove, api from chembl uncurated
```
