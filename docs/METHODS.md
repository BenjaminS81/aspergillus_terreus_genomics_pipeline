# Recorded methods and parameters

The source archive is the authority for implementation details. This summary separates recorded runs from proposed extensions; it does not silently upgrade historical tools or models.

| Analysis | Recorded settings | Scope / limit |
|---|---|---|
| Exact-66 single markers | IQ-TREE 2; `-m MFP -B 1000 -alrt 1000 -bnni -nt 1 -seed 20261111`; display outgroup `LOCAL_SRR10507650` | BenA 515, CaM 547, RPB2 1,014 retained columns |
| Exact-66 combinations | Partition file; `-m MFP+MERGE -B 1000 -alrt 1000 -bnni -nt 1 -seed 20261111` | Seven overlapping combinations are not seven independent gene votes |
| Complete 230-tip tree | IQ-TREE 2.4.0; partitioned `-m MFP+MERGE -B 1000 -alrt 1000 -bnni -wbt -nt 1 -seed 20261118`; outgroup `REF_436_Aspergillus_flavipes_NRRL302` | All three marker regions represented; occupancy alone does not prove single-copy orthology |
| Stored core | 4,661-locus concatenation, 66 tips | No comparable core bootstrap evidence in the original run |
| RF | Unrooted nonterminal split symmetric difference divided by total split count | Exact common leaves; additional short-edge sensitivity at >1.1e-6 |
| Clade recovery | Exact member-set separation by an unrooted split | Singleton N/A; recovery in the tree used to define the group is not independent validation |
| Projected bootstrap | Split frequency among 1,000 stored marker UFBoot trees restricted to 66 labels | Not equivalent to full reference/type-clade support |
| Quartets | Exhaustive 720,720 four-genome subsets; four-point unit-edge distances | Unresolved cases kept separate; literal binary short-edge resolutions can affect agreement |
| Tanglegrams | Rectangular phylograms; common numerical scale per pair; fixed core order | Rotations preserve unrooted splits and patristic distances; connector crossings are not biological evidence |
| Lightweight audit workers | One CPU, one numerical thread, Windows Job Object 2 GiB cap | Applies only to those workers, not all historical scripts |
| Saved BUSCO policy | 20 GB memory ceiling, no swap, one-core CPU quota in recorded service policy | No change to the separate Panaroo 100 GB cap; inspect current service state before future jobs |

For extraction, inspect the BLAST scripts, genomic rescue scripts, orthogroup mappings and retained-window checks together. Do not replace a duplicate locus with an arbitrary highest-scoring match. For absence, distinguish absent genomic sequence, assembly fragmentation, masking, model failure, insufficient similarity and true gene loss.

For additional branchwise evidence, IQ-TREE documents [gCF/sCF](https://iqtree.github.io/doc/Concordance-Factor) and [constrained topology/AU tests](https://iqtree.github.io/doc/Advanced-Tutorial#tree-topology-tests). These are proposed next analyses, not results manufactured during repository publication. [chewBBACA CDS input](https://chewbbaca.readthedocs.io/en/latest/user/getting_started/quick_start.html) is a candidate allele-calling route after fungal CDS validation. [OrthoFinder](https://pmc.ncbi.nlm.nih.gov/articles/PMC6857279/) supplies orthology/duplication evidence; gene-model copy counts alone are insufficient.
