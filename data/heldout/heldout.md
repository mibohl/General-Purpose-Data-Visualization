# Held-out RT datasets (not seen by the DeepLC 4.5 default model)

Purpose: compare TabPFN in-context RT prediction with DeepLC 4.5's own
multi-head calibration on LC setups that DeepLC 4.5.0 (released 2026-09-11,
6,543 training setups) cannot have seen.

## How they were selected

- Source: PRIDE Archive v3 search API, about 40 keywords (proteomics,
  phosphorylation, acetylation, TMT, HLA, DIA-NN, MaxQuant, timsTOF, ...),
  sorted by publication date. That gave 1,198 projects published 2026-09-21 to
  2026-11-20. None of their PXD accessions occur in DeepLC 4.5's task names
  (`task_names.txt`, 958 PXDs). Each selected PXD was grepped again: 0 hits.
- Each project needed a PSM- or precursor-level table with RTs (<= 300 MB to
  download), and >= 3,000 unique peptidoforms at 1% FDR in one run.
- One run per dataset: the raw file with the most unique peptidoforms at 1%
  FDR. RT is the median over that run's identifications of each peptidoform
  (sequence + modifications), so no cross-run alignment is involved.
- Conversion: `python3 -I rtpoc/heldout_prepare.py --check` (one reader per
  format: MaxQuant evidence.txt, FragPipe psm.tsv, DIA-NN report.parquet, PEAKS
  db.psms.csv, Proteome Discoverer PSMs.txt). Output `data/heldout/<PXD>.csv`
  has the columns `seq,modifications,tr` in the Zenodo PEPREC convention:
  `pos|Name` with Unimod names, 0 = N-term, 1-based residues, -1 = C-term,
  empty if unmodified, fixed Carbamidomethyl (and fixed TMTpro on K) written
  explicitly, and tr in minutes.
- Check: for every dataset, 3,000 converted peptidoforms (all of them where
  there were fewer) were rebuilt as ProForma with psm_utils 1.5.5. Their
  theoretical masses all match the mass reported by the search engine (100%
  within 0.02 Da). So the Unimod mapping, including fixed mods, is right.
- Mod mapping: no peptidoform had to be dropped for an unmappable
  modification or a non-standard residue in any of the 8 chosen runs (0 / 0
  everywhere).
- FDR filters: MaxQuant uses evidence.txt (1% PSM FDR), keeping MS/MS-identified
  rows only (`*MSMS` types: no match-between-runs transfers, no second
  peptides) and dropping reverse hits. FragPipe psm.tsv is already filtered
  at 1%. DIA-NN keeps precursors with run-level Q.Value <= 0.01 and
  Lib.Q.Value <= 0.01. PEAKS uses the export at 1% peptide FDR. PD keeps
  PSMs with Percolator q <= 0.01.
- Column definitions: "% modified" is the share of peptidoforms carrying any
  modification other than Carbamidomethyl. Mod counts are the number of
  peptidoforms carrying at least one of that modification.

## Summary

| PXD | Label | Peptidoforms | % modified (non-CAM) | Modification types (peptidoforms) | Chromatography |
|---|---|---|---|---|---|
| PXD076643 | easy | 12,944 | 4.3% | Carbamidomethyl 2,345; Acetyl (prot N-term) 502; Oxidation 69 | nano RP C18 (bioZen Peptide XB-C18, 25 cm), 60 min, 0.1% FA |
| PXD079751 | easy | 92,685 | 8.2% | Carbamidomethyl 10,447; Oxidation 7,638 | nano RP C18 (Vanquish Neo), 120 min, 0.1% FA |
| PXD062983 | easy | 26,941 | 19.3% | Oxidation 5,012; Carbamidomethyl 3,196; Acetyl (N-term) 285 | nano RP C18 (25 cm, 1.5 µm ReproSil Saphir), 60 min, 0.1% FA |
| PXD050641 | hard | 14,275 | 79.4% | Crotonyl 11,231; Carbamidomethyl 1,043; Oxidation 668; Acetyl 121 | nano RP (nanoElute), 60 min; anti-Kcr enriched |
| PXD074962 | hard | 4,006 | 95.5% (70.0% excluding TMTpro) | TMTpro 3,400; Phospho 2,265; Oxidation 1,021; Carbamidomethyl 503; Acetyl 16 | nano RP (one of 12 high-pH RP fractions); Fe-NTA phospho, TMTpro 18-plex |
| PXD061731 | hard | 5,021 | 5.8% | Oxidation 176; Carbamidomethyl 99; Deamidated 46; Acetyl (K) 44; Methyl 21; Cysteinyl 15 | nano RP (not described in PRIDE); HLA-I peptides |
| PXD084949 | hard | 3,575 | 5.2% | Oxidation 185; Carbamidomethyl 1 | nano RP C18 (10 cm, 1.9 µm ReproSil-Pur), about 85 min; HLA-I peptides |
| PXD085246 | hard | 5,800 | 6.1% | Acetyl (prot N-term) 274; Oxidation 104 | nano RP (EASY-nLC 1200, not described further); endogenous peptidome |

3 easy, 5 hard. None of them uses non-reversed-phase online chromatography,
because none was available (see "Rejected"). The hard sets are hard because of
their chemistry: PTM-enriched, TMT-labeled, or non-tryptic.

Download: 352 MB in total for the kept files. For PXD050641, only
`txt/evidence.txt` (plus parameters/summary) was pulled from the 589 MB
`txt_LPCr.zip` with HTTP range requests; that member is 72 MB compressed and
261 MB once extracted. About 45 MB more was downloaded and then deleted
(listed under "Rejected"). On disk, `data/heldout/raw/` is 541 MB.

## Datasets

### PXD076643 (easy)
- Title: Comparison of arsenite and cadmium: Metal homeostasis and cellular stress response in human bronchial cells (Link et al., Int J Mol Sci 2026)
- Published: 2026-09-21. Lab: Food Chemistry and Toxicology, KIT Karlsruhe.
- File: `220822_ML_Cd_evidence.txt` (86.4 MB; MaxQuant 2.4.2.0)
- Run: `220822_ML_Cd_2_1` (12,944 peptidoforms; 12 runs, next best 12,732)
- Chromatography: UltiMate 3000 RSLCnano; Phenomenex bioZen Peptide XB-C18, 2.6 µm, 250 mm x 75 µm, 40 °C, 300 nL/min; 0.1% FA / ACN; 2 to 34% B in 60 min (80 min method). Observed tr 13.2 to 71.8 min.
- Instrument: Q Exactive Plus, DDA Top20. Search: MaxQuant/Andromeda; fixed CAM; variable Met-ox and protein N-term acetyl; tryptic.
- Peptidoforms: 12,944. Modified (non-CAM): 4.3%. Mods: Carbamidomethyl 2,345; Acetyl 502; Oxidation 69.
- Label **easy**: standard tryptic, unlabeled human cell-line proteome on acidic C18, common mods only.

### PXD079751 (easy)
- Title: Cold-induced structural remodeling of the translation machinery in a thermophilic eukaryote
- Published: 2026-10-08. Lab: ZMBH Core Facility, Heidelberg University.
- File: `report.parquet` (86.3 MB; DIA-NN 2.2.0)
- Run: `EC-DIA120-25-72-04` (92,685 peptidoforms; 8 runs, next best 92,641)
- Chromatography: Vanquish Neo nanoLC, RP C18 (column per PMID 40846850, not named in PRIDE); 120 min gradient (4 to 32% B over 100 min, to 49% over 20 min). Observed tr 0.5 to 127.5 min.
- Instrument: Orbitrap Eclipse, DIA (40 variable windows). Search: DIA-NN 2.2.0 library-free; fixed CAM; variable Met-ox; tryptic.
- Peptidoforms: 92,685. Modified (non-CAM): 8.2%. Mods: Carbamidomethyl 10,447; Oxidation 7,638.
- Label **easy**: standard tryptic, unlabeled acidic C18 with a long gradient. The organism (*Chaetomium thermophilum*) is unusual, but the chemistry is standard.

### PXD062983 (easy)
- Title: Synthesis of oocyte-specific proteins during mouse embryogenesis and biological significance thereof (Nolte et al., BMC Genomics 2026)
- Published: 2026-10-09. Lab: MPI for Molecular Biomedicine, Münster.
- File: `report.parquet` (33.4 MB; DIA-NN 2.0.2)
- Run: `20250312_1350_Michele_Oocytes_saphir25cm_dia6L_01_Slot2-73_1_6162` (26,941 peptidoforms; 12 runs, next best 26,821)
- Chromatography: nanoElute 2; 25 cm x 75 µm, home-packed 1.5 µm ReproSil Saphir 100 C18 (CoAnn); 300 nL/min; 0.1% FA; 4 to 27% B over 60 min. Observed tr 0.8 to 63.6 min.
- Instrument: timsTOF Pro 2, diaPASEF. Search: DIA-NN 2.0.2 library-free; fixed CAM; variable Met-ox and protein N-term acetyl; tryptic.
- Peptidoforms: 26,941. Modified (non-CAM): 19.3%. Mods: Oxidation 5,012; Carbamidomethyl 3,196; Acetyl 285.
- Label **easy**: standard tryptic, unlabeled acidic C18, common mods only. It adds a Bruker/timsTOF platform and a short 60 min gradient.

### PXD050641 (hard)
- Title: Fluctuation of lysine crotonylation confers azole tolerance in *Candida auris*
- Published: 2026-10-08. Lab: Changzheng Hospital, Shanghai (PTM Bio-style antibody enrichment).
- File: `txt/evidence.txt` from `txt_LPCr.zip`. The zip is 588.5 MB; only this member was range-downloaded (72.1 MB compressed, 261 MB uncompressed). MaxQuant 1.6.6.0.
- Run: `WA083LPCr_F2852A_1_Slot2-4_1_988` (14,275 peptidoforms; 24 runs, next best 13,708)
- Chromatography: nanoElute UHPLC, RP (column not named), 450 nL/min, 0.1% FA. Gradient: 7 to 24% B over 0 to 42 min, 24 to 32% B to 54 min, then 80% B to 60 min. Observed tr 5.1 to 60.0 min.
- Instrument: timsTOF Pro, DDA-PASEF. Search: MaxQuant; fixed CAM; variable Met-ox, protein N-term acetyl, and "Crotony(K)". This name maps to Unimod Crotonyl (+68.026); the mass check confirms it. Trypsin/P.
- Peptidoforms: 14,275. Modified (non-CAM): 79.4%. Mods: Crotonyl 11,231; Carbamidomethyl 1,043; Oxidation 668; Acetyl 121.
- Label **hard**: anti-crotonyl-lysine enriched, so 79% of peptidoforms carry Kcr. This PTM is nearly absent from DeepLC's training data (only prosit_ptm_2020 has about 300 crotonyl peptides). Each Kcr also blocks tryptic cleavage, adding internal K.
- Alternatives in the same project, not used: acetyl, succinyl, malonyl and ubiquitin enrichments (`txt_LPAc/LPSc/LPMa/LPUb.zip`; for example, LPAc's evidence.txt is 55 MB compressed).

### PXD074962 (hard)
- Title: AAV-MYO mediated overexpression of RBM20 in mice
- Published: 2026-11-19. Lab: EMBL Heidelberg (Genome Biology / Proteomics Core).
- File: `P4096_Set1_PP_psm.tsv` (49.8 MB; FragPipe 23 / MSFragger 4.3)
- Run: `Zelda_251106_P4096_JS_JK_Set1_TMT16_PhosphoP_R1_FR07` (4,006 peptidoforms; 12 fractions, next best 3,879)
- Chromatography: phosphopeptides were enriched with Fe-NTA IMAC and labeled with TMTpro 18-plex. They were then fractionated offline by high-pH RP (Gemini C18, pH 10, 48 fractions pooled into 12 by concatenation). The chosen run is one fraction on online nano RP. The online column and gradient are not described in the PRIDE record (the protocol text is truncated); observed tr 11.5 to 87.5 min.
- Instrument: Orbitrap Fusion Lumos, DDA. Search: MSFragger; fixed CAM and TMTpro (K); variable TMTpro (peptide N-term), Met-ox, protein N-term acetyl, and Phospho (STY); tryptic.
- Peptidoforms: 4,006. Modified (non-CAM): 95.5% (70.0% excluding TMTpro). Mods: TMTpro 3,400; Phospho 2,265; Oxidation 1,021; Carbamidomethyl 503; Acetyl 16.
- Label **hard**: TMTpro-labeled phosphopeptides, 57% phosphorylated. TMTpro sits on every K (fixed) and on most peptide N-termini (variable): 85% of peptidoforms carry TMTpro, which shifts their hydrophobicity. It is also a single high-pH fraction (concatenated, so it still spans the gradient).

### PXD061731 (hard)
- Title: H3K27M Mutation Reshapes the Antigenic Landscape in Diffuse Midline Glioma (Shamekhi et al., Cancers 2026)
- Published: 2026-10-26. Lab: Hudson Institute / Monash University (Faridi lab).
- File: `db.psms.csv` (13.1 MB; PEAKS 11)
- Run: `EXP1TL20240319_Tima_DDA_noFAIMS_BPT_II.raw` (5,021 peptidoforms; 12 runs, next best 4,596)
- Chromatography: HLA-I peptides (W6/32 immunoprecipitation, C18 SPE cleanup) on nano RP. The LC column and gradient are not described in the PRIDE record; observed tr 16.1 to 111.1 min.
- Instrument: Orbitrap Exploris 480, DDA, no FAIMS. Search: PEAKS 11 (de novo-assisted database search), unspecific digestion, 1% peptide FDR.
- Peptidoforms: 5,021 (length 8 to 25, median 9). Modified (non-CAM): 5.8%. Mods: Oxidation 176; Carbamidomethyl 99; Deamidated 46; Acetyl (K) 44; Methyl (K/R) 21; Cysteinyl 15. Cysteines are not uniformly alkylated: CAM is written only where PEAKS reports it.
- Label **hard**: non-tryptic HLA class I immunopeptides (short, no C-terminal K/R constraint).

### PXD084949 (hard)
- Title: The *Mycobacterium tuberculosis* ESX-1 secretion system controls antigen presentation on MHC-I in human macrophages
- Published: 2026-09-30. Lab: MIT Biological Engineering.
- File: `070224_hMDM_pathMHCI_WT_DDA_PSMs.txt` (5.7 MB; Proteome Discoverer, Sequest HT + Percolator)
- Run: `070224_hMDM_pathMHCI_WT_DDA.raw` (the only run in this file; 3,575 peptidoforms at Percolator q <= 0.01)
- Chromatography: UltiMate 3000 RSLCnano, in-house 10 cm column of 1.9 µm ReproSil-Pur C18 with an integrated emitter, column oven. Gradient: 6 to 25% B (70% ACN, 0.1% FA) over 75 min, then 25 to 45% over 5 min, 45 to 100% over 5 min (about 88 min). Observed tr 4.4 to 113.0 min.
- Instrument: Orbitrap Exploris 480, DDA. Search: PD Sequest HT, unspecific; Met-ox variable.
- Peptidoforms: 3,575 (length 8 to 13, median 9). Modified (non-CAM): 5.2%. Mods: Oxidation 185; Carbamidomethyl 1.
- Label **hard**: non-tryptic HLA-I peptides from primary human macrophages. It comes from a different lab, LC and search engine than PXD061731, but it is close to the size threshold.

### PXD085246 (hard)
- Title: SRB-induced peptidomic profiling of *Hermetia illucens* larvae for antimicrobial peptide discovery
- Published: 2026-11-09. Lab: Shandong University.
- File: `evidence.txt` (5.5 MB; MaxQuant 1.6.15.0)
- Run: `XA05698LQPe_experiment` (5,800 peptidoforms; 2 runs, next best 4,812)
- Chromatography: EASY-nLC 1200, RP. The column and gradient are not described in the PRIDE record; observed tr 0.9 to 89.4 min.
- Instrument: Q Exactive HF-X, DDA. Search: MaxQuant unspecific search (8 to 25 aa); no fixed mods (no alkylation, so C is left unmodified); variable Met-ox and protein N-term acetyl.
- Peptidoforms: 5,800 (length 8 to 25, median 16). Modified (non-CAM): 6.1%. Mods: Acetyl 274; Oxidation 104.
- Label **hard**: endogenous, undigested insect peptidome (non-tryptic ends, unalkylated Cys). The sequences are far from DeepLC's tryptic training distribution.

## Considered and rejected

- **No online non-RP setup found.** I searched all 1,198 project protocols
  for HILIC, PGC/Hypercarb, ERLIC, SAX/SCX, ZIC, CZE and amide columns. Every
  hit used HILIC/SCX/SAX only offline (glycopeptide enrichment, fractionation,
  bead cleanup) before RP nanoLC, or was not a bottom-up peptide dataset:
  - PXD078897: HILIC peptide mapping of a few biopharmaceuticals, xlsx exports only
  - PXD076819: nanoLC x CZE top-down proteoforms
  - PXD080737, PXD074760: ZIC-HILIC glyco-enrichment, Byonic/xlsx outputs
  - PXD045462: offline HILIC fractions, pep.xml/mzXML only
  - PXD080891: HILIC protein purification
- **PXD073597** (Ubi-SCAPE, GG): report.parquet is 858 MB and report.tsv is
  6 GB, both over 300 MB. It is also a protein-level TUBE enrichment, so GG
  peptides are a small minority.
- **PXD024429** (USP8, GG): this is a FLAG/HA IP, not diGly enrichment
  (GlyGly sites table is 0.7 MB). txt.zip is 617 MB.
- **PXD084146** (mouse liver, Q Exactive): its only real run has 2,247
  peptidoforms, under 3,000. evidence.txt (1.1 MB) was downloaded and then
  deleted.
- **PXD084949** fractionated donor files (`100423_Donor_N_MHC_I_*_PSMs.txt`):
  each of the 12 fractions has at most 1,457 peptidoforms. The mock file
  (43 MB) was downloaded and then deleted. The single-shot DDA file was used
  instead.
- **PXD078769** (report.tsv 701 MB) and **PXD071507** (report.tsv 1.24 GB):
  DIA-NN reports over 300 MB.
- **PXD065531** (Exploris 480 DIA, report.parquet 85 MB), **PXD083668**
  (bovine liver phospho, MaxQuant evidence 46.5 MB), **PXD075865** (tomato
  fruit peptidomics, FragPipe psm.tsv 8 MB, timsTOF Ultra) and **PXD080009**
  (TMT16 proteome, psm.tsv 21 MB): all viable. They were not needed after
  the 8-dataset cap and were not downloaded. They are the first backups
  (PXD075865 for non-tryptic, PXD083668 for phospho).
- **PXD076213** (bacteriophage acetylome): only a site-level table, with no
  PSM-level RTs.
- **PXD070942, PXD070850** (IP interactomes, evidence.txt about 3 MB): too
  small to give 3,000 peptidoforms in one run.
- **PXD061731 / PXD084949 / PXD085246** all lack some LC details in PRIDE
  (column or gradient). They were kept because only within-run RTs are used.
  The observed RT range is given instead.
