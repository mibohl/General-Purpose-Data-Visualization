# LC metadata coverage

## 1. DeepLC 4.5 training collection (PRIDE protocol text, automatic extraction)

958 PRIDE projects behind 5,553 of 6,543 DeepLC 4.5 setup heads. The other 990 heads (MassIVE and named datasets) were not looked up and count as missing below. Publications were not mined for these projects, although 89% have one, so these numbers are a lower bound on what is recorded anywhere.

| field                     | projects with value   | heads with value   | extraction precision (hand check)   |
|:--------------------------|:----------------------|:-------------------|:------------------------------------|
| column chemistry          | 59%                   | 58%                | ~100%                               |
| column length             | 50%                   | 54%                | ~95%                                |
| inner diameter            | 33%                   | 36%                | ~88%                                |
| particle size             | 35%                   | 44%                | ~88%                                |
| gradient length           | 60%                   | 59%                | ~95%                                |
| gradient composition (%B) | 37%                   | 46%                | ~65% (not used)                     |
| flow rate                 | 42%                   | 52%                | ~85%                                |
| column temperature        | 11%                   | 7%                 | ~85%                                |
| additive / ion pairing    | 53%                   | 54%                | ~85%                                |

Fields found per project (of 9): 0: 268, 1: 70, 2: 50, 3: 66, 4: 47, 5: 56, 6: 150, 7: 129, 8: 100, 9: 22

- No LC parameter at all: 268 projects (28%).
- Usable (chemistry + gradient length + at least 2 of length, ID, particle, flow; rule fixed before counting): 415 projects (43%), 3,340 heads (51%).
- Usable and one unambiguous gradient length (the training set of Competitor B): 335 projects (35%), 1,245 heads.

## 2. Datasets with observed RTs in this repo (curated by hand from PRIDE and publications)

| role                           | dataset                      | pxd       | fields_filled   | usable_for_B   | provenance                                                                                 |
|:-------------------------------|:-----------------------------|:----------|:----------------|:---------------|:-------------------------------------------------------------------------------------------|
| Zenodo (Part 1 / DeepLC paper) | hela_hf_psms_aligned         | PXD006932 | 1 / 9           | False          | DeepLC Data availability ref 42 (HeLa hf); Zenodo 5171792 README table maps hela_hf to PXD |
| Zenodo (Part 1 / DeepLC paper) | mod_fixed_mods               | PXD000612 | 6 / 9           | False          | PRIDE PXD000612 sampleProcessingProtocol; DeepRT+ GitHub data/README.md (mod.txt = HeLa, S |
| Zenodo (Part 1 / DeepLC paper) | ATLANTIS_SILICA_fixed_mods   | nan       | 6.5 / 9         | False          | DeepRT+ GitHub data/README.md (source = Spicer et al.); Spicer & Krokhin 2018 J Chromatogr |
| Zenodo (Part 1 / DeepLC paper) | LUNA_SILICA_fixed_mods       | nan       | 6.5 / 9         | False          | DeepRT+ GitHub data/README.md (source = Spicer et al.); Spicer & Krokhin 2018 J Chromatogr |
| Zenodo (Part 1 / DeepLC paper) | LUNA_HILIC_fixed_mods        | nan       | 6.5 / 9         | False          | DeepRT+ GitHub data/README.md (source = Spicer et al.); Spicer & Krokhin 2018 J Chromatogr |
| Zenodo (Part 1 / DeepLC paper) | Xbridge_fixed_mods           | nan       | 6.5 / 9         | False          | DeepRT+ GitHub data/README.md (source = Spicer et al.); Spicer & Krokhin 2018 J Chromatogr |
| Zenodo (Part 1 / DeepLC paper) | SCX_fixed_mods               | nan       | 6.5 / 9         | False          | Gussakovsky 2017 abstract (Polysulfoethyl A; offline SCX-RPLC of S. cerevisiae; Europe PMC |
| Zenodo (Part 1 / DeepLC paper) | pancreas_psms_aligned        | PXD010154 | 8 / 9           | True           | Wang 2019 Methods 'On-line LC-MS/MS' (PMC6379049, via NCBI efetch); PRIDE PXD010154; DeepL |
| Zenodo (Part 1 / DeepLC paper) | yeast_120min_psms_aligned    | PXD003472 | 8 / 9           | True           | Jarnuczak 2016 'Mass Spectrometry Data Acquisition and Processing' and Table 1 (Yeast_HF_1 |
| Zenodo (Part 1 / DeepLC paper) | yeast_60min_psms_aligned     | PXD003472 | 8 / 9           | True           | Jarnuczak 2016 'Mass Spectrometry Data Acquisition and Processing' and Table 1 (Yeast_HF_6 |
| Zenodo (Part 1 / DeepLC paper) | arabidopsis_psms_aligned     | PXD008812 | 8 / 9           | False          | PRIDE PXD008812 sampleProcessingProtocol; DeepLC ref 58; Zenodo 5171792 README.md table +  |
| Zenodo (Part 1 / DeepLC paper) | dia_fixed_mods               | PXD000953 | 8 / 9           | True           | Rosenberger 2014 Methods 'DDA mass spectrometry for spectral library generation' (PMC43225 |
| Zenodo (Part 1 / DeepLC paper) | PXD005573_mcp                | PXD005573 | 8 / 9           | True           | Bruderer 2017 Methods 'Mass Spectrometric Acquisition' (PMC5724188) and supplemental PDF ' |
| Zenodo (Part 1 / DeepLC paper) | unmod_fixed_mods             | nan       | 8.5 / 9         | True           | Nagaraj 2012 Experimental Procedures 'LC-MS/MS' (PMC3316726); DeepRT+ GitHub data/README.m |
| Zenodo (Part 1 / DeepLC paper) | hela_lumos_1h_psms_aligned   | PXD013477 | 9 / 9           | True           | PRIDE PXD013477 sampleProcessingProtocol (EASY-nLC 1200; PicoFrit 75 um x 50 cm; ReproSil- |
| Zenodo (Part 1 / DeepLC paper) | hela_lumos_2h_psms_aligned   | PXD013477 | 9 / 9           | True           | PRIDE PXD013477 sampleProcessingProtocol (EASY-nLC 1200; PicoFrit 75 um x 50 cm; ReproSil- |
| Zenodo (Part 1 / DeepLC paper) | plasma_lumos_1h_psms_aligned | PXD013477 | 9 / 9           | True           | PRIDE PXD013477 sampleProcessingProtocol (EASY-nLC 1200; PicoFrit 75 um x 50 cm; ReproSil- |
| Zenodo (Part 1 / DeepLC paper) | plasma_lumos_2h_psms_aligned | PXD013477 | 9 / 9           | True           | PRIDE PXD013477 sampleProcessingProtocol (EASY-nLC 1200; PicoFrit 75 um x 50 cm; ReproSil- |
| Zenodo (Part 1 / DeepLC paper) | proteometools_library        | PXD004732 | 9 / 9           | True           | Zolg 2017 Online Methods 'Nanoscale liquid chromatography' (PMC5868332, via NCBI efetch);  |
| Zenodo (Part 1 / DeepLC paper) | prosit_ptm_2020              | PXD009449 | 9 / 9           | True           | Zolg 2018 Experimental Procedures (PMC6126394, PMC HTML); PRIDE PXD009449; DeepLC ref 30   |
| held-out                       | PXD085246                    | PXD085246 | 0 / 9           | False          | PRIDE PXD085246 sampleProcessingProtocol only; no paper found (Europe PMC, Crossref, Nov 2 |
| held-out                       | PXD050641                    | PXD050641 | 4 / 9           | False          | PRIDE PXD050641 sampleProcessingProtocol only; no paper found (Europe PMC, Crossref, web s |
| held-out                       | PXD074962                    | PXD074962 | 6 / 9           | True           | PRIDE PXD074962 (enrichment, labeling, offline fractionation); online LC from companion pa |
| held-out                       | PXD084949                    | PXD084949 | 6 / 9           | True           | PRIDE PXD084949 dataProcessingProtocol 'DDA MS analyses'; same text in the lab's earlier p |
| held-out                       | PXD079751                    | PXD079751 | 7 / 9           | True           | PRIDE PXD079751 sampleProcessingProtocol (gradient; 'as described before [PMID 40846850]') |
| held-out                       | PXD062983                    | PXD062983 | 8 / 9           | True           | PRIDE PXD062983 sampleProcessingProtocol; confirmed by paper Methods 'Proteome analysis of |
| held-out                       | PXD061731                    | PXD061731 | 8 / 9           | True           | Shamekhi 2026 section 2.6 'LC-MS/MS Data Acquisition' (PMC13664579); PRIDE PXD061731 (IP o |
| held-out                       | PXD076643                    | PXD076643 | 9 / 9           | True           | PRIDE PXD076643 sampleProcessingProtocol; confirmed by paper section 3.6.2 'Shotgun Proteo |

Held-out setups usable for Competitor B: 6 of 8. PXD074962's online-LC values come only from a companion paper of the same core facility (no source for this dataset itself); PXD050641 has no column information; PXD085246 has no LC parameter at all, so Competitor B is untestable there.
