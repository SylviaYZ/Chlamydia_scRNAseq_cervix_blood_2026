# Cervical and blood single-cell RNA-seq analysis of Chlamydia infection

R scripts for single-cell RNA-seq analysis of cervical immune cells and circulating T cells during Chlamydia infection.

## Folder contents

- **`cervix/`** — Data processing, quality control, normalization, integration, broad cell-type annotation, differential expression, pathway enrichment, and CollecTRI/VIPER transcription factor analysis.
- **`cervix/CD8/`** — Cervical CD8 T-cell processing, normalization, integration, clustering, subtype annotation, composition, differential expression, pathway enrichment, and trajectory analysis.
- **`cervix/CD4/`** — Cervical CD4 T-cell processing, normalization, integration, clustering, subtype annotation, composition, differential expression, pathway enrichment, and trajectory analysis.
- **`blood/`** — Blood data processing, annotation, cell composition, and differential expression analysis.
- **`figures/`** — Code for manuscript and supplementary figures.
- **`data/`** — Cytokine values used in Figure 5 (`cytokines_range_adjusted_public.csv`), Chlamydia burden data used in Figure 9 (`BurdenforscRNA.xlsx`), and the STRING network SVG exports.

STRING network analysis was performed on the [STRING website](https://string-db.org/). Networks were exported as SVG files can be found under `data/`.

Comments in the analysis scripts indicate which exported outputs correspond to Supplementary Tables 3–12. Supplementary Table 1 combines QC outputs with study metadata, and Supplementary Table 2 contains a curated marker list. The manuscript workbooks also include manual formatting.

## Usage

Run the analysis scripts in numerical order within each folder, followed by the corresponding figure scripts. Run from the repository root so `data/` paths resolve. Update the other input and output paths consistently for your computer; each generated input has a comment naming its producing script. The two `01_read_10x.R` scripts read the GEO H5 files. Download Supplemental Table 1 from the manuscript supplementary materials and save it as `data/Supplemental Table 1.xlsx` to map GEO masked IDs to cohort IDs.

Figure 5 uses the public cytokine CSV: run `figures/05a` before `05b`; `05c` reads that CSV directly. Run `figures/S01_S02` before `S03_S04` to generate the participant plotting key. The STRING script formats the supplied SVG exports in `data/STRING/`; network selection and layout were performed in STRING.

Install the packages loaded by the selected scripts.

The dataset is available on GEO under accession [GSE334218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334218).

