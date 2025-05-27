# Welcome to CloEvoScore! <img src="www/logo2.png" width="200" height="200" align="right">








CloevoScore is an interactive, open-source tool for analyzing the evolutionary trajectories of patients based on gene copy number values. It provides two complementary approaches:

- The **Shiny app**, designed for small-scale, user-friendly analysis of up to **4 patients** with approximately **2000 genes** each.
- The **pipeline**, implemented in Docker, optimized for analyzing larger **cohorts of patients**, ensuring scalability and efficiency in processing high-throughput genomic data.

This repository includes:

- The [**Shiny app**](https://gmazz.shinyapps.io/CloEvoScore/) code for interactive analysis.
- The **pipeline** with all functions necessary for large-scale genomic evolution studies.
- Example **datasets** to help users get started.

With this tool, you can:

- Upload gene copy number data.
- Compute evolutionary trajectories and derive evolutionary scores.
- Visualize results through interactive plots and summary tables.
- Download the results for further analysis.

The **Shiny app** is ideal for exploratory analysis of small datasets, while the **pipeline** is recommended for comprehensive studies involving larger cohorts.

For more details on the pipeline, please refer to the reference paper:  
[Preprint available here](https://i.pinimg.com/originals/88/f5/aa/88f5aac99900cc1d075be33a06285db6.jpg).

**IMPORTANT:** This example uses genomic profiles of 4 multiple myeloma samples created solely for demonstrating the online version of this tool. 🧬

---

## Input Files Description 📂

### Gene Call Table
This file contains detailed information on the copy number values of individual genes for each patient in the study. These values represent the number of gene copies present in the genome and are essential for understanding genetic alterations associated with the disease being analyzed. This file serves as the foundation for identifying patterns of genomic amplification or deletion across the patient cohort.

#### Required columns:
- **chromosome:** The chromosome where the gene is located.
- **start:** Start position of the gene.
- **end:** End position of the gene.
- **gene_names:** Gene symbol.
- At least two columns per sample with copy number values at two time points (e.g., `PT_001_D`, `PT_001_R`).

#### Example:
| chromosome | start  | end    | gene_names | PT_001_D | PT_001_R |
|------------|--------|--------|------------|----------|----------|
| 1          | 34567  | 34890  | TP53       | 2.1      | 1.8      |
| 2          | 56789  | 57234  | BRCA1      | 3.0      | 2.5      |

---

### Ploidy and Purity Table (Optional)
This file provides two critical parameters used to refine the interpretation of copy number values:
- **Ploidy** is the correction factor multiplied by 100 to be applied to the baseline region of the copy number profile. It can also be achieved with another of our tools, [BOBaFIT](https://github.com/bioinformatic-seragnoli/BOBaFIT).
- **Purity** indicates the percetage of cancer cells within a given sample compared to normal cells.

Correcting the raw copy number values using ploidy and purity ensures a more accurate reflection of true genetic alterations.

#### Required columns:
- **pt_name:** Patient’s code (must match the Gene Call Table).
- **Ploidy_D** and **Ploidy_R:** Ploidy values for each time point.
- **Purity_D** and **Purity_R:** Purity values for each time point.

#### Example:
| pt_name  | Ploidy_D | Ploidy_R | Purity_D | Purity_R |
|----------|---------|---------|---------|---------|
| PT_001   | 50     | 4     | 80    | 75    |
| PT_002   | 0     | 10   | 85    | 78    |

---

### Disease-Related Gene Regions (Optional)
This file contains a curated list of genes that are particularly relevant to the disease under investigation. These genes are typically identified through computational analyses using tools such as **GISTIC2.0**. In the analysis, if these genes are not associated with any specific cluster, they will still be considered individually and contribute to the final score.

#### Required columns:
- **gene_names:** Gene symbol.
- **chromosome:** Chromosome where the gene is located.
- **start:** Start position of the gene.
- **end:** End position of the gene.
- **CNA:** Copy number alteration detected by GISTIC (e.g., `A` for Amplification, `D` for Deletion).

#### Example:
| gene_names | chromosome | start  | end    | CNA |
|-----------|----|--------|--------|-----|
| TP53      | 17 | 7565097 | 7590856 | D   |
| MYC       | 8  | 128748315 | 128753680 | A   |

### NB: All files must be provided in CSV, TSV, or TXT format.
---

## Running CloevoScore with Docker 🐳

You can run CloevoScore using Docker with the following commands:
```sh
docker pull magaiaa/cloevoscore
```
with test files:
```sh
docker run --rm -v $(pwd)/outputs:/cloevoscore/outputs -it magaiaa/cloevoscore \
--genes /cloevoscore/data/example/gene_table.txt \
--pp /cloevoscore/data/example/test_pp.txt \
--rel /cloevoscore/data/focal_loci_hg19.txt
```
with your files:
```sh
docker run --rm -v $(pwd):/cloevoscore/tmp -v $(pwd)/outputs/:/cloevoscore/outputs/ -it magaiaa/cloevoscore\
--genes /cloevoscore/tmp/user_genes.txt \
--pp /cloevoscore/tmp/user_ploidy_purity.tsv \
--rel /cloevoscore/tmp/user_disease-rel-genes.txt 
```

This command will:
- Mount the `outputs` directory to store results.
- Run the analysis using test files or user's files
- Compute evolutionary scores and generate results in the specified output directory.

The analysis could also be performed without correcting for purity and ploidy and/or including disease-related genes.

---

## License
This project is open-source and available under the [GPLv3](LICENSE).

---

## Citation
If you use CloevoScore in your research, please cite the reference paper linked above.

---

## Contact
For any questions or issues, please open an issue in this repository or email us at **bioinformatic.seragnoli@gmail.com**. 📩

