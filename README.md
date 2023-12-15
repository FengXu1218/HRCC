# README

## Abstract

In each cell cycle, cells perform complicated processes to ensure the accurate replication of DNA and partitioning of different types of material. It is important to recognize this cycle as a multiphase progress and mark each of the phases to allow detailed investigations. Advancements in RNA sequencing and analysis of individual cells give us the opportunity to explore unprecedented details of gene transcription during the cell cycle. Here, we show that partitioning the cell cycle into 9 or more phases can provide more information on the cycle progress than the typical method, which only uses G1, S, G2, and M phases. In scRNAseq data, we identified 152 periodic marker genes and found that the temporal order of their transcriptional waves was highly consistent across various cell types. Clustering and further analysis of the transcriptions suggested that groups of individual cells fell into different high-resolution phases of the cell cycle. The different high-resolution phases have transcriptional signatures that are different from each other and are characterized by essential cell cycle genes, such as *CCNA2*, *CCNB1*, *CCNE1*, and *PCNA*. We therefore developed multiphase frameworks to further partition the cell cycle to achieve higher resolution and obtained efficient classifiers to determine the high-resolution phases of individual cells by machine learning.

## DataSets

There were a total of 12 datasets involved in this study, encompassing various cell lines and experimental conditions (Supplementary Table S1). They consisted of HeLa cells of the HeLa CCL2 cell line at the 9th passage (HeLaCCL2_2019a) and the 14th passage (HeLaCCL2_2019b), HeLa S3 cells with AGO2 knockout (HeLaS3_2020KO) and the corresponding control (HeLaS3_2020a, HeLaS3_2020b) , human fibroblasts (hFib_2022), human neural progenitor cells (hNPC_2021), and human cells of embryoid bodies exposed to nicotine (hESC_2019b) and the corresponding control (hESC_2019a). Furthermore, three additional selected datasets were obtained for the phase labels of the cell cycle by using fluorescent ubiquitination-based cell cycle indicators (FUCCI51), i.e., human embryonic stem cell lines (hESC_2015), human osteosarcoma cell lines (hU2OS_2020)52, and induced pluripotent stem cell lines (hiPSC_2019).

## R codes

All the code for data analysis and figure generation is accessible in the `./R` directory.

## Rdata

Partial intermediate data generated during the analysis process, as well as essential supplementary information, can be accessed in the `./data` directory.
