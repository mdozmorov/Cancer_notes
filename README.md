# A continually expanding collection of cancer genomics notes and data

Issues with suggestions and pull requests are welcome!


* [Survival analysis](#suvival-analysis)
* [Cancer driver genes](#cancer-driver-genes)
* [Drugs](#drugs)
* [Consortia](#consortia)
* [Tools](#tools)
* [Misc](#misc)

## Survival analysis

- `cBioPortal` - The cBioPortal for Cancer Genomics provides visualization, analysis and download of large-scale cancer genomics data sets. OncoPrint mutation plots, differential expression, coexpression, survival. Compare gene expression with copy number variation. http://www.cbioportal.org/



- `R2` - Genomics Analysis and Visualization Platform. Gene-centric, survival analysis, collection of preprocessed microarray studies. http://hgserver1.amc.nl/

- `KM plotter` - Gene-centric, customizable survival analysis for breast, ovarian, lung, gastric cancers. http://kmplot.com/
    - Györffy, Balazs, Andras Lanczky, Aron C. Eklund, Carsten Denkert, Jan Budczies, Qiyuan Li, and Zoltan Szallasi. “An Online Survival Analysis Tool to Rapidly Assess the Effect of 22,277 Genes on Breast Cancer Prognosis Using Microarray Data of 1,809 Patients.” Breast Cancer Research and Treatment 123, no. 3 (October 2010): 725–31. https://doi.org/10.1007/s10549-009-0674-9.

- `The Human Protein Atlas` - Gene- and protein expression data in multiple cancer tissues, cell lines. Easy one-gene search, summary of tissue-specific expression, survival significance. http://www.proteinatlas.org/
    - Uhlen, Mathias, Cheng Zhang, Sunjae Lee, Evelina Sjöstedt, Linn Fagerberg, Gholamreza Bidkhori, Rui Benfeitas, et al. “A Pathology Atlas of the Human Cancer Transcriptome.” Science (New York, N.Y.) 357, no. 6352 (August 18, 2017). doi:10.1126/science.aan2507. http://science.sciencemag.org/content/357/6352/eaan2507. Rich downloadable data - tissue-specific gene expression in cancer and normal, isoform expression, protein expression. http://www.proteinatlas.org/about/download. Supplementary material http://science.sciencemag.org/content/suppl/2017/08/16/357.6352.eaan2507.DC1  
        - `Table S2` - summary of tissue specific expression for each gene, in normal and cancer tissues.  
        - `Table S6` - summary of survival prognostic value, with a simple "favorable/unfavorable" label for each gene. Each worksheet corresponds to a different cancer.  
        - `Table S8` - per-gene summary, in which cancers it is prognostic of survival.  

- `PRECOG` - PREdiction of Clinical Outcomes from Genomic Profiles. Gene-centric, quick overview of survival effect of a gene across all cancers, KM plots. https://precog.stanford.edu
    - Gentles, Andrew J., Aaron M. Newman, Chih Long Liu, Scott V. Bratman, Weiguo Feng, Dongkyoon Kim, Viswam S. Nair, et al. “The Prognostic Landscape of Genes and Infiltrating Immune Cells across Human Cancers.” Nature Medicine 21, no. 8 (August 2015): 938–45. https://doi.org/10.1038/nm.3909. - TCGA pan-cancer survival analysis PRECOG, CIBERSORT. 39 cancers. Intro into heterogeneity. Z-score description. Batch effect does not significantly affect z-scores. 2/3 prognostic genes shared across cancers. AutoSOME clustering method

- `GEPIA` - single- and multiple-gene analyses of TCGA data. Gene expression in different tumor-normal comparisons, differentially expressed genes, correlation analysis, similar genes, survival analysis. http://gepia.cancer-pku.cn/
    - Zefang Tang et al., “GEPIA: A Web Server for Cancer and Normal Gene Expression Profiling and Interactive Analyses,” Nucleic Acids Research 45, no. W1 (July 3, 2017): W98–102, https://doi.org/10.1093/nar/gkx247. - TCGA and GTEX web interface. Classical analyses - differential expression analysis, profiling plotting, correlation analysis, patient survival analysis, similar gene detection and dimensionality reduction analysis. http://gepia.cancer-pku.cn/

<!--
- `PrognoScan`, Gene-centric, survival effect of a gene in cancer studies from GEO. http://dna00.bio.kyutech.ac.jp/PrognoScan/
-->

- `UALCAN` - Gene-centric, tumor-normal expression, survival analusis, TCGA cancers. http://ualcan.path.uab.edu/
    - Chandrashekar DS, Bashel B, Balasubramanya SAH, Creighton CJ, Rodriguez IP, Chakravarthi BVSK and Varambally S. UALCAN: A portal for facilitating tumor subgroup gene expression and survival analyses. Neoplasia. 2017 Aug;19(8):649-658. doi: 10.1016/j.neo.2017.05.002 [PMID:28732212]

- `Project Betastasis` - Gene-centric, survival analysis, gene expression, select cancer studies. http://www.betastasis.com/

- `OncoLnc` - Gene-centric, survival analysis in any TCGA cancer. http://www.oncolnc.org/


## Cancer driver genes

- Cancer Gene Census (CGC), download [COSMIC](http://cancer.sanger.ac.uk/cosmic/download)
    - Hudson, T. J. et al. International network of cancer genome projects. Nature 464, 993–8 (2010).
    - `data/Census_all*.csv` - [The cancer Gene Census](http://cancer.sanger.ac.uk/census)
    - `data/COSMIC_genes.txt` - Genes sorted by the number of records associated with them. Obtained using `zcat <CosmicCompleteTargetedScreensMutantExport.tsv.gz | sed '1d' | cut -f1 | sort | uniq -c | sort -k1 -r -n > COSMIC_genes.txt`
    - `data/CosmicCodingMuts.vcf.gz` - VCF file of all coding mutations in the current release (release v83, 7th November 2017).

- Tumor suppressor gene database (TSGene), https://bioinfo.uth.edu/TSGene/
    - Zhao, M., Sun, J. & Zhao, Z. TSGene: a web resource for tumor suppressor genes. Nucleic Acids Res, 41(Database issue), D970–6 (2013).
    - Download various lists of tumor suppressor genes, https://bioinfo.uth.edu/TSGene/download.cgi

- `OncoScape` - Genes with oncogenic/tumor suppressor/combined scores as a sum contribution from gene expression, somatic mutations, DNA copy-number and methylation as well as data from shRNA knock-down screens. http://oncoscape.nki.nl/
    - Schlicker, Andreas, Magali Michaut, Rubayte Rahman, and Lodewyk F. A. Wessels. “OncoScape: Exploring the Cancer Aberration Landscape by Genomic Data Fusion.” Scientific Reports 6 (20 2016): 28103. https://doi.org/10.1038/srep28103.

- `data/Bailey_2018_cancer_genes.xlsx` - Table S1, consensus list of cancer driver genes.
	- Bailey, Matthew H., Collin Tokheim, Eduard Porta-Pardo, Sohini Sengupta, Denis Bertrand, Amila Weerasinghe, Antonio Colaprico, et al. “Comprehensive Characterization of Cancer Driver Genes and Mutations.” Cell 173, no. 2 (April 5, 2018): 371-385.e18. https://doi.org/10.1016/j.cell.2018.02.060. - Pan-Cancer mutation analysis. Combined use of 26 tools (https://www.cell.com/cell/fulltext/S0092-8674(18)30237-X#secsectitle0075, description of each tool in Methods) on harmonized data. 299 cancer driver genes, >3,400 putative missense driver mutations. Table S6 - excluded TCGA samples.

- `data/TARGET_db_v3_02142015.xlsx` - TARGET (tumor alterations relevant for genomics-driven therapy) is a database of genes that, when somatically altered in cancer, are directly linked to a clinical action. TARGET genes may be predictive of response or resistance to a therapy, prognostic, and/or diagnostic. https://software.broadinstitute.org/cancer/cga/target

- `data/Tokheim_2016_cancer_driver_genes.xlsx` - Dataset S2: Predicted driver genes by various number of methods
    - Tokheim, Collin J., Nickolas Papadopoulos, Kenneth W. Kinzler, Bert Vogelstein, and Rachel Karchin. “Evaluating the Evaluation of Cancer Driver Genes.” Proceedings of the National Academy of Sciences 113, no. 50 (December 13, 2016): 14330–35. https://doi.org/10.1073/pnas.1616440113. - 20/20+ machine learning method, ratiometric approach to predict cancer driver genes. Performance comparison of other methods, 20/20+, TUSON, OncodriveFML and MutsigCV are the top performers. https://github.com/KarchinLab/2020plus


## Drugs

- `GDA` - Genomics and Drugs integrated Analysis. The Genomics and Drugs integrated Analysis portal (GDA) is a web-based tool that combines NCI60 uniquely large number of drug sensitivity data with CCLE and NCI60 gene mutation and expression profiles. Gene-to-drug and reverse analysis. http://gda.unimore.it/

- `OncoKB` - OncoKB cancer gene database, different levels of evidence, fully downloadable. http://oncokb.org
    - Chakravarty, Debyani, Jianjiong Gao, Sarah M. Phillips, Ritika Kundra, Hongxin Zhang, Jiaojiao Wang, Julia E. Rudolph, et al. “OncoKB: A Precision Oncology Knowledge Base.” JCO Precision Oncology 2017 (July 2017).

- `CancerRxGene` - Drug-gene targets. Lots of drug sensitivity information. http://www.cancerrxgene.org/
    - Yang, Wanjuan, Jorge Soares, Patricia Greninger, Elena J. Edelman, Howard Lightfoot, Simon Forbes, Nidhi Bindal, et al. “Genomics of Drug Sensitivity in Cancer (GDSC): A Resource for Therapeutic Biomarker Discovery in Cancer Cells.” Nucleic Acids Research 41, no. Database issue (January 2013): D955-961. https://doi.org/10.1093/nar/gks1111.

- `CTRP` - The Cancer Therapeutics Response Portal (CTRP) links genetic, lineage, and other cellular features of cancer cell lines to small-molecule sensitivity with the goal of accelerating discovery of patient-matched cancer therapeutics. https://portals.broadinstitute.org/ctrp/


## Consortia

- `PCAGW` - The PCAWG study is an international collaboration to identify common patterns of mutation in more than 2,800 cancer whole genomes from the International Cancer Genome Consortium. The project produced large amount data with many types including simple somatic mutations (SNVs, MNVs and small INDELs), large-scale somatic structural variations, copy number alterations, germline variations, RNA expression profiles, gene fusions, and phenotypic annotations etc. PCAWG data have been imported, processed and made available in the following four major online resources for download and exploration by the cancer researchers worldwide. http://docs.icgc.org/pcawg/
    - Goldman, Mary, Junjun Zhang, Nuno A. Fonseca, Qian Xiang, Brian Craft, Elena Piñeiro, Brian O’Connor, et al. “Online Resources for PCAWG Data Exploration, Visualization, and Discovery.” BioRxiv, October 18, 2017. https://doi.org/10.1101/163907. https://www.biorxiv.org/content/early/2017/10/18/163907

- Data from PanCancer publications. Clinical annotations, RNA-seq counts, RPPA, Methylation, miRNA, copy number, mutations in .maf format. https://gdc.cancer.gov/about-data/publications/pancanatlas

- Zehir, Ahmet, Ryma Benayed, Ronak H Shah, Aijazuddin Syed, Sumit Middha, Hyunjae R Kim, Preethi Srinivasan, et al. “Mutational Landscape of Metastatic Cancer Revealed from Prospective Clinical Sequencing of 10,000 Patients.” Nature Medicine 23, no. 6 (May 8, 2017): 703–13. https://doi.org/10.1038/nm.4333. - MSK-IMPACT study. Deep sequencing of 341-410 genes in 10,000 samples in multiple cancers. Focus on mutations, copy number alterations, fusions. Data at http://www.cbioportal.org/study?id=msk_impact_2017#summary, downloadable, includes clinical data for survival analysis.


## Tools

- `TNBCtype` tool to classify triple negative breast cancer samples (microarray gene expression) into six subtypes, http://cbc.mc.vanderbilt.edu/tnbc/index.php

- `genefu` R package for PAM50 classification and survival analysis. https://www.bioconductor.org/packages/release/bioc/html/genefu.html

- NCI Genomics Data Commons API. https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/ - docs. https://github.com/Bioconductor/GenomicDataCommons - R package
    - Shane Wilson, Michael Fitzsimons, Martin Ferguson, Allison Heath, Mark Jensen, Josh Miller, Mark W. Murphy, James Porter, Himanso Sahni, Louis Staudt, Yajing Tang, Zhining Wang, Christine Yu, Junjun Zhang, Vincent Ferretti and Robert L. Grossman. "Developing Cancer Informatics Applications and Tools Using the NCI Genomic Data Commons API." DOI: 10.1158/0008-5472.CAN-17-0598 Published November 2017 http://cancerres.aacrjournals.org/content/77/21/e15


## Misc

- `MEXPRESS` - Gene-centric methylation and correlation with clinical parameters. http://mexpress.be/

- Zhang, Zhuo, Hao Li, Shuai Jiang, Ruijiang Li, Wanying Li, Hebing Chen, and Xiaochen Bo. “A Survey and Evaluation of Web-Based Tools/Databases for Variant Analysis of TCGA Data.” Briefings in Bioinformatics, March 29, 2018. https://doi.org/10.1093/bib/bby023. - The most comprehensive review of TCGA-related tools. Table 3 - List of Web servers and databases. 

- The Pan-Cancer analysis by TCGA consortium, all papers. https://www.cell.com/pb-assets/consortium/pancanceratlas/pancani3/index.html

- TCGA mutation data, MAF files, open and controlled access. https://gdc.cancer.gov/about-data/publications/mc3-2017



