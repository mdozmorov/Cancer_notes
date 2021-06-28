# A continually expanding collection of cancer genomics tools and data

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Tools](#tools)
  - [Preprocessing](#preprocessing)
  - [Purity](#purity)
  - [Immune cell deconvolution](#immune-cell-deconvolution)
  - [BRCA](#brca)
  - [OvCa](#ovca)
  - [SCLC](#sclc)
  - [TCGA](#tcga)
  - [Integrative](#integrative)
- [Deep Learning](#deep-learning)
  - [Image analysis](#image-analysis)
- [Clonal analysis](#clonal-analysis)
- [Survival analysis](#survival-analysis)
  - [Methods to find best cutoff for survival](#methods-to-find-best-cutoff-for-survival)
- [Cancer driver genes](#cancer-driver-genes)
- [Cancer driver mutations](#cancer-driver-mutations)
- [Drugs](#drugs)
- [Data](#data)
  - [TCGA PanCancer](#tcga-pancancer)
  - [Pediatric](#pediatric)
- [Methylation](#methylation)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Tools

### Preprocessing

- [The Trinity Cancer Transcriptome Analysis Toolkit (CTAT)](https://github.com/NCIP/Trinity_CTAT/wiki) aims to provide tools for leveraging RNA-Seq to gain insights into the biology of cancer transcriptomes. Bioinformatics tool support is provided for mutation detection, fusion transcript identification, de novo transcript assembly of cancer-specific transcripts, lncRNA classification, and foreign transcript detection (viruses, microbes). 

- [cacao](https://github.com/sigven/cacao) - Callable Cancer Loci - assessment of sequencing coverage for actionable and pathogenic loci in cancer, example of QC report. Data: BED files for cancer loci from ClinVar, CIViC, cancerhotspots. 

### Purity

- `ABSOLUTE` - infers tumor purity, ploidy from SNPs, CNVs. Also detects subclonal heterogeneity.
    - Carter, Scott L., Kristian Cibulskis, Elena Helman, Aaron McKenna, Hui Shen, Travis Zack, Peter W. Laird, et al. “Absolute Quantification of Somatic DNA Alterations in Human Cancer.” Nature Biotechnology 30, no. 5 (May 2012): 413–21. https://doi.org/10.1038/nbt.2203.
    - Aran, Dvir, Marina Sirota, and Atul J. Butte. “Systematic Pan-Cancer Analysis of Tumour Purity.” Nature Communications 6, no. 1 (December 2015). https://doi.org/10.1038/ncomms9971. - TCGA tumor purity estimation using four methods: ESTIMATE, ABSOLUTE, LUMP, IHC, and a median consensus purity estimation. Gene expression correlates with purity and may affect correlation and differential expression detection analyses - big confounding effect.
        - [data/ABSOLUTE_scores.xlsx](data/ABSOLUTE_scores.xlsx) - Supplementary Data 1: Tumor purity estimates according to four methods and the consensus method for all TCGA samples with available data. [Source](https://media.nature.com/original/nature-assets/ncomms/2015/151204/ncomms9971/extref/ncomms9971-s2.xlsx)

- `ESTIMATE` (Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data) is a tool for predicting tumor purity, and the presence of infiltrating stromal/immune cells in tumor tissues using gene expression data. ESTIMATE algorithm is based on single sample Gene Set Enrichment Analysis and generates three scores: stromal score (that captures the presence of stroma in tumor tissue), immune score (that represents the infiltration of immune cells in tumor tissue), and estimate score (that infers tumor purity). http://bioinformatics.mdanderson.org/main/ESTIMATE:Overview. R package http://bioinformatics.mdanderson.org/estimate/rpackage.html
    - Yoshihara, Kosuke, Maria Shahmoradgoli, Emmanuel Martínez, Rahulsimham Vegesna, Hoon Kim, Wandaliz Torres-Garcia, Victor Treviño, et al. “Inferring Tumour Purity and Stromal and Immune Cell Admixture from Expression Data.” Nature Communications 4 (2013): 2612. https://doi.org/10.1038/ncomms3612. - ESTIMATE - tumor-stroma purity detection. 141 immune and stromal genes. single-sample GSEA analysis. ESTIMATE score as a combination of immune and stromal scores. [Supplementary data](https://www.nature.com/articles/ncomms3612#supplementary-information).
- `data/ESTIMATE_signatures.xlsx` - A gene list of stromal and immune signatures. [Source](https://media.nature.com/original/nature-assets/ncomms/2013/131011/ncomms3612/extref/ncomms3612-s2.xlsx)
- `data/ESTIMATE_scores.xlsx` - A list of stromal, immune, and ESTIMATE scores in TCGA data sets. All cancers, all gene expression plaforms. [Source](https://media.nature.com/original/nature-assets/ncomms/2013/131011/ncomms3612/extref/ncomms3612-s3.xlsx)

- [ISOpureR](https://CRAN.R-project.org/package=ISOpureR) - Deconvolution of Tumour Profiles to purify tumor samples. Regression-based, uses purified tumor profile to estimate the proportion of tumor samples. Discussion of overfitting due to overparametrization
    - Quon, Gerald, Syed Haider, Amit G Deshwar, Ang Cui, Paul C Boutros, and Quaid Morris. “[Computational Purification of Individual Tumor Gene Expression Profiles Leads to Significant Improvements in Prognostic Prediction](https://doi.org/10.1186/gm433).” Genome Medicine 5, no. 3 (2013)

### Immune cell deconvolution

- [Tumor Immune Single-cell Hub (TISCH)](http://tisch.comp-genomics.org/home/) is a scRNA-seq database focusing on tumor microenvironment (TME). TISCH provides detailed cell-type annotation at the single-cell level, enabling the exploration of TME across different cancer types. [Tweet](https://twitter.com/XShirleyLiu/status/1408472207152533506?s=20)

- [quanTIseq](https://icbi.i-med.ac.at/software/quantiseq/doc/) - quantification of the Tumor Immune cell proportions from human RNA-seq data. Input - FASTQ files (Trimmomatic, Kallisto to TPM, normalization), or TPM matrix. Deconvolution into 10 cell types (B, macrophages M1 and M2, Monocytes, Neutrophils, NK, CD8 T, CD4 T, Dendritic cells), and uncharacterized fraction (TIL10 signature). Custom processing of 51 datasets to generate TIL10. Compared with CIBERSORT, TIMER, EPIC on simulated and real-life data. Command line, Docker/Singularity implementation, no GitHub.
    - Finotello, Francesca, Clemens Mayer, Christina Plattner, Gerhard Laschober, Dietmar Rieder, Hubert Hackl, Anne Krogsdam, et al. “[Molecular and Pharmacological Modulators of the Tumor Immune Contexture Revealed by Deconvolution of RNA-Seq Data](https://doi.org/10.1186/s13073-019-0638-6).” Genome Medicine 11, no. 1 (May 24, 2019)

- [CIBERSORT](https://cibersort.stanford.edu/) - cell type identification using Support Vector Regression. p-value for the overall goodness of deconvolution (H0 - no cell types are present in a given gene expression profile), also Pearson and RMSE for estimating goodness of fit. References to six GEP deconvolution methods: linear least-squares regression (LLSR), quadratic programming (QP), perturbation model for gene expression deconvolution (PERT), robust linear regression (RLR), microarray microdissection with analysis of differences (MMAD) and digital sorting algorithm (DSA). References to datasets for benchmarking.
    - Newman, Aaron M., Chih Long Liu, Michael R. Green, Andrew J. Gentles, Weiguo Feng, Yue Xu, Chuong D. Hoang, Maximilian Diehn, and Ash A. Alizadeh. “[Robust Enumeration of Cell Subsets from Tissue Expression Profiles](https://doi.org/10.1038/nmeth.3337).” Nature Methods 12, no. 5 (May 2015)

- [DeconRNAseq](http://bioconductor.org/packages/DeconRNASeq/) - deconvolution of RNA-seq datasets into cell proportions using cell signatures. Non-negative decomposition algorithm (X = AS) solved using quadratic programming
    - Gong, Ting, and Joseph D. Szustakowski. “[DeconRNASeq: A Statistical Framework for Deconvolution of Heterogeneous Tissue Samples Based on MRNA-Seq Data](https://doi.org/10.1093/bioinformatics/btt090).” Bioinformatics (Oxford, England) 29, no. 8 (April 15, 2013)

- `Immunophenogram` - partitioning immune cell types in cancer. https://github.com/mui-icbi/Immunophenogram, https://tcia.at/home
    - Charoentong, Pornpimol, Francesca Finotello, Mihaela Angelova, Clemens Mayer, Mirjana Efremova, Dietmar Rieder, Hubert Hackl, and Zlatko Trajanoski. “Pan-Cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade.” BioRxiv, 2016, 056101. - https://tcia.at/ - immune cells in cancers. Estimated using functional GSEA enrichment, and CIBERSORT. Immunophenogram generation: https://github.com/MayerC-imed/Immunophenogram

- `ImmQuant` - Deconvolution of immune cell lineages. http://csgi.tau.ac.il/ImmQuant/downloads.html
    - Frishberg, Amit, Avital Brodt, Yael Steuerman, and Irit Gat-Viks. “ImmQuant: A User-Friendly Tool for Inferring Immune Cell-Type Composition from Gene-Expression Data.” Bioinformatics 32, no. 24 (December 15, 2016): 3842–43. https://doi.org/10.1093/bioinformatics/btw535.

- `TIMER` - a resource to estimate the abundance of immune infiltration of six cell types, the effect on survival. Accounting for tumor purity using CHAT R package. Linear regression to estimate immune cell abundance. Macrophage infiltration predicts worse outcome, including BRCA. Signature genes from microarray, merged with RNA-seq using ComBat for batch removal, filtered. All TCGA processed. All data at http://cistrome.org/TIMER/download.html, tool at https://cistrome.shinyapps.io/timer/
    - Li, Bo, Eric Severson, Jean-Christophe Pignon, Haoquan Zhao, Taiwen Li, Jesse Novak, Peng Jiang, et al. “Comprehensive Analyses of Tumor Immunity: Implications for Cancer Immunotherapy.” Genome Biology 17, no. 1 (22 2016): 174. https://doi.org/10.1186/s13059-016-1028-7.

- [TIDE](http://tide.dfci.harvard.edu/) - Tumor Immune Dysfunction and Exclusion, a gene expression biomarker to predict the clinical response to immune checkpoint blockade using patient-specific gene expression. http://tide.dfci.harvard.edu/

- [TRUST4](https://github.com/liulab-dfci/TRUST4) - infers tumor-infiltrating TCR and BCR repertores from bulk and scRNA-seq (5' 10X Genomics kit). https://github.com/liulab-dfci/TRUST4

### BRCA

- [HRDetect](https://github.com/eyzhao/hrdetect-pipeline) (homologous recombination-repair deficiency) classification of TNBC patients. Whole-genome sequencing-based. HRDetect-high have higher chemosensitivity, better survival profiles
    - Staaf, Johan, Dominik Glodzik, Ana Bosch, Johan Vallon-Christersson, Christel Reuterswärd, Jari Häkkinen, Andrea Degasperi, et al. “[Whole-Genome Sequencing of Triple-Negative Breast Cancers in a Population-Based Clinical Study](https://doi.org/10.1038/s41591-019-0582-4).” Nature Medicine, September 30, 2019

- `TNBCtype` tool to classify triple negative breast cancer samples (microarray gene expression) into six subtypes, http://cbc.mc.vanderbilt.edu/tnbc/index.php

- `genefu` R package for PAM50 classification and survival analysis. https://www.bioconductor.org/packages/release/bioc/html/genefu.html

### OvCa

- [PrOTYPE](https://ovcare.shinyapps.io/PrOType/) - Ovarian Cancer subtype prediction. Model trained on gene expression from 1650 tumors (specimens from the Ovarian Tumor Tissue Analysis consortium), validated on NanoString data on 3829 tumors. 55 genes, predict with >95% accuracy, also associated with age, stage, residual disease, TILs, outcome. Review of previous studies classifying into 5 outcomes. Random Forest, cross-validation. [Supplemental Table SC7](https://clincancerres.aacrjournals.org/content/early/2020/06/17/1078-0432.CCR-20-0103.figures-only) lists all predictor genes. PrOTYPE web tool to classify NanoString OvCa samples, https://ovcare.shinyapps.io/PrOType/
    - Talhouk, Aline, Joshy George, Chen Wang, Timothy Budden, Tuan Zea Tan, Derek S Chiu, Stefan Kommoss, et al. “Development and Validation of the Gene-Expression Predictor of High-Grade-Serous Ovarian Carcinoma Molecular SubTYPE (PrOTYPE).” Clinical Cancer Research, June 17, 2020, clincanres.0103.2020. https://doi.org/10.1158/1078-0432.CCR-20-0103.

### SCLC

- [SCLC-CellMiner](https://discover.nci.nih.gov/SclcCellMinerCDB/) - an analysis portal for Small Cell Lung Cancer. Integrating drug sensitivity and methylome (Illumina 80K), CNV (ChAMP, from methylation), and transcriptome data from 118 SCLC cell lines. Integrated data from multiple resources, uniformly processed. Confirms NEUROD1, ASCL1, POU2F3, and YAP1 (NAPY) classification. Table 1 - types of analyses. [Shiny web interface](https://discover.nci.nih.gov/SclcCellMinerCDB/) and [data download](https://zenodo.org/record/3959142#.X9-lf2RKgoI)
    - Tlemsani, Camille, Lorinc Pongor, Fathi Elloumi, Luc Girard, Kenneth E. Huffman, Nitin Roper, Sudhir Varma, et al. “[SCLC-CellMiner: A Resource for Small Cell Lung Cancer Cell Line Genomics and Pharmacology Based on Genomic Signatures](https://doi.org/10.1016/j.celrep.2020.108296).” Cell Reports 33, no. 3 (October 2020)



### TCGA

- Google Cloud Pilot RNA-Sequencing for CCLE and TCGA, https://osf.io/gqrz9/

- PanCancer atlas RNA-seq, RPPA, Methylation, miRNA, Copy Number, Mutation, Clinical data, includes ABSOLUTE purity estimates, https://gdc.cancer.gov/about-data/publications/pancanatlas

- Normalized TCGA RNA-seq data, from Yu, K., Chen, B., Aran, D., Charalel, J., Yau, C., Wolf, D.M., van ‘t Veer, L.J., Butte, A.J., Goldstein, T., and Sirota, M. (2019). Comprehensive transcriptomic analysis of cell lines as models of primary tumors across 22 tumor types. Nature Communications 10. https://www.synapse.org/#!Synapse:syn18685536/files/

- Zhang, Zhuo, Hao Li, Shuai Jiang, Ruijiang Li, Wanying Li, Hebing Chen, and Xiaochen Bo. “A Survey and Evaluation of Web-Based Tools/Databases for Variant Analysis of TCGA Data.” Briefings in Bioinformatics, March 29, 2018. https://doi.org/10.1093/bib/bby023. - The most comprehensive review of TCGA-related tools. Table 3 - List of Web servers and databases. 

- NCI Genomics Data Commons API. https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/ - docs. https://github.com/Bioconductor/GenomicDataCommons - R package
    - Shane Wilson, Michael Fitzsimons, Martin Ferguson, Allison Heath, Mark Jensen, Josh Miller, Mark W. Murphy, James Porter, Himanso Sahni, Louis Staudt, Yajing Tang, Zhining Wang, Christine Yu, Junjun Zhang, Vincent Ferretti and Robert L. Grossman. "Developing Cancer Informatics Applications and Tools Using the NCI Genomic Data Commons API." DOI: 10.1158/0008-5472.CAN-17-0598 Published November 2017 http://cancerres.aacrjournals.org/content/77/21/e15

### Integrative

- Review of tools and methods for the integrative analysis of multiple omics data, cancer-oriented. [Table 1](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/bbia/2020/bbia_14/1177932219899051/20200130/images/large/10.1177_1177932219899051-table1.jpeg) - multi-omics data repositories (TCGA, CPTAC, ICGC, CCLE, METABRIC, TARGET, Omics Discovery Index). Three broad areas of multi-omics analysis: 1. Disease subtyping and classification based on multi-omics profiles; 2. Prediction of biomarkers for various applications including diagnostics and driver genes for diseases; 3. Deriving insights into disease biology. [Table 2](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/bbia/2020/bbia_14/1177932219899051/20200130/images/large/10.1177_1177932219899051-table2.jpeg) - software categorized by use case (PARADIGM, iClusterPlus, PSDF, BCC, MDI, SNF, PFA, PINSPlus, NEMO, mixOmics, moCluster, MCIA, JIVE, MFA, sMBPLS, T-SVD, Joint NMF). Brief description of each tool, links, exemplary publications. [Table 3](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/bbia/2020/bbia_14/1177932219899051/20200130/images/large/10.1177_1177932219899051-table3.jpeg) - visualization portals (cBioPortal, Firebrowse, UCSC Xena, LinkedOmics, 3Omics, NetGestalt, OASIS, Paintomics, MethHC). Description of each, data types, analysis examples.
    - Subramanian, Indhupriya, Srikant Verma, Shiva Kumar, Abhay Jere, and Krishanpal Anamika. “[Multi-Omics Data Integration, Interpretation, and Its Application](https://doi.org/10.1177/1177932219899051).” Bioinformatics and Biology Insights, January 31, 2020 

- [LinkedOmics](http://www.linkedomics.org/login.php) - clinical, gennomic (expression, SNPs, CNVs, methylation, miRNAs) and protein expression data from TCGA, CPTAC. Three analysis modules: LinkFinder finds associations between molecular and clinical attributes; LinkCompare compares the associations. LinkInterpreter maps associations to pathways and networks.
    - Vasaikar, Suhas V, Peter Straub, Jing Wang, and Bing Zhang. “[LinkedOmics: Analyzing Multi-Omics Data within and across 32 Cancer Types](https://doi.org/10.1093/nar/gkx1090).” Nucleic Acids Research 46, no. D1 (January 4, 2018)


## Deep Learning

- Review of drug response prediction methods in cancer. Overview of deep learning architectures, Table 1 - DL terminology, Table 2 - data resources (NCI-60, CCLE, GDSC1 etc.), Table 3 - drug combination screening studies. DREAM challenges addressing drug sensitivity. Approaches for building and evaluating drug response prediction models, single-drug (Table 4 - studies, models) and drug combination (Table 6) approaches. DeepSynergy, NCI-ALMANAC dataset, AstraZeneca-Sanger Drug Combination DREAM challenge paper. Evaluation of model ensembles. Deep learning for drug repurposing. Methods for model evaluation, improving performance, increase interpretability. Deep reading with lots of references.
    - Baptista, Delora, Pedro G Ferreira, and Miguel Rocha. “[Deep Learning for Drug Response Prediction in Cancer](https://doi.org/10.1093/bib/bbz171),” Briefings in Bioinformatics, 17 January 2020

- [DrugCell](https://github.com/idekerlab/DrugCell) - interpretable neural network, directly mapping neurons of a deep neural network to Gene Ontology hierarchy. Input - genotype (binary gene indicator) and drugs, output - response of cell to drug. Performance is the same as unconstrained network. Trained on 1235 tumor cell lines and 684 drugs. Allows for the detection of synergistic combinations. [Perspective](https://www.sciencedirect.com/science/article/pii/S1535610820305456)
    - Kuenzi, Brent M., Jisoo Park, Samson H. Fong, Kyle S. Sanchez, John Lee, Jason F. Kreisberg, Jianzhu Ma, and Trey Ideker. “[Predicting Drug Response and Synergy Using a Deep Learning Model of Human Cancer Cells](https://doi.org/10.1016/j.ccell.2020.09.014).” Cancer Cell 38, no. 5 (November 2020)

- Drug response prediction from gene expression data. Deep Neural Network (DNN, H2O.ai framework) compared with Elastic Net, Random Forest. Trained on highly variable (by MAD) gene expression in 1001 cell lines and 251 drugs pharmacogenomic dataset (GDSC. CCLP) to predict IC50. Hyper-parameter optimization using 5-fold cross-validation and minimizing Mean Square Error. Batch correction between the datasets Tested on unseen patient cohorts (OCCAMS, MD Anderson, TCGA, Multiple Myeloma Consortium) to predict IC50 and test low, medium, high IC50 groups for survival differences. [RDS files data](https://genome.med.nyu.edu/public/tsirigoslab/deep-drug-response/), [R code](https://genome.med.nyu.edu/public/tsirigoslab/deep-drug-response/)
    - Sakellaropoulos, Theodore, Konstantinos Vougas, Sonali Narang, Filippos Koinis, Athanassios Kotsinas, Alexander Polyzos, Tyler J. Moss, et al. “[A Deep Learning Framework for Predicting Response to Therapy in Cancer](https://doi.org/10.1016/j.celrep.2019.11.017).” Cell Reports 29, no. 11 (December 2019)


### Image analysis

- [High content image analysis with CellProfiler (2D and 3D)](https://github.com/hmbotelho/SPAOM2020-ws1-high-content-screening) - detailed Python- and R-based workshop

- Prediction of chromosomal instability (CIN, fraction of genome altered, binarized into high/low) in breast cancer using deep learning on histology images. Tested different convnet architectures, Densenet-121 worked best. TensorFlow2. [Code for CIN](https://github.com/eipm/CIN)
    - Xu, Zhuoran, Akanksha Verma, Uska Naveed, Samuel Bakhoum, Pegah Khosravi, and Olivier Elemento. “[Using Histopathology Images to Predict Chromosomal Instability in Breast Cancer: A Deep Learning Approach](https://doi.org/10.1101/2020.09.23.20200139).” Preprint. Genetic and Genomic Medicine, September 24, 2020. 

- `DeepPATH` - Lung cancer image classification using deep convolutional neural network. Classification by tumor type, mutation type. Refs to other image classification studies that use deep learning. GoogleNet inception v3 architecture. Training, validation, testing cohorts (70%, 15%, 15%). Details on image processing. https://github.com/ncoudray/DeepPATH
    - Coudray, Nicolas, Paolo Santiago Ocampo, Theodore Sakellaropoulos, Navneet Narula, Matija Snuderl, David Fenyö, Andre L. Moreira, Narges Razavian, and Aristotelis Tsirigos. “Classification and Mutation Prediction from Non–Small Cell Lung Cancer Histopathology Images Using Deep Learning.” Nature Medicine 24, no. 10 (October 2018): 1559–67. https://doi.org/10.1038/s41591-018-0177-5.


- `IHCount` - IHC-image analysis workflow, https://github.com/mui-icbi/IHCount

- `pathology_learning` - Using traditional machine learning and deep learning methods to predict stuff from TCGA pathology slides. [https://github.com/millett/pathology_learning](https://github.com/millett/pathology_learning)

## Clonal analysis

- `Awesome-CancerEvolution` - list of papers and tools for studying cancer evolution. https://github.com/iron-lion/Awesome-CancerEvolution

- [timescape](https://bioconductor.org/packages/timescape/) - an R package for visualizing temporal clonal evolution data. [Examples](https://bioconductor.org/packages/release/bioc/vignettes/timescape/inst/doc/timescape_vignette.html), [GitHub](https://github.com/shahcompbio/timescape)

- `clonevol` R package, Inferring and visualizing clonal evolution in multi-sample cancer sequencing. https://github.com/hdng/clonevol
    - Dang, H. X., B. S. White, S. M. Foltz, C. A. Miller, J. Luo, R. C. Fields, and C. A. Maher. “ClonEvol: Clonal Ordering and Visualization in Cancer Sequencing.” Annals of Oncology: Official Journal of the European Society for Medical Oncology 28, no. 12 (December 1, 2017): 3076–82. https://doi.org/10.1093/annonc/mdx517.

- `ape` - R package, Analyses of Phylogenetics and Evolution, https://cran.r-project.org/web/packages/ape/index.html

- `DeconstructSig` - Contribution of known SNP cancer mutation signatures to tumor samples. Data from Alexandrov, COSMIC, others. https://github.com/raerose01/deconstructSigs 

- `E-scape` - cancer evolution visualization. Timescape - time series analysis, http://bioconductor.org/packages/release/bioc/html/timescape.html, MapScape - spatial distribution,http://bioconductor.org/packages/release/bioc/html/mapscape.html, CellScape - single-cell phylogenetic, http://bioconductor.org/packages/release/bioc/html/cellscape.html
    - Smith, Maia A., Cydney B. Nielsen, Fong Chun Chan, Andrew McPherson, Andrew Roth, Hossein Farahani, Daniel Machev, Adi Steif, and Sohrab P. Shah. “E-Scape: Interactive Visualization of Single-Cell Phylogenetics and Cancer Evolution.” Nature Methods 14, no. 6 (30 2017): 549–50. https://doi.org/10.1038/nmeth.4303.

- `fishplot` - Create timecourse "fish plots" that show changes in the clonal architecture of tumors. https://github.com/chrisamiller/fishplot

- `MACHINA` - Metastatic And Clonal History INtegrative Analysis. https://github.com/raphael-group/machina 

- `MEDICC` - Minimum Event Distance for Intra-tumour Copy number Comparisons. https://bitbucket.org/rfs/medicc

- `oncoNEM` - A package for inferring clonal lineage trees from single-cell somatic single nucleotide variants. https://bitbucket.org/edith_ross/onconem

- `PyClone` - inferring the cellular prevalence of point mutations from deeply sequenced data. http://compbio.bccrc.ca/software/pyclone/

- `SciClone` - number and genetic composition of tumor subclones by analyzing the variant allele frequencies of somatic mutations. Excludes CNV regions. https://github.com/genome/sciclone

- `treeomics` - Decrypting somatic mutation patterns to reveal the evolution of cancer. https://github.com/johannesreiter/treeomics


## Survival analysis

- [www.tcgaportal.org](http://tcgaportal.org/index.html) - web server for survival analysis using TCGA data

- `cBioPortal` - The cBioPortal for Cancer Genomics provides visualization, analysis and download of large-scale cancer genomics data sets. OncoPrint mutation plots, differential expression, coexpression, survival. Compare gene expression with copy number variation. http://www.cbioportal.org/

- `R2` - Genomics Analysis and Visualization Platform. Gene-centric, survival analysis, collection of preprocessed microarray studies. http://hgserver1.amc.nl/

- `KM plotter` - Gene-centric, customizable survival analysis for breast, ovarian, lung, gastric cancers. http://kmplot.com/
    - Györffy, Balazs, Andras Lanczky, Aron C. Eklund, Carsten Denkert, Jan Budczies, Qiyuan Li, and Zoltan Szallasi. “An Online Survival Analysis Tool to Rapidly Assess the Effect of 22,277 Genes on Breast Cancer Prognosis Using Microarray Data of 1,809 Patients.” Breast Cancer Research and Treatment 123, no. 3 (October 2010): 725–31. https://doi.org/10.1007/s10549-009-0674-9.

- `The Human Protein Atlas: Pathology atlas` - Gene- and protein expression data in multiple cancer tissues, cell lines. Easy one-gene search, summary of tissue-specific expression, survival significance. http://www.proteinatlas.org/
    - Uhlen, Mathias, Cheng Zhang, Sunjae Lee, Evelina Sjöstedt, Linn Fagerberg, Gholamreza Bidkhori, Rui Benfeitas, et al. “[A Pathology Atlas of the Human Cancer Transcriptome](http://science.sciencemag.org/content/357/6352/eaan2507).” Science (August 18, 2017). [Data download](http://www.proteinatlas.org/about/download) - tissue-specific gene expression in cancer and normal, isoform expression, protein expression. [Supplementary material](http://science.sciencemag.org/content/suppl/2017/08/16/357.6352.eaan2507.DC1)  
        - [Table S2](https://science.sciencemag.org/highwire/filestream/698233/field_highwire_adjunct_files/0/Supplementary-Tables.zip) - summary of tissue specific expression for each gene, in normal and cancer tissues.
        - [Table S6](https://science.sciencemag.org/highwire/filestream/698233/field_highwire_adjunct_files/0/Supplementary-Tables.zip) - summary of survival prognostic value, with a simple "favorable/unfavorable" label for each gene. Each worksheet corresponds to a different cancer.  
        - [Table S8](https://science.sciencemag.org/highwire/filestream/698233/field_highwire_adjunct_files/0/Supplementary-Tables.zip) - per-gene summary, in which cancers it is prognostic of survival.  

- [Breast Cancer Gene-Expression Miner v4.4](http://bcgenex.centregauducheau.fr/BC-GEM/GEM-Accueil.php?js=1) - gene expression, correlation, and survival analysis in different microarray (e.g., METABRIC) and RNA-seq (e.g., TCGA) datasets

- [G-2-O, Genotype to Outcome](http://www.g-2-o.com/) -  web-server linking mutation (or CNV) of a gene to clinical outcome (survival) by utilizing next generation sequencing and gene chip data. For Breast and Lung cancer

- `PRECOG` - PREdiction of Clinical Outcomes from Genomic Profiles. Gene-centric, quick overview of survival effect of a gene across all cancers, KM plots. https://precog.stanford.edu
    - Gentles, Andrew J., Aaron M. Newman, Chih Long Liu, Scott V. Bratman, Weiguo Feng, Dongkyoon Kim, Viswam S. Nair, et al. “The Prognostic Landscape of Genes and Infiltrating Immune Cells across Human Cancers.” Nature Medicine 21, no. 8 (August 2015): 938–45. https://doi.org/10.1038/nm.3909. - TCGA pan-cancer survival analysis PRECOG, CIBERSORT. 39 cancers. Intro into heterogeneity. Z-score description. Batch effect does not significantly affect z-scores. 2/3 prognostic genes shared across cancers. AutoSOME clustering method

- `GEPIA` - single- and multiple-gene analyses of TCGA data. Gene expression in different tumor-normal comparisons, differentially expressed genes, correlation analysis, similar genes, survival analysis. http://gepia.cancer-pku.cn/
    - Zefang Tang et al., “GEPIA: A Web Server for Cancer and Normal Gene Expression Profiling and Interactive Analyses,” Nucleic Acids Research 45, no. W1 (July 3, 2017): W98–102, https://doi.org/10.1093/nar/gkx247. - TCGA and GTEX web interface. Classical analyses - differential expression analysis, profiling plotting, correlation analysis, patient survival analysis, similar gene detection and dimensionality reduction analysis. http://gepia.cancer-pku.cn/
- `GEPIA2` - isoform-level TCGA analysis. Cancer subtype-specific analyses. Eight types of expression analyses, and additional Cancer Subtype Classifier and Expression Comparison. Python package for API access. http://gepia2.cancer-pku.cn
    - Tang, Zefang, Boxi Kang, Chenwei Li, Tianxiang Chen, and Zemin Zhang. “GEPIA2: An Enhanced Web Server for Large-Scale Expression Profiling and Interactive Analysis.” Nucleic Acids Research, May 22, 2019. https://doi.org/10.1093/nar/gkz430.



<!--
- `PrognoScan`, Gene-centric, survival effect of a gene in cancer studies from GEO. http://dna00.bio.kyutech.ac.jp/PrognoScan/
-->

- `UALCAN` - Gene-centric, tumor-normal expression, survival analusis, TCGA cancers. http://ualcan.path.uab.edu/
    - Chandrashekar DS, Bashel B, Balasubramanya SAH, Creighton CJ, Rodriguez IP, Chakravarthi BVSK and Varambally S. UALCAN: A portal for facilitating tumor subgroup gene expression and survival analyses. Neoplasia. 2017 Aug;19(8):649-658. doi: 10.1016/j.neo.2017.05.002 [PMID:28732212]

- `Project Betastasis` - Gene-centric, survival analysis, gene expression, select cancer studies. http://www.betastasis.com/

- `OncoLnc` - Gene-centric, survival analysis in any TCGA cancer. http://www.oncolnc.org/

### Methods to find best cutoff for survival

- `KMplotter` - the Kaplan Meier plotter is capable to assess the effect of 54,675 genes on survival using 18,674 cancer samples. These include 5,143 breast, 1,816 ovarian, 2,437 lung, 364 liver, 1,065 gastric cancer patients with relapse-free and overall survival data. The miRNA subsystems include additional 11,456 samples from 20 different cancer types. Primary purpose of the tool is a meta-analysis based biomarker assessment.
    - Györffy, Balazs, Andras Lanczky, Aron C. Eklund, Carsten Denkert, Jan Budczies, Qiyuan Li, and Zoltan Szallasi. “An Online Survival Analysis Tool to Rapidly Assess the Effect of 22,277 Genes on Breast Cancer Prognosis Using Microarray Data of 1,809 Patients.” Breast Cancer Research and Treatment 123, no. 3 (October 2010): 725–31. https://doi.org/10.1007/s10549-009-0674-9. - cutoff selection for survival by scanning gene expression range.

- `ctree` function for automatic cutoff finding and building a regression tree out of multiple covariates. `partykit::ctree()`. 
    - Hothorn, Torsten, Kurt Hornik, and Achim Zeileis. “Ctree: Conditional Inference Trees.” The Comprehensive R Archive Network, 2015, 1–34.

- `Cutoff Finder` - web tool for finding optimal dichotomization with respect to an outcome or survival variable. Five methods. http://molpath.charite.de/cutoff/
    - Budczies, Jan, Frederick Klauschen, Bruno V. Sinn, Balázs Győrffy, Wolfgang D. Schmitt, Silvia Darb-Esfahani, and Carsten Denkert. “Cutoff Finder: A Comprehensive and Straightforward Web Application Enabling Rapid Biomarker Cutoff Optimization.” PloS One 7, no. 12 (2012): e51862. https://doi.org/10.1371/journal.pone.0051862.


## Cancer driver genes

- [The list of cancer-related genes from the Bushman Lab](http://www.bushmanlab.org/assets/doc/allOnco_May2018.tsv). [allOnco_May2018.tsv](data/allOnco_May2018.tsv)

- Integrative pathway and network analysis of 2583 cancers (27 tumor types) identified 87 driver Pathway Implicated Driver (PID) genes with coding variants (PID-C) and 93 drivers with noncoding variants (PID-N). These gene classes are associated with different biological processes. Six pathway databases, seven pathway and network methods, data references in Methods. Non-Coding Added Value (NCVA) score to identify genes with noncoding variants increasing the overall significance.
    - PCAWG Drivers and Functional Interpretation Working Group, PCAWG Consortium, Matthew A. Reyna, David Haan, Marta Paczkowska, Lieven P. C. Verbeke, Miguel Vazquez, et al. “[Pathway and Network Analysis of More than 2500 Whole Cancer Genomes](https://doi.org/10.1038/s41467-020-14367-0).” Nature Communications 11, no. 1 (December 2020)
    - [Supplementary Data 2: PID-C genes. List of 87 genes](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14367-0/MediaObjects/41467_2020_14367_MOESM4_ESM.txt), [PCAWG_2020_PID_C_87_genes.txt](data/PCAWG_2020_PID_C_87_genes.txt)
    - [Supplementary Data 3: PID-N genes. List of 93  genes](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14367-0/MediaObjects/41467_2020_14367_MOESM5_ESM.txt), [PCAWG_2020_PID_N_93_genes.txt](data/PCAWG_2020_PID_N_93_genes.txt)
    - More lists in [Supplementary data](https://www.nature.com/articles/s41467-020-14367-0#Sec23), [description](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14367-0/MediaObjects/41467_2020_14367_MOESM2_ESM.pdf). 

- `MoonlightR` - integrative analysis of TCGA data to predict cancer driver genes. http://bioconductor.org/packages/release/bioc/vignettes/MoonlightR/inst/doc/Moonlight.html, https://github.com/ibsquare/MoonlightR
    - [Supplementary Data 5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-13803-0/MediaObjects/41467_2019_13803_MOESM8_ESM.xlsx) - Cancer Driver Genes for TCGA BRCA molecular subtypes
    - [Supplementary Data 6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-13803-0/MediaObjects/41467_2019_13803_MOESM9_ESM.xlsx) - Moonlight’s oncogenic mediators in 18 cancer types
    - Other supplementary data - oncogenic mediators overlapping with methylation, chromatin accessibility, copy number changes, mutations, survival.
    - Colaprico, Antonio, Catharina Olsen, Matthew H. Bailey, Gabriel J. Odom, Thilde Terkelsen, Tiago C. Silva, André V. Olsen, et al. “Interpreting Pathways to Discover Cancer Driver Genes with Moonlight.” Nature Communications 11, no. 1 (December 2020): 69. https://doi.org/10.1038/s41467-019-13803-0.


- `CancerGeneNet` - CancerGeneNet is a resource that aims at linking genes that are frequently mutated in cancers to cancer phenotypes. The resource takes advantage of a curation effort aimed at embedding a large fraction of the gene products that are found altered in cancers in the cell network of causal protein relationships. Graph algorithms, in turn, allow to infer likely paths of causal interactions linking cancer associated genes to cancer phenotypes thus offering a rational framework for the design of strategies to revert disease phenotypes. CancerGenNet bridges two interaction layers by connecting proteins whose activities are affected by cancer gene products to proteins that impact on cancer phenotypes. This is achieved by implementing graph algorithms that allow searching for graph path that link any gene of interest to the “hallmarks of cancer". https://signor.uniroma2.it/CancerGeneNet/


- Cancer Gene Census (CGC), download [COSMIC](http://cancer.sanger.ac.uk/cosmic/download)
    - Hudson, T. J. et al. International network of cancer genome projects. Nature 464, 993–8 (2010).
    - [data/Census_all.csv](data/Census_allThu_Dec_3_19_52_03_2020.csv) - [The cancer Gene Census](http://cancer.sanger.ac.uk/census). Updated 2020-12-03
    - [data/COSMIC_genes.txt](data/COSMIC_genes.txt) - Genes sorted by the number of records associated with them. Obtained using `cat Census_allThu_Dec_3_19_52_03_2020.csv  | sed '1d' | cut -f1 -d, | sort | uniq  > COSMIC_genes.txt`. Updated 2020-12-03
    - `data/CosmicCodingMuts.vcf.gz` - VCF file of all coding mutations in the current release (release v83, 7th November 2017).

- Oncology Model Fidelity Score based on the Hallmarks of Cancer, an R/Shiny app to check cancer samples from preprocessed or user-supplied gene expression data for the presence of these hallmarks. https://github.com/tedgoldstein/hallmarks

- Tumor suppressor gene database (TSGene), https://bioinfo.uth.edu/TSGene/
    - Zhao, M., Sun, J. & Zhao, Z. TSGene: a web resource for tumor suppressor genes. Nucleic Acids Res, 41(Database issue), D970–6 (2013).
    - Download various lists of tumor suppressor genes, https://bioinfo.uth.edu/TSGene/download.cgi

- `OncoScape` - Genes with oncogenic/tumor suppressor/combined scores as a sum contribution from gene expression, somatic mutations, DNA copy-number and methylation as well as data from shRNA knock-down screens. http://oncoscape.nki.nl/
    - Schlicker, Andreas, Magali Michaut, Rubayte Rahman, and Lodewyk F. A. Wessels. “OncoScape: Exploring the Cancer Aberration Landscape by Genomic Data Fusion.” Scientific Reports 6 (20 2016): 28103. https://doi.org/10.1038/srep28103.

- [data/Bailey_2018_cancer_genes.xlsx](Bailey_2018_cancer_genes.xlsx) - Table S1, consensus list of cancer driver genes.
	- Bailey, Matthew H., Collin Tokheim, Eduard Porta-Pardo, Sohini Sengupta, Denis Bertrand, Amila Weerasinghe, Antonio Colaprico, et al. “Comprehensive Characterization of Cancer Driver Genes and Mutations.” Cell 173, no. 2 (April 5, 2018): 371-385.e18. https://doi.org/10.1016/j.cell.2018.02.060. - Pan-Cancer mutation analysis. Combined use of 26 tools (https://www.cell.com/cell/fulltext/S0092-8674(18)30237-X#secsectitle0075, description of each tool in Methods) on harmonized data. 299 cancer driver genes, >3,400 putative missense driver mutations. Table S6 - excluded TCGA samples.

- [data/TARGET_db_v3_02142015.xlsx](data/TARGET_db_v3_02142015.xlsx) - TARGET (tumor alterations relevant for genomics-driven therapy) is a database of genes that, when somatically altered in cancer, are directly linked to a clinical action. TARGET genes may be predictive of response or resistance to a therapy, prognostic, and/or diagnostic. https://software.broadinstitute.org/cancer/cga/target

- [data/Tokheim_2016_cancer_driver_genes.xlsx](data/Tokheim_2016_cancer_driver_genes.xlsx) - Dataset S2: Predicted driver genes by various number of methods
    - Tokheim, Collin J., Nickolas Papadopoulos, Kenneth W. Kinzler, Bert Vogelstein, and Rachel Karchin. “Evaluating the Evaluation of Cancer Driver Genes.” Proceedings of the National Academy of Sciences 113, no. 50 (December 13, 2016): 14330–35. https://doi.org/10.1073/pnas.1616440113. - 20/20+ machine learning method, ratiometric approach to predict cancer driver genes. Performance comparison of other methods, 20/20+, TUSON, OncodriveFML and MutsigCV are the top performers. https://github.com/KarchinLab/2020plus

## Cancer driver mutations

- Resources / databases for clinical interpretation of cancer variants, by Malachi Griffith, https://www.biostars.org/p/403117/

- [sigminer](https://github.com/ShixiangWang/sigminer) - an R package for SNP, CNV, DBS, InDel signature extraction from whole-exome data. NMF-based. Tested on tumor-notmal prostate cancer data. https://github.com/ShixiangWang/sigminer
    - Wang, Shixiang, Huimin Li, Minfang Song, Zaoke He, Tao Wu, Xuan Wang, Ziyu Tao, Kai Wu, and Xue-Song Liu. “Copy Number Signature Analyses in Prostate Cancer Reveal Distinct Etiologies and Clinical Outcomes.” Preprint. Genetic and Genomic Medicine, April 29, 2020. https://doi.org/10.1101/2020.04.27.20082404.

- `clinvar` -  tools to convert ClinVar data into a tab-delimited flat file, and also provides that resulting tab-delimited flat file. https://github.com/macarthur-lab/clinvar

- `CANCERSIGN` - identifies 3-mer and 5-mer mutational signatures, cluster samples by signatures. Based on Alexandrov method, Non-negative matrix factorization, explanation. Other tools - SomaticSignatures, SigneR, deconstructSigs, compared in Table 1. https://github.com/ictic-bioinformatics/CANCERSIGN
    - Bayati, Masroor, Hamid Reza Rabiee, Mehrdad Mehrbod, Fatemeh Vafaee, Diako Ebrahimi, Alistair Forrest, and Hamid Alinejad-Rokny. “CANCERSIGN: A User-Friendly and Robust Tool for Identification and Classification of Mutational Signatures and Patterns in Cancer Genomes.” BioRxiv, January 1, 2019, 424960. https://doi.org/10.1101/424960.


## Drugs

- Drug combination screen, synergy. Statistics for the analysis of large-scale drug screens, 108 drugs, 40 cell lines. Bliss independence model description. Bliss-based linear model to evaluate viabilities for individual drugs. Web-interface: http://www.cmtlab.org:3000/combo_app.html. Code to reproduce the analysis: https://github.com/arnaudmgh/synergy-screen. Raw data: https://raw.githubusercontent.com/arnaudmgh/synergy-screen/master/data/rawscreen.csv
    - Amzallag, Arnaud, Sridhar Ramaswamy, and Cyril H. Benes. “Statistical Assessment and Visualization of Synergies for Large-Scale Sparse Drug Combination Datasets.” BMC Bioinformatics 20, no. 1 (December 2019). https://doi.org/10.1186/s12859-019-2642-7.

- `DGB` - Drug Gene Budger, small molecule prioritization using LINCS L1000, CMap, GEO, CREEDS. Output - drugs that up- or downregulated the selected gene, stratified per database. http://amp.pharm.mssm.edu/DGB/
    - Wang, Zichen, Edward He, Kevin Sani, Kathleen M Jagodnik, Moshe C Silverstein, and Avi Ma’ayan. “Drug Gene Budger (DGB): An Application for Ranking Drugs to Modulate a Specific Gene Based on Transcriptomic Signatures.” Edited by Jonathan Wren. Bioinformatics 35, no. 7 (April 1, 2019): 1247–48. https://doi.org/10.1093/bioinformatics/bty763.

- `CellMinerCDB` - genomics (gene expression, mutations, copy number, methylation, and protein expression) and pharmacogenomics (drug responses and genomics interplay) analyses of cancer cell lines. Integrates NCI-60, GDSC, CCLE, CTRP, and NCI-SCLC databases built on top of `rcellminer` R package. Correlation and multivariate analyses. Tissue-specific analysis. https://discover.nci.nih.gov/cellminercdb/, 10m video tutorial https://youtu.be/XljXazRGkQ8.
    - Rajapakse, Vinodh N., Augustin Luna, Mihoko Yamade, Lisa Loman, Sudhir Varma, Margot Sunshine, Francesco Iorio, et al. “CellMinerCDB for Integrative Cross-Database Genomics and Pharmacogenomics Analyses of Cancer Cell Lines.” IScience 10 (December 2018): 247–64. https://doi.org/10.1016/j.isci.2018.11.029.

- `DGIdb` - drug-gene interaction database integrating 30 sources. API access http://dgidb.org/api. Downloads in text format http://dgidb.org/downloads. http://www.dgidb.org/
    - Cotto, Kelsy C, Alex H Wagner, Yang-Yang Feng, Susanna Kiwala, Adam C Coffman, Gregory Spies, Alex Wollam, Nicholas C Spies, Obi L Griffith, and Malachi Griffith. “DGIdb 3.0: A Redesign and Expansion of the Drug–Gene Interaction Database.” Nucleic Acids Research 46, no. D1 (January 4, 2018): D1068–73. https://doi.org/10.1093/nar/gkx1143.

- `CARE` - biomarker identification from interactions of drug target genes with other genes. Multivariate linear modeling with interaction term. Illustrative example of interaction of BRAF mutation and EGFR expression. Sample separation by gene expression correlation with CARE score better predicts survival. Comparison with correlation, elastic net, support vector regression. http://care.dfci.harvard.edu/, download page http://care.dfci.harvard.edu/download/, nls_logsig tool to compute AUC for dose curves.
    - Jiang, Peng, Winston Lee, Xujuan Li, Carl Johnson, Jun S. Liu, Myles Brown, Jon Christopher Aster, and X. Shirley Liu. “Genome-Scale Signatures of Gene Interaction from Compound Screens Predict Clinical Efficacy of Targeted Cancer Therapies.” Cell Systems 6, no. 3 (March 2018): 343-354.e5. https://doi.org/10.1016/j.cels.2018.01.009.

- `OncoKB` - OncoKB cancer gene database, different levels of evidence, fully downloadable. http://oncokb.org
    - Chakravarty, Debyani, Jianjiong Gao, Sarah M. Phillips, Ritika Kundra, Hongxin Zhang, Jiaojiao Wang, Julia E. Rudolph, et al. “OncoKB: A Precision Oncology Knowledge Base.” JCO Precision Oncology 2017 (July 2017).

- `DepMap` - Large-scale RNAi screen for cancer vulnerability genes in 501 cell lines from 20 cancers, shRNA silencing \~17,000 genes. DEMETER - Modeling and removal of shRNA off-target effects. 6 sigma cutoff of DEMETER scores to identify 769 differential gene dependencies. ATLANTIS model to predict other genes - MDPs, marker dependency pairs. https://www.broadinstitute.org/news/mapping-cancers-vulnerabilities, http://news.harvard.edu/gazette/story/2017/07/first-draft-of-genome-wide-cancer-dependency-map-released/, https://depmap.org/rnai/index, https://depmap.org/portal/depmap/
    - Tsherniak, Aviad, Francisca Vazquez, Phil G. Montgomery, Barbara A. Weir, Gregory Kryukov, Glenn S. Cowley, Stanley Gill, et al. “Defining a Cancer Dependency Map.” Cell 170, no. 3 (July 2017): 564-576.e16. https://doi.org/10.1016/j.cell.2017.06.010.

- `DSigDB` - drug-gene signature database. D1 (approved drugs), D2 (kinase inhibitors), D3 (perturbagent signatures), D4 (computational predictions). Download, online. http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/download.html
    - Yoo, Minjae, Jimin Shin, Jihye Kim, Karen A. Ryall, Kyubum Lee, Sunwon Lee, Minji Jeon, Jaewoo Kang, and Aik Choon Tan. “DSigDB: Drug Signatures Database for Gene Set Analysis: Fig. 1.” Bioinformatics 31, no. 18 (September 15, 2015): 3069–71. https://doi.org/10.1093/bioinformatics/btv313.

- `Drug Repurposing Hub` - drugs with targets, manually curated, experimentally validated. Data, drugs and targets, are at https://clue.io/repurposing#download-data
    - Corsello, Steven M, Joshua A Bittker, Zihan Liu, Joshua Gould, Patrick McCarren, Jodi E Hirschman, Stephen E Johnston, et al. “The Drug Repurposing Hub: A next-Generation Drug Library and Information Resource.” Nature Medicine 23, no. 4 (April 2017): 405–8. https://doi.org/10.1038/nm.4306.

- `CancerRxGene` - Drug-gene targets. Lots of drug sensitivity information. http://www.cancerrxgene.org/
    - Yang, Wanjuan, Jorge Soares, Patricia Greninger, Elena J. Edelman, Howard Lightfoot, Simon Forbes, Nidhi Bindal, et al. “Genomics of Drug Sensitivity in Cancer (GDSC): A Resource for Therapeutic Biomarker Discovery in Cancer Cells.” Nucleic Acids Research 41, no. Database issue (January 2013): D955-961. https://doi.org/10.1093/nar/gks1111.

- `GDA` - Genomics and Drugs integrated Analysis. The Genomics and Drugs integrated Analysis portal (GDA) is a web-based tool that combines NCI60 uniquely large number of drug sensitivity data with CCLE and NCI60 gene mutation and expression profiles. Gene-to-drug and reverse analysis. http://gda.unimore.it/

- `CTRP` - The Cancer Therapeutics Response Portal (CTRP) links genetic, lineage, and other cellular features of cancer cell lines to small-molecule sensitivity with the goal of accelerating discovery of patient-matched cancer therapeutics. https://portals.broadinstitute.org/ctrp/



## Data

- [nstd186 (NCBI Curated Common Structural Variants)](https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd186/) - hg19 genomic coordinates of variant regions

- LINCS, Breast Cancer Profiling Project, Gene Expression 1: Baseline mRNA sequencing on 35 breast cell lines, downloadable matrix of RPKM values. http://lincs.hms.harvard.edu/db/datasets/20348/main

- `MiOncoCirc` - cancer-oriented circRNA database. Exome-capture RNA-seq protocol achieves better enrichment for circRNAs than Ribo-Zero and Rnase R protocols. CIRCExplorer pipeline. Data: https://mioncocirc.github.io/download/
    - Vo, Josh N., Marcin Cieslik, Yajia Zhang, Sudhanshu Shukla, Lanbo Xiao, Yuping Zhang, Yi-Mi Wu, et al. “The Landscape of Circular RNA in Cancer.” Cell 176, no. 4 (February 2019): 869-881.e13. https://doi.org/10.1016/j.cell.2018.12.021.

- `UCSCXenaTools` - An R package downloading and exploring data from UCSC Xena data hubs. CRAN: https://cran.r-project.org/web/packages/UCSCXenaTools/, GitHub: https://github.com/ShixiangWang/UCSCXenaTools, https://shixiangwang.github.io/UCSCXenaTools/. The accompanying package is a Shiny app, https://github.com/openbiox/XenaShiny.
    - Wang, Shixiang, and Xuesong Liu. “The UCSCXenaTools R Package: A Toolkit for Accessing Genomics Data from UCSC Xena Platform, from Cancer Multi-Omics to Single-Cell RNA-Seq.” Journal of Open Source Software 4, no. 40 (August 5, 2019): 1627. https://doi.org/10.21105/joss.01627.

- BREAST CANCER LANDSCAPE RESOURCE, A web portal to proteomics, transcriptomics, genomics and metabolomics of breast cancer. Download http://www.breastcancerlandscape.org/
    - Consortia Oslo Breast Cancer Research Consortium (OSBREAC), Henrik J. Johansson, Fabio Socciarelli, Nathaniel M. Vacanti, Mads H. Haugen, Yafeng Zhu, Ioannis Siavelis, et al. “Breast Cancer Quantitative Proteome and Proteogenomic Landscape.” Nature Communications 10, no. 1 (December 2019): 1600. https://doi.org/10.1038/s41467-019-09018-y. - Proteogenomics of breast cancer subtypes. \~10K proteins by LS-MS/MS. 9 samples for each of the five PAM50 subtypes. Protein expression partially recapitulates PAM50 subtypes, their own consensus clustering. High correlation with mRNA, less so for CNV. Correlation of 290 proteins that are FDA-approved drug targets. Online tool, http://www.breastcancerlandscape.org/, supplementary Data 1 has the full protein expression matrix,https://www.nature.com/articles/s41467-019-09018-y#Sec15

- `Refine.bio` harmonizes petabytes of publicly available biological data into ready-to-use datasets for cancer researchers and AI/ML scientists. https://www.refine.bio/. Documentation, http://docs.refine.bio/en/latest/, GitHub, https://github.com/AlexsLemonade/refinebio.

- Zehir, Ahmet, Ryma Benayed, Ronak H Shah, Aijazuddin Syed, Sumit Middha, Hyunjae R Kim, Preethi Srinivasan, et al. “Mutational Landscape of Metastatic Cancer Revealed from Prospective Clinical Sequencing of 10,000 Patients.” Nature Medicine 23, no. 6 (May 8, 2017): 703–13. https://doi.org/10.1038/nm.4333. - MSK-IMPACT study. Deep sequencing of 341-410 genes in 10,000 samples in multiple cancers. Focus on mutations, copy number alterations, fusions. Data at http://www.cbioportal.org/study?id=msk_impact_2017#summary, downloadable, includes clinical data for survival analysis.

- Gendoo, Deena M.A., Michael Zon, Vandana Sandhu, Venkata Manem, Natchar Ratanasirigulchai, Gregory M. Chen, Levi Waldron, and Benjamin Haibe-Kains. “MetaGxData: Clinically Annotated Breast, Ovarian and Pancreatic Cancer Datasets and Their Use in Generating a Multi-Cancer Gene Signature,” November 12, 2018. https://doi.org/10.1101/052910. - MetaGxData package containing breast and ovarian cancer data, microarray- and RNA-seq gene expression and clinical annotations. Scripts to conduct genome-wide survival analysis for all genes. https://github.com/bhklab/MetaGxData 

- `DepMap` - Large-scale RNAi screen for cancer vulnerability genes in 501 cell lines from 20 cancers, shRNA silencing \~17,000 genes. DEMETER - Modeling and removal of shRNA off-target effects. 6 sigma cutoff of DEMETER scores to identify 769 differential gene dependencies. ATLANTIS model to predict other genes - MDPs, marker dependency pairs. Main data portal: https://depmap.org/portal/download/
    - Tsherniak, Aviad, Francisca Vazquez, Phil G. Montgomery, Barbara A. Weir, Gregory Kryukov, Glenn S. Cowley, Stanley Gill, et al. “Defining a Cancer Dependency Map.” Cell 170, no. 3 (July 2017): 564-576.e16. https://doi.org/10.1016/j.cell.2017.06.010. [Supplemental tables](https://www.sciencedirect.com/science/article/pii/S0092867417306517?via%3Dihub#app2), `DepMap_TableS3_DependencyCorrelation.csv` - Table S3. Gene Dependency-Dependency Correlations, pairs of genes essential for proliferation/viability. Columns: Gene symbol 1, Gene symbol 2, correlation (r), z_score. [Source](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306517-mmc3.csv)

- Pooled shRNA screen of 77 breast cancer cell lines. siMEM algorithm to improve identification of susceptibility/driver genes. Known and novel genes. [Processed gene expression and proteomics data](http://neellab.github.io/bfg/), [Supplementary data](https://www.cell.com/cell/fulltext/S0092-8674(15)01624-4#supplementaryMaterial) with BRCA genes
    - Marcotte, Richard, Azin Sayad, Kevin R. Brown, Felix Sanchez-Garcia, Jüri Reimand, Maliha Haider, Carl Virtanen, et al. “[Functional Genomic Landscape of Human Breast Cancer Drivers, Vulnerabilities, and Resistance](https://doi.org/10.1016/j.cell.2015.11.062).” Cell 164, no. 1–2 (January 14, 2016)
    - [Table S2](data/Marcotte_2016_BRCA/Marcotte_2016_BRCA_mmc3.xls) - BRCA general essential and HER2+ specific genes.
    - [Table S3](data/Marcotte_2016_BRCA/Marcotte_2016_BRCA_mmc4.xls) - subtype-specific essential genes
    - [Table S4](data/Marcotte_2016_BRCA/Marcotte_2016_BRCA_mmc5.xls) - Genes Synthetic Lethal with Specific CNAs, MYCN-interacting genes
    - [Table S5](data/Marcotte_2016_BRCA/Marcotte_2016_BRCA_mmc6.xls) - Genes Associated with Functional Clusters and Drug Response
    - [Table S6](data/Marcotte_2016_BRCA/Marcotte_2016_BRCA_mmc7.xls) - Genes Associated with Expression or Copy Number Loss

- `CCLE2 data` - CCLE characterization using sequencing technologies. Data described: RNA splicing, DNA methylation, Histone modification, miRNA expression, RPPA for 1072 cells. Data availability: https://portals.broadinstitute.org/ccle/data, https://depmap.org/portal/download/
    - Ghandi, Mahmoud, Franklin W. Huang, Judit Jané-Valbuena, Gregory V. Kryukov, Christopher C. Lo, E. Robert McDonald, Jordi Barretina, et al. “Next-Generation Characterization of the Cancer Cell Line Encyclopedia.” Nature, May 8, 2019. https://doi.org/10.1038/s41586-019-1186-3.


### TCGA PanCancer

- Pan-cancer analysis of somatic noncoding driver mutations. Raw data: https://docs.icgc.org/pcawg/data/, Processed data: https://dcc.icgc.org/releases/PCAWG/drivers, significantly recurring breakpoints and juxtapositions http://www.svscape.org/. Extended data figures and tables should be considered individually https://www.nature.com/articles/s41586-020-1965-x#additional-information
    - PCAWG Drivers and Functional Interpretation Working Group, PCAWG Structural Variation Working Group, PCAWG Consortium, Esther Rheinbay, Morten Muhlig Nielsen, Federico Abascal, Jeremiah A. Wala, et al. “Analyses of Non-Coding Somatic Drivers in 2,658 Cancer Whole Genomes.” Nature 578, no. 7793 (February 2020): 102–11. https://doi.org/10.1038/s41586-020-1965-x.

- Pan-Cancer atlas of alternative splicing events, called using SplAdder. Data in GFF3, HDF5, TXT formats: https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Splicing-2018
    - Kahles, André, Kjong-Van Lehmann, Nora C. Toussaint, Matthias Hüser, Stefan G. Stark, Timo Sachsenberg, Oliver Stegle, et al. “Comprehensive Analysis of Alternative Splicing Across Tumors from 8,705 Patients.” Cancer Cell 34, no. 2 (August 2018): 211-224.e6. https://doi.org/10.1016/j.ccell.2018.07.001.

- ATAC-seq data in 410 tumor samples from TCGA (23 cancer types). Correlation with gene expression predicts distal interactions. 18 clusters by cancer type. Data: hg19 coordinates of pan-cancer and BRCA-specific ATAC-seq peaks (Data S2), eQTLs (Data S5), peak-to-gene and enhancer-to-gene links (Data S7), and more https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
    - Corces, M. Ryan, Jeffrey M. Granja, Shadi Shams, Bryan H. Louie, Jose A. Seoane, Wanding Zhou, Tiago C. Silva, et al. “The Chromatin Accessibility Landscape of Primary Human Cancers.” Edited by Rehan Akbani, Christopher C. Benz, Evan A. Boyle, Bradley M. Broom, Andrew D. Cherniack, Brian Craft, John A. Demchok, et al. Science 362, no. 6413 (2018). https://doi.org/10.1126/science.aav1898.

- Papers and supplementary data from PanCancer publications. Clinical annotations, RNA-seq counts, RPPA, Methylation, miRNA, copy number, mutations in .maf format. https://gdc.cancer.gov/about-data/publications/pancanatlas
    - Ding, Li, Matthew H. Bailey, Eduard Porta-Pardo, Vesteinn Thorsson, Antonio Colaprico, Denis Bertrand, David L. Gibbs, et al. “Perspective on Oncogenic Processes at the End of the Beginning of Cancer Genomics.” Cell 173, no. 2 (April 5, 2018): 305-320.e10. https://doi.org/10.1016/j.cell.2018.03.033. - An overview of PanCancer Atlas.

- The Pan-Cancer analysis by TCGA consortium, all papers. https://www.cell.com/pb-assets/consortium/pancanceratlas/pancani3/index.html

- TCGA MC3 variant calling project. Eight variant callers. Protocols for filtering samples, variants. Public and controlled access MAF files at https://gdc.cancer.gov/about-data/publications/mc3-2017
    - Ellrott, Kyle, Matthew H. Bailey, Gordon Saksena, Kyle R. Covington, Cyriac Kandoth, Chip Stewart, Julian Hess, et al. “Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines.” Cell Systems 6, no. 3 (March 2018): 271-281.e7. https://doi.org/10.1016/j.cels.2018.03.002.

- `PCAGW` - The PCAWG study is an international collaboration to identify common patterns of mutation in more than 2,800 cancer whole genomes from the International Cancer Genome Consortium. The project produced large amount data with many types including simple somatic mutations (SNVs, MNVs and small INDELs), large-scale somatic structural variations, copy number alterations, germline variations, RNA expression profiles, gene fusions, and phenotypic annotations etc. PCAWG data have been imported, processed and made available in the following four major online resources for download and exploration by the cancer researchers worldwide. http://docs.icgc.org/pcawg/
    - Goldman, Mary, Junjun Zhang, Nuno A. Fonseca, Qian Xiang, Brian Craft, Elena Piñeiro, Brian O’Connor, et al. “Online Resources for PCAWG Data Exploration, Visualization, and Discovery.” BioRxiv, October 18, 2017. https://doi.org/10.1101/163907. https://www.biorxiv.org/content/early/2017/10/18/163907

### Pediatric

- `PedcBioPortal` - childhood cancer genomics portal. 261 pediatric PDX models, 37 pediatric malignancies, CNS and rhabdoid tumors, extracranial solid tumors, hematological malignancies. WES (includes segmentation, copy number estimates, high breakpoint density), RNA-seq. Description of individual tumor types. hg19-mm10 alignment. Raw data: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001437.v1.p1, Processed data: https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147, Code at https://github.com/marislab, Online cBioPortal interface: https://pedcbioportal.kidsfirstdrc.org/
    - Rokita, Jo Lynne, Komal S. Rathi, Maria F. Cardenas, Kristen A. Upton, Joy Jayaseelan, Katherine L. Cross, Jacob Pfeil, et al. “Genomic Profiling of Childhood Tumor Patient-Derived Xenograft Models to Enable Rational Clinical Trial Design.” Cell Reports 29, no. 6 (November 2019): 1675-1689.e9. https://doi.org/10.1016/j.celrep.2019.09.071.

- TARGET pediatric tumors RNA-sequencing dataset: https://ocg.cancer.gov/programs/target/data-matrix

## Methylation

- `MEXPRESS` - Gene-centric methylation and correlation with clinical parameters. http://mexpress.be/

- `Pancan-meQTL` database of meQTLs across 23 TCGA cancer types. Cis-, trans-meQTLs, pancancer-meQTLs, survival meQTLs. SNP-, gene-, CpG-centric search for each cancer. Visualization, KM plots for survival. Download. http://bioinfo.life.hust.edu.cn/Pancan-meQTL/
    - Gong, Jing, Hao Wan, Shufang Mei, Hang Ruan, Zhao Zhang, Chunjie Liu, An-Yuan Guo, Lixia Diao, Xiaoping Miao, and Leng Han. “Pancan-MeQTL: A Database to Systematically Evaluate the Effects of Genetic Variants on Methylation in Human Cancer.” Nucleic Acids Research, September 7, 2018. https://doi.org/10.1093/nar/gky814.

- [Wanderer](http://gattaca.imppc.org:3838/wanderer/index.html) - An interactive viewer to explore DNA methylation and gene expression data in human cancer


## Misc

- [CancerSubtypes](http://bioconductor.org/packages/release/bioc/html/CancerSubtypes.html) - an R package implementing four methods for clustering/subtype identification. Includes consensus nonnegative matrix factorization (CNMF), iCluster, SNF and its variations. Better distinguish survival groups, single file format
    - Xu, Taosheng, Thuc Duy Le, Lin Liu, Ning Su, Rujing Wang, Bingyu Sun, Antonio Colaprico, Gianluca Bontempi, and Jiuyong Li. “[CancerSubtypes: An R/Bioconductor Package for Molecular Cancer Subtype Identification, Validation and Visualization](https://doi.org/10.1093/bioinformatics/btx378).” Bioinformatics 33, no. 19 (October 1, 2017)

- Oncology Model Fidelity Score based on the Hallmarks of Cancer. Quantify fidelity of PDX models in recapitulating human cancers. https://github.com/tedgoldstein/hallmarks

-  SMAP is a pipeline for the process of xenografts sequencing data. It takes FASTQ as input and outputs specie-specific BAM or gene counts. https://github.com/cit-bioinfo/SMAP

