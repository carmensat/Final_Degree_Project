# **Deciphering the role of Dan and Danr in the identity of Lamina Neurons in the Drosophila Visual System**


In this repository, you will find all the code, images, and files related to my final degree project, which consisted in generating a Gene Regulatory Network for the Lamina cells in the *Drosophila Melanoagaster*


# Abstract


**Motivation:** Dan and Danr play a role in different mechanisms that help orchestrate the neuronal identity in the Drosophila Visual System, and yet their role in the lamina cells remains unclear. 


**Results:** We employed confocal imaging, scRNA-seq, snATAC-seq, and a gene regulatory network (GRN) inference analysis to investigate the regulatory function of these genes in the lamina.
Danr emerged as a top regulator, repressing many targets, including another known Terminal Selector, while Dan showed minimal connectivity, primarily repressing Danr. 
These findings suggest Danr may act not as a Terminal Selector, but as a key modulator of neuronal maturation and subtype restriction.


**Supplementary information:** Supplementary data are available here.
If you have any questions or troubles while executing, don't hesitate to get in touch with me through carmen.samedi@alum.esci.upf.edu.

All single-cell RNA sequencing and single-nucleus ATAC sequencing analysed in this project were obtained from the publicly available Gene Expression Omnibus (GEO) repository, accession numbers [GSE167266](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) and [GSE163697](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163697), respectively.


# Repository Files
In the repository, you find many files used for either the execution or preparation of the GRN. You can also find the analysis run on other files, the results, and the files necessary for execution. 
For the Network inference, I chose to use the Inferelator 3 (v0.5.6), a robust and widely adopted framework for inferring gene regulatory networks from high-dimensional transcriptomic data. It combines prior knowledge of transcription factor binding with gene expression to estimate transcription factor activity (TFA) and predict direct regulatory relationships. This model excels at handling sparse data and supports multiple regression techniques, including Bayesian Best Subset Regression (BBSR), which I used for its accuracy and interpretability. Its modular design also facilitated integration with custom motif-based priors, making it well suited for analyzing single-cell and multiomic datasets like those in this project.


The following explains each file and folder. 


The folder `prior_inferelator` contains the shell file that executes the command to find the Transcription Factor binding motifs with the help of FIMO scanning ±10 kb around each gene in the *Drosophila melanogaster* genome (BDGP6.54.114). 
The Inferelator version was Inferelator-Prior v0.2.3, which, in addition to Python dependencies, this package also requires STAR, sra-tools, bedtools, samtools, and fimo. More information regarding the installation and version required at this [link](https://github.com/flatironinstitute/inferelator-prior). 
All the required files for the execution are :
- The `.bed` file, which was extracted from the snATAC-seq and contains the peaks from the ATAC-seq.  
- The `.meme` file, which contains known transcription factor binding motifs of *Drosophila*.  
- The `.fasta` file containing the actual DNA sequences of the desired species, in this case, *Drosophila melanogaster*.  
- The `.gtf` file containing the gene annotation of the species.  
- And finally, the prior matrix, previously built from the peaks file.


The next folder is the `network_inferelator`; this folder contains the inferelator code in the file named `lamina_inferelator.py`. The code specifies all the fine-tunings the Inferelator needs for proper execution. 
The other file in this folder is the `run_inferelator.sh`, which does as its name states and helps run the code using shell. 


The `dan_danr.Rmd` is an RMarkdown file, containing the preliminary analysis run on an available CisTopic representing the ATAC-seq. This file lacked metadata, restricting its utility; we therefore did not pursue further analysis beyond this one. 


The file `viusalisation.py` is the Python script used to visualize the inferred GRN as an interactive network. This code was developed based on scripts generously shared by Giuseppe Saldi, whose guidance was instrumental in its implementation.

The file `scrna_seq.r` is an R script that performs single-cell RNA-seq of Drosophila lamina neurons. It includes cell clustering, marker identification, and gene ontology enrichment.


