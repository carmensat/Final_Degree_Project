# Deciphering the role of pipsqueak family genes Dan/Danr in the identity of Lamina Neurons in the Drosophila Visual System 



In this repository, you will find all the code, images, and files related to my final degree project, which consisted of analyzing the gene expression of two specific genes, Dan and Danr, in the Lamina cells of *Drosophila Melanogaster*. This is an archival research repository; code is not intended for reuse and may require lab-specific data and environment.


# Abstract


**Motivation:** Our single-cell transcriptomics data show that Dan/Danr genes are widely expressed in the neural stem cells and neuronal populations in the Drosophila central brain and optic lobe. 
Yet, our understanding of the role of these genes in brain development remains unclear.  
We analysed recently published single-cell transcriptomic data and observed expression of Dan/Danr in the lamina neurons of the Drosophila optic lobe.  
However, whether Dan/Danr orchestrate the neuronal identity in lamina neurons of the Drosophila visual system remains unclear.  


**Results:** We used immunofluorescence microscopy to understand the expression of Dan/Danr in lamina neurons. This was followed by single-cell multiome, gene regulatory network (GRN) inference, and Drosophila genetics to investigate the role of Dan/Danr in the neuronal specification of lamina neurons. 
Danr emerged as a top regulator, repressing many targets, including another known Terminal Selector, while Dan showed minimal connectivity, primarily repressing Danr. 
These findings suggest Danr may act not as a Terminal Selector, but as a key modulator of neuronal maturation and subtype restriction.


**Supplementary information:** Supplementary data are available here.
If you have any questions or issues while executing, please don't hesitate to contact me at cjs9968@nyu.edu

All single-cell RNA sequencing and single-nucleus ATAC sequencing analysed in this project were obtained from the publicly available Gene Expression Omnibus (GEO) repository, accession numbers [GSE167266](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167266) and [GSE163697](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163697), respectively.


# Repository Contents
In the repository, you find many files used for either the execution or preparation of the GRN. You can also find the analysis run on other files, the results, and the files necessary for execution. 
For the Network Inference, I chose to use the Inferelator 3 (v0.5.6), a robust and widely adopted framework for inferring gene regulatory networks from high-dimensional transcriptomic data. It combines prior knowledge of transcription factor binding with gene expression to estimate transcription factor activity (TFA) and predict direct regulatory relationships. This model excels at handling sparse data and supports multiple regression techniques, including Bayesian Best Subset Regression (BBSR), which I used due to its accuracy and interpretability. Its modular design also facilitated integration with custom motif-based priors, making it well-suited for analyzing single-cell and multiomic datasets like those in this project.


The following explains each file and folder. 


The folder [`figures`](figures) contains a PDF file with all the figures that are in the main manuscript

The folder [`prior_inferelator`](prior_inferelator) contains the shell file that executes the command to find the Transcription Factor binding motifs with the help of FIMO scanning ±10 kb around each gene in the *Drosophila melanogaster* genome (BDGP6.54.114) 
Almost all the necessary files are found in that folder as well. 
The Inferelator version was Inferelator-Prior v0.2.3, which, in addition to Python dependencies, this package also requires STAR, sra-tools, bedtools, samtools, and fimo. More information regarding the installation and version required at this [link](https://github.com/flatironinstitute/inferelator-prior). 
All the required files for the execution are :
- The [`.bed`](prior_inferelator/GSE163697_consensus_peaks_500pb.bed.gz) file, which was extracted from the snATAC-seq and contains the peaks from the ATAC-seq.  
- The [`.meme`](prior_inferelator/CisBPDrosophilaALL.meme.zip) file, which contains known transcription factor binding motifs of *Drosophila*.  
- The `.fasta` file containing the actual DNA sequences of the desired species, in this case, *Drosophila melanogaster*.  
- The [`.gtf`](prior_inferelator/Drosophila_melanogaster.BDGP6.54.114.zip) file containing the gene annotation of the species.  
- And finally, the prior matrix, previously built from the peaks file.


The next folder is the [`network_inferelator`](network_inferelator); this folder contains the inferelator code in the file named [`lamina_inferelator.py`](network_inferelator/lamina_inferelator.py). The code specifies all the fine-tunings the Inferelator needs for proper execution. 
The other file in this folder is the [`run_inferelator.sh`](network_inferelator/run_inferelator.sh), which does as its name states and helps run the code using shell.
There is also [`FlyBaseTFs.txt`](network_inferelator/FlyBaseTFs.txt), which contains a list of known or candidate Transcription Factors (TFs) from FlyBase, to help the inferelator define the "regulators" better.


The folder [`results`](figures) contains the different HTML files with the interactive network built with the code `visualisation.py`

Additional files may be found in the folder [`archive`](figures). 
