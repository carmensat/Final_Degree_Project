# Final_Degree_Project

In this repository, you will find all the code, images, and files related to my final degree project, which consisted in generating a Gene Regulatory Network fo the Lamina cells in the *Drosophila Melanoagaster*:

**Deciphering the role of Dan and Danr in the identity of Lamina Neurons in the Drosophila Visual System**


Carmen Samedi¹


Scientific director: Dr. Asif Bakshi¹, Dr Claude Desplan² 


¹Center for Genomics and Systems Biology (CGSB), New York University Abu Dhabi, PO Box 129188, Abu Dhabi, United Arab Emirates


²Department of Biology, New York University, 100 Washington Place, New York, NY 10003, USA

# Abstract


**Motivation:** Dan and Danr play a role in different mechanisms that help orchestrate the neuronal identity in the Drosophila Visual System, and yet their role in the lamina cells remains unclear. 


**Results:** We used confocal imaging, scRNA-seq, snATAC-seq, and a gene regulatory network (GRN) inference to investigate their regulatory function in the lamina cells.
Danr emerged as a top regulator, repressing many targets, including another known Terminal Selector, while Dan showed minimal connectivity, primarily repressing Danr. 
These findings suggest Danr may act not as a Terminal Selector, but as a key modulator of neuronal maturation and subtype restriction.


**Supplementary information:** Supplementary data are available here.
If you have any questions or troubles while executing, please contact me through carmen.samedi@alum.esci.upf.edu.

All single-cell RNA sequencing and single-nucleus ATAC sequencing analysed in this project were obtained from the publicly available Gene Expression Omnibus (GEO) repository, accession numbers [GSE167266](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) and [GSE163697](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163697), respectively.


# Repository Files
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


