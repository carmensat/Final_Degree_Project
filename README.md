# Final_Degree_Project

In this repository, you will find all the code and supplementary images related to this project. 
If you have any questions or troubles while executing, please contact me through carmen.samedi@alum.esci.upf.edu.

The folder prior_inferelator contains the shell file that executes the command to find the Transcription Factor binding motifs with the help of FIMO scanning ±10 kb around each gene in the Drosophila melanogaster genome (BDGP6.54.114). 
All the required files for the execution are :
- The peaks file, which was extracted from the snATAC-seq.  
- The `.meme` file, which contains known transcription factor binding motifs of *Drosophila*.  
- The `.fasta` file containing the actual DNA sequences of the desired species, in this case, *Drosophila melanogaster*.  
- The `.gtf` file containing the gene annotation of the species.  
- And finally, the prior matrix, previously built from the peaks file.

The next folder is the network_inferelator; this folder contains the inferelator code in the file named lamina_inferelator.py. The code specifies all the fine-tunings the Inferelator needs for proper execution. 
The other file in this folder is the run_inferelator.sh, which does as its name states and helps run the code using shell. 

The dan_danr.Rmd is an RMarkdown file, containing the preliminary analysis run on an available CisTopic representing the ATAC-seq. This file was restricted and contained no metadata; we therefore did not pursue further analysis beyond this one. 

The file viusalisation.py is the code allowing to visualize the resulting matrix from the inferelator execution. This code was developed based on scripts generously shared by Giuseppe Saldi, whose guidance was instrumental in its implementation.


