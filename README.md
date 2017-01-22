# 2017_ms_immunostim
Analyze the Atlantic Salmon and sea lice immunostimulant project 

### Requirements:
limma `https://bioconductor.org/packages/release/bioc/html/limma.html`    
ReadqPCR `https://bioconductor.org/packages/release/bioc/html/ReadqPCR.html`    
NormqPCR `https://www.bioconductor.org/packages/release/bioc/html/NormqPCR.html`    
stringi    
tidyr    

### Setup:
Put raw data files (block and channel separated) into 02_raw_data  
See file `00_archive/raw_data_requirements.txt` for information on where to obtain raw data files and where to put them    

## Overview:
1. Louse transcriptome analysis
2. Louse qPCR analysis
3. Salmon physiology and lice counts
4. Salmon qPCR analysis

## 1. Louse transcriptome analysis 
With this section, you can take the microarray files and normalize, visualize, quality filter, extract the expressed genes list, save probes to correlate with qPCR data (step 2), and perform differential expression analysis.

## 2. Louse qPCR analysis
With this section, you take non-normalized qPCR data, calculate the best normalizers with geNORM, normalize data, and correlate qPCR data with probes from step 1.

## 3. Salmon physiology and lice counts    
With this section, you can evaluate fish weights before the experiment, during the experiment, generate figures and statistical models of lice per fish for each group over time, and determine the proportions of each lice stage and sex for each group.     

## 4. Salmon qPCR analysis
With this section, you can transform qPCR data, generate the skin qPCR plot in the paper, and generate all of the plots and models used to analyze the salmon gene expression response for genes il-1, il-8, mmp9 and tlr9.    
