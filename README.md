# Differential-Gene-Expression-Analysis
## 📝 Background
This project was inspired by [this study from ScienceSignaling](https://www.science.org/doi/10.1126/scisignal.adf1947#sec-4) that identified a long noncoding RNA (lncRNA) named LETS1, which amplifies the TGF-β signaling pathway, known for driving epithelial-to-mesenchymal transition (EMT) and cancer progression. LETS1 stabilizes the TGF-β type I receptor (TβRI) by interacting with NFAT5, leading to the suppression of SMAD7, a protein that typically degrades TβRI. This creates a feedback loop that enhances TGF-β signaling, promoting cancer cell migration and invasion. LETS1's role was confirmed in breast and lung cancer cells, suggesting it may contribute to cancer metastasis in patients.

I then performed a differential gene expression analysis (in R using DESeq2) on the RNA-seq counts extracted from [NCBI GEO (Gene Experssion Omnibus) - GSE203119](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203119). 

## 📊 Visualizations
### Dispersion Estimates
Visualize the variability of gene expression across samples. 
![de](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/Dispersion%20Estimates.png)
### Histogram of p-values
The peak near 0 indicates a high number of differentially expressed genes. 
![h](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/p-value%20Histogram.png)
### Volcano Plot
Visualize the relationship between the magnitude of change (log fold change) and statistical significance (p-value) of genes, as well as up- and down-regulation. 
![vp](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/Volcano%20Plot.png)
### PCA Plot
Reduce dimensionality of data and visualize sample relationships.
![pca](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/PCA%20Plot.png)
### Heatmap
Visualize expression levels of log transformed normalized counts using top 10 genes. 
![hmp](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/Heatmap%20(top%2010%20genes).png)
