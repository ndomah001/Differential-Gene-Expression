# Differential Gene Expression Analysis
## ❓ What is Differential Gene Expression?
Differential gene expression (DGE) is a method used to identify changes in gene expression levels between different conditions or groups. It helps researchers understand how genes are regulated under various circumstances, such as between healthy and diseased tissues, or between different developmental stages.

In a typical DGE analysis, you compare gene expression profiles from two or more conditions (e.g., treated vs. untreated samples). The goal is to determine which genes are significantly upregulated or downregulated in one condition relative to another. This can reveal insights into the biological processes underlying the conditions being studied and identify potential biomarkers or therapeutic targets. The most common way to perform an analysis is using the DESeq2 package in R.

## 📝 Background
This project was inspired by [this study from ScienceSignaling](https://www.science.org/doi/10.1126/scisignal.adf1947#sec-4) that identified a long noncoding RNA (lncRNA) named LETS1, which amplifies the TGF-β signaling pathway, known for driving epithelial-to-mesenchymal transition (EMT) and cancer progression. LETS1 stabilizes the TGF-β type I receptor (TβRI) by interacting with NFAT5, leading to the suppression of SMAD7, a protein that typically degrades TβRI. This creates a feedback loop that enhances TGF-β signaling, promoting cancer cell migration and invasion. LETS1's role was confirmed in breast and lung cancer cells, suggesting it may contribute to cancer metastasis in patients.

I then performed a differential gene expression analysis (in R using DESeq2) on the RNA-seq counts extracted from [NCBI GEO (Gene Experssion Omnibus) - GSE203119](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203119). 

## 📊 Results
### Dispersion Estimates
Dispersion reflects the variability of gene expression, and accurately estimating dispersion is crucial for reliable differential expression analysis.

- The general trend shows that as the mean expression increases, the dispersion decreases. This pattern is expected in RNA-seq data, where lowly expressed genes typically have higher variability.
- The fitted red line represents the dispersion model's fit, and the blue dots indicate shrinkage applied to the gene-wise dispersion estimates to stabilize them, especially for genes with low counts.
- The shrinkage helps control type I error rates in downstream differential expression testing by borrowing information across genes, which leads to more reliable results.


![de](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/Dispersion%20Estimates.png)

### p-value Distribution
- The majority of the adjusted p-values are clustered near 0. This suggests that a large number of genes show statistically significant differential expression (with padj-values close to 0) when applying multiple testing corrections.
- The rest of the distribution shows a much smaller number of genes with adjusted p-values spread across the 0.1 to 1 range. This suggests that a smaller number of genes are not significantly differentially expressed.
- The concentration near 0 suggests that there may be a significant number of genes with strong differential expression between the conditions.
- Since the study focuses on the role of LETS1 and its involvement in TGF-β signaling and EMT, genes that are differentially expressed with low padj-values may include key players in these pathways.


![h](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/p-value%20Histogram.png)

### Volcano Plot
The volcano plot shows the results of the differential expression analysis between TGF-β treated and control samples. Each point represents a gene, with the x-axis representing the log2 fold change (log2FC) and the y-axis representing the negative log10 of the adjusted p-value (Padj).
- **Significance Threshold**: The plot highlights significantly differentially expressed genes with Padj < 0.05. Blue points represent downregulated genes, and red points represent upregulated genes.
- **Upregulated Genes**: Genes on the right side (positive log2FC) are upregulated in TGF-β treated samples compared to controls. This includes genes potentially involved in TGF-β signaling or downstream effects like EMT.
- **Downregulated Genes**: Genes on the left side (negative log2FC) are downregulated in response to TGF-β treatment, suggesting a possible suppression of pathways that inhibit EMT or cancer progression.
- **Central Black Dots**: These represent genes with no significant change in expression, suggesting they are not impacted by TGF-β treatment within the set significance threshold.
- **Biological Implications**: The plot underscores the impact of TGF-β signaling in modulating gene expression, corroborating its role in driving processes like EMT and cancer progression. This aligns with the hypothesis that LETS1 and other genes influenced by TGF-β can have profound effects on cellular behavior, potentially contributing to metastasis.


![vp](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/Volcano%20Plot.png)

### PCA Plot
The Principal Component Analysis (PCA) plot visually summarizes the variance between the control and TGF-β-treated samples.
- **PC1 and PC2**: The plot shows two principal components (PC1 and PC2), with PC1 accounting for 99% of the variance. This indicates that the primary difference between the samples can be captured along the PC1 axis.
- The clear separation of orange (control) and blue (TGF-β-treated) samples suggests that the expression profiles are distinctly different between these two groups, further supporting the impact of TGF-β on gene expression.
- **Clustering**: The tight clustering of control samples (in orange) and treated samples (in blue) confirms that the experimental conditions (TGF-β treatment vs. control) drive the majority of the variance in the dataset.


![pca](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/PCA%20Plot.png)

### Heatmap
The heatmap presents the expression levels of the top 10 differentially expressed genes across different samples (GSM identifiers).
- **Genes**: The rows represent specific genes that are most differentially expressed in the dataset. These include ELF3, LBH, MARCHF4, TNS1, F2RL1, TGFB1, PCDH1, ADAM19, SERPINE1, and TGFBR1.
  - Notably, TGFBR1 and TGFB1, which are part of the TGF-β signaling pathway, show elevated expression in some samples.
  - SERPINE1 is known to be involved in cancer progression and EMT, which aligns with the study's focus on TGF-β-driven EMT.
- **Samples**: The columns represent different RNA-seq samples, categorized by their respective groups (likely control and TGF-β treated, based on the PCA plot).
  - There seems to be a clear distinction in expression patterns between samples, particularly for genes like TGFB1 and SERPINE1, which show higher expression levels in specific samples (indicated by red color).
- **Color Scale**: The color gradient indicates the level of gene expression, with red showing higher expression and blue showing lower expression. For instance, samples GSM6106021 and GSM6106022 exhibit higher expression of TGFB1, SERPINE1, and TGFBR1 compared to others, suggesting these genes are upregulated in these samples.


![hmp](https://github.com/ndomah001/Differential-Gene-Expression-Analysis/blob/main/Heatmap%20(top%2010%20genes).png)
