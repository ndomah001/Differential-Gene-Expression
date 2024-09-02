Differential Gene Expression Analysis
================

# Setup

## Load libraries

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
    ##     tapply, union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(pheatmap)
library(RColorBrewer)
```

## Read count data

``` r
count_table <- read.csv('raw_counts.tsv', sep='\t', row.names=1)
head(count_table)
```

    ##             GSM6160812 GSM6160813 GSM6160814 GSM6160821 GSM6160822 GSM6160823
    ## DDX11L1              8         11          9         10          6          6
    ## WASH7P             329        377        345        350        324        317
    ## MIR6859-1            5          7          2          6          6          7
    ## MIR1302-2HG          1          1          0          0          0          0
    ## MIR1302-2            0          1          0          0          0          0
    ## FAM138A              0          0          0          0          0          0

``` r
dim(count_table)
```

    ## [1] 39374     6

## Read the sample information

``` r
sample_info <- read.csv('design.tsv', sep='\t', row.names=1)
sample_info
```

    ##               Group
    ## GSM6160812  control
    ## GSM6160813  control
    ## GSM6160814  control
    ## GSM6160821 tgf-beta
    ## GSM6160822 tgf-beta
    ## GSM6160823 tgf-beta

``` r
dim(sample_info)
```

    ## [1] 6 1

## Set factor levels

``` r
factors <- factor(sample_info$Group)
groups <- unique(sample_info$Group)

sample_info$Group
```

    ## [1] "control"  "control"  "control"  "tgf-beta" "tgf-beta" "tgf-beta"

``` r
groups
```

    ## [1] "control"  "tgf-beta"

## We want the control group to be after tgf-beta

``` r
groups <- rev(groups)
groups
```

    ## [1] "tgf-beta" "control"

``` r
sample_info$Group <- factors
sample_info$Group
```

    ## [1] control  control  control  tgf-beta tgf-beta tgf-beta
    ## Levels: control tgf-beta

## Create DESeq object

``` r
dds <- DESeqDataSetFromMatrix(countData=count_table, colData=sample_info, design=~Group)
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

## Set the reference for the Group factor

``` r
dds$Group <- relevel(dds$Group, ref='control')
```

## Filter out low gene counts

Keep genes with at least N counts \>= 10, where N = size of smallest
group

``` r
keep <- rowSums(counts(dds)>=10) >= min(table(sample_info$Group))
dds <- dds[keep,]
```

## Perform statistical tests and get result

``` r
dds <- DESeq(dds, test='Wald', sfType='poscount')
```

    ## estimating size factors

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## final dispersion estimates

    ## fitting model and testing

``` r
deseq_result <- results(dds)
deseq_result <- as.data.frame(deseq_result)
head(deseq_result)
```

    ##               baseMean log2FoldChange     lfcSE       stat     pvalue
    ## WASH7P       340.27254    -0.04113513 0.1122334 -0.3665141 0.71398146
    ## LOC100996442  43.85160    -0.57873959 0.2881260 -2.0086339 0.04457597
    ## LOC729737     46.67727    -0.57202894 0.2806266 -2.0383985 0.04151009
    ## DDX11L17      19.31418    -0.30600789 0.4342199 -0.7047302 0.48097813
    ## WASH9P       388.81481    -0.04842712 0.1021837 -0.4739224 0.63555524
    ## LOC100132287  72.80562    -0.38181973 0.2341107 -1.6309366 0.10290371
    ##                    padj
    ## WASH7P       0.80895518
    ## LOC100996442 0.09527542
    ## LOC729737    0.08963263
    ## DDX11L17     0.61729290
    ## WASH9P       0.74867638
    ## LOC100132287 0.19015899

``` r
dim(deseq_result)
```

    ## [1] 18314     6

## Add ‘GeneName’ column, save to tsv file

``` r
deseq_result$GeneName <- row.names(deseq_result)

deseq_result <- subset(deseq_result,
                       select=c('GeneName', 'padj', 'pvalue', 'lfcSE', 'stat',  'log2FoldChange', 'baseMean'))

head(deseq_result)
```

    ##                  GeneName       padj     pvalue     lfcSE       stat
    ## WASH7P             WASH7P 0.80895518 0.71398146 0.1122334 -0.3665141
    ## LOC100996442 LOC100996442 0.09527542 0.04457597 0.2881260 -2.0086339
    ## LOC729737       LOC729737 0.08963263 0.04151009 0.2806266 -2.0383985
    ## DDX11L17         DDX11L17 0.61729290 0.48097813 0.4342199 -0.7047302
    ## WASH9P             WASH9P 0.74867638 0.63555524 0.1021837 -0.4739224
    ## LOC100132287 LOC100132287 0.19015899 0.10290371 0.2341107 -1.6309366
    ##              log2FoldChange  baseMean
    ## WASH7P          -0.04113513 340.27254
    ## LOC100996442    -0.57873959  43.85160
    ## LOC729737       -0.57202894  46.67727
    ## DDX11L17        -0.30600789  19.31418
    ## WASH9P          -0.04842712 388.81481
    ## LOC100132287    -0.38181973  72.80562

``` r
write.table(deseq_result, file='deseq_result.all.tsv', row.names=F, sep='\t')
```

## Extract DE genes with padj \< 0.05 and log2foldchange \<= -1 or \>= 1

``` r
deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange)>=1)
dim(deg)
```

    ## [1] 1491    7

``` r
dim(deseq_result)
```

    ## [1] 18314     7

## Order by padj ascending

``` r
deg <- deg[order(deg$padj),]
head(deg)
```

    ##         GeneName padj pvalue      lfcSE      stat log2FoldChange  baseMean
    ## ELF3        ELF3    0      0 0.04840507 -39.07177      -1.891272  4685.371
    ## LBH          LBH    0      0 0.04444850  41.78416       1.857243  7779.323
    ## MARCHF4  MARCHF4    0      0 0.06907873  46.65602       3.222939  2192.817
    ## TNS1        TNS1    0      0 0.06923313  46.68812       3.232365  2183.254
    ## F2RL1      F2RL1    0      0 0.04238652  46.13813       1.955635 11797.268
    ## TGFBI      TGFBI    0      0 0.03618620  44.19816       1.599363 51441.326

## Write to tsv file

``` r
write.table(deg, file='deseq_deg.tsv', row.names=F, sep='\t')
```

# Gene Expression Data Visualization

## Plot dispersion estimates

``` r
plotDispEsts(dds, main='GSE203159 Dispersion Estimates')
```

![](Differential-Gene-Expression-Analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

## Create histogram of p-values

``` r
hist(deseq_result$padj, breaks=seq(0,1,length=21), col='grey', border='white', xlab='', ylab='', ylim=c(0,8000), main='GSE203159 Frequencies of padj-values')
```

![](Differential-Gene-Expression-Analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

## Volcano plot

### Set colors

``` r
old.pal <- palette(c('#00BFFF', '#FF3030'))
```

### Set margin size

``` r
par(mar=c(4,4,2,1), cex.main=1.5)
```

### Plot values

Add legend for up- and down-regulation

``` r
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), main='tgf-beta vs control', xlab='log2FC', ylab='-log10(Padj)', pch=20, cex=0.5)

with(subset(deseq_result, padj<0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange)+3)/2, cex=1))

legend('bottomleft', title=paste('Padj<', 0.05, sep=''), legend=c('down', 'up'), pch=20, col=1:2)
```

![](Differential-Gene-Expression-Analysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

## PCA plot

### Variance stabilizing transformation

``` r
vsd <- vst(dds, blind=FALSE)
```

### Use transformed values to generate PCA plot

``` r
plotPCA(vsd, intgroup=c('Group'))
```

    ## using ntop=500 top features by variance

![](Differential-Gene-Expression-Analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

*Note: I couldn’t find a fix to make the PCA plot’s y-axis longer. In
the future, I could manually plot the points but for the sake of this
project I will continue.*

## Heatmaps

Heatmap of log transformed normalized counts using the top 10 genes

``` r
normalized_counts <- counts(dds, normalized=T)
head(normalized_counts)
```

    ##              GSM6160812 GSM6160813 GSM6160814 GSM6160821 GSM6160822 GSM6160823
    ## WASH7P        323.27898  370.58129  341.44322  359.25993  324.58495  322.48689
    ## LOC100996442   53.06099   52.09763   52.45360   33.87308   39.07041   32.55388
    ## LOC729737      53.06099   60.94440   53.44329   36.95245   41.07402   34.58850
    ## DDX11L17       17.68700   21.62543   24.74226   18.47622   14.02528   19.32887
    ## WASH9P        390.09652  397.12159  398.84527  383.89490  397.71674  365.21386
    ## LOC100132287   89.41759   86.50173   71.25772   63.64033   74.13360   51.88275

``` r
transformed_counts <- log2(normalized_counts + 1)
head(transformed_counts)
```

    ##              GSM6160812 GSM6160813 GSM6160814 GSM6160821 GSM6160822 GSM6160823
    ## WASH7P         8.341092   8.537534   8.419721   8.492894   8.346890   8.337563
    ## LOC100996442   5.756516   5.730576   5.740215   5.124042   5.324465   5.068408
    ## LOC729737      5.756516   5.952902   5.766682   5.246121   5.394858   5.153339
    ## DDX11L17       4.223963   4.499874   4.686067   4.283642   3.909320   4.345458
    ## WASH9P         8.611381   8.637065   8.643298   8.588321   8.639220   8.516543
    ## LOC100132287   6.498532   6.451240   6.175080   6.014363   6.231386   5.724725

``` r
top_hits <- row.names(deg[1:10, ])
top_hits
```

    ##  [1] "ELF3"     "LBH"      "MARCHF4"  "TNS1"     "F2RL1"    "TGFBI"   
    ##  [7] "PCDH1"    "ADAM19"   "SERPINE1" "TGFBR1"

``` r
top_hits <- transformed_counts[top_hits,]
head(top_hits)
```

    ##         GSM6160812 GSM6160813 GSM6160814 GSM6160821 GSM6160822 GSM6160823
    ## ELF3     12.801244  12.860544  12.886295   10.96987   10.94365   10.96359
    ## LBH      11.702620  11.735953  11.712123   13.62275   13.54636   13.55138
    ## MARCHF4   8.765913   8.644172   8.782745   12.03060   11.91587   11.90663
    ## TNS1      8.769165   8.798648   8.573803   11.99948   11.90761   11.93164
    ## F2RL1    12.233558  12.251646  12.234462   14.17586   14.16373   14.24506
    ## TGFBI    14.653403  14.620772  14.645269   16.25472   16.21870   16.24396

``` r
pheatmap(top_hits, cluster_rows=FALSE, cluster_cols=FALSE)
```

![](Differential-Gene-Expression-Analysis_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->
