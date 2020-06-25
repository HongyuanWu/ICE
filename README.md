
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICE

<!-- badges: start -->

<!-- badges: end -->

ICE (**I**mpute **C**ircRNA **E**xpression) is designed for imputing
circRNA expression using PCGs (protein-coding genes).

## Installation

> Currently bundled example data may be huge in size. One could try use
> proxy and pull it first.

``` r
remotes::install_github("bioinformatist/ICE")
```

## Quick start

### `library()` package and load example data

``` r
library(ICE)
```

> Currently data used for testing is too large to make a bundle. You
> should follow [example dataset](#example-dataset) part to prepare
> them.

## Performance

### ML and DL on 871 circRNAs

``` r
library(ICE)
library(ggplot2)
library(data.table)

cv.full.pearson <- list(
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', method = 'KNN'),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', method = 'RF', ncores = 56),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', method = 'EN', ncores = 56),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', method = 'lasso', ncores = 56),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', method = 'SVM', ncores = 56))
names(cv.full.pearson) <- c('KNN', 'RF', 'EN', 'lasso', 'SVM')
cv.full.pearson <- lapply(names(cv.full.pearson), function(x) data.table(cv.full.pearson[[x]], method = x))
cv.full.pearson <- do.call(rbind, cv.full.pearson)
ggplot() + geom_boxplot(data = cv.full.pearson, aes(x = method, y = CC, color = method))
ggplot() + geom_density(data = cv.full.pearson[method == 'RF'], aes(CC))
ggplot() + geom_density(data = cv.full.pearson[method == 'RF'], aes(CC)) + scale_x_continuous(breaks = seq(0.5, 0.8, 0.025), labels = seq(0.5, 0.8, 0.025), limits = c(0.4, 0.825))
all.CC <- rbind(cv.full.pearson[, .(CC, method)], data.table(DL, method = 'DL'))

require("reticulate")

source_python("pickle_reader.py")
DL <- read_pickle_file("cc_pancancer.pkl")
names(DL) <- 'CC'
all.CC <- rbind(cv.full.pearson[, .(CC, method)], data.table(DL, method = 'DL'))
ggplot() + geom_density(data = all.CC, aes(CC, color = method))

DL <- read_pickle_file("cc_pancancer_30L_77e_val_mse.pkl")
names(DL) <- 'CC'
all.CC <- rbind(all.CC, data.table(DL, method = 'DL_val_mse'))
ggplot() + geom_density(data = all.CC, aes(CC, color = method))

DL <- read_pickle_file("final.pkl")
names(DL) <- 'CC'
all.CC <- rbind(all.CC, data.table(DL, method = 'DL_final'))
ggplot() + geom_density(data = all.CC, aes(CC, color = method))
```

### ML and DL on 2672 circRNAs

``` r
cv.full.40.pearson <- list(
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', cutoff = 40, method = 'KNN'),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', cutoff = 40, method = 'RF', ncores = 56),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', cutoff = 40, method = 'EN', ncores = 56),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', cutoff = 40, method = 'lasso', ncores = 56),
  ICE_cv_entire(train.pcg = mionco.pcg, train.circ = mionco.circ, cor.method = 'pearson', cutoff = 40, method = 'SVM', ncores = 56))
names(cv.full.40.pearson) <- c('KNN', 'RF', 'EN', 'lasso', 'SVM')
cv.full.40.pearson <- lapply(names(cv.full.40.pearson), function(x) data.table(cv.full.40.pearson[[x]], method = x))
cv.full.40.pearson <- do.call(rbind, cv.full.40.pearson)
DL <- read_pickle_file("final_40.pkl")
names(DL) <- 'CC'
all.CC.40 <- rbind(cv.full.40.pearson[, .(CC, method)], data.table(DL, method = 'DL_40'))
ggplot() + geom_density(data = all.CC, aes(CC, color = method))
```

## Example dataset

### Pan-cancer dataset

To perform further analysis on circRNAs, we use convert genome region
(**Chromosome**:**Start**-**End**) as circRNA identifiers.

One pair of circRNA and PCGs matrices is from paper [*The Landscape of
Circular RNA in
Cancer*](https://www.sciencedirect.com/science/article/pii/S0092867418316350).

circRNA expression matrix was downloaded from [this
website](https://mioncocirc.github.io/download/) with the metadata and
pre-processed by below code:

``` r
library(data.table)
mionco.circ <- fread('v0.1.release.txt')[, `:=`(symbol = NULL, release = NULL, ID = paste0(chr, ':', start, '-', end))]
mionco.circ <- dcast(mionco.circ, formula = sample ~ ID, value.var = 'reads', fun.aggregate = max)
mionco.circ <- as.matrix(mionco.circ, rownames = 'sample')

mionco.pcg <- fread('fpkm_matrix.csv', drop = 1)
mionco.pcg[, Row.names := tstrsplit(Row.names, '\\.')[[1]]]
mionco.pcg[, V2 := NULL]

library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
pcg.id <- getBM(attributes = 'ensembl_gene_id', filters = "biotype", values = "protein_coding", mart = mart)
setnames(mionco.pcg, 'Row.names', 'ensembl_gene_id')
mionco.pcg <- mionco.pcg[pcg.id, on = 'ensembl_gene_id', nomatch = NULL]
mionco.pcg <- as.matrix(mionco.pcg, rownames = 'ensembl_gene_id')
mionco.pcg <- t(mionco.pcg)

sample.used <- intersect(rownames(mionco.circ), rownames(mionco.pcg))
mionco.circ <- mionco.circ[sample.used, ]
mionco.pcg <- mionco.pcg[sample.used, ]
```

### CCLE dataset

The CCLE circRNAs expression was downloaded from [CircRic
database](https://hanlab.uth.edu/cRic/download), while the RNAseq gene
expression as well as sample (cell line) annotation was downloaded from
[CCLE portal of The Broad
Institute](https://portals.broadinstitute.org/ccle/data). As CCLE data
is processed with hg19, we also need prepare a [chain
file](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz)
for
liftOver.

``` r
ccle.circ <- fread('circRNA_expression.csv')[, `:=`(ccl_id = NULL, cancer_type = NULL, circ_gene = NULL)]
anno <- fread('Cell_lines_annotations_20181226.txt')
anno <- anno[, .(CCLE_ID, Name)]
ccle.circ <- merge(ccle.circ, anno, by.x = 'ccl_name', by.y = 'Name')[, ccl_name := NULL]
splitted <- tstrsplit(ccle.circ[[1]], '_|\\|')
temp.19 <- GRanges(seqnames = paste0('chr', splitted[[1]]), ranges = IRanges(start = as.numeric(splitted[[2]]), end = as.numeric(splitted[[3]])), name = ccle.circ[['CCLE_ID']], count = ccle.circ[['number_bsReads']])
chain <- import.chain("hg19ToHg38.over.chain")
temp.38 <- liftOver(temp.19, chain)
temp.38 <- as.data.table(unlist(temp.38))
ccle.circ <- data.table(ID = paste0(temp.38[['seqnames']], ':', temp.38[['start']], '-', temp.38[['end']]), name = temp.38[['name']], count = temp.38[['count']])
ccle.circ <- dcast(ccle.circ, formula = name ~ ID, value.var = 'count', fun.aggregate = max)
ccle.circ <- as.matrix(ccle.circ, rownames = 'name')

ccle.pcg <- fread('CCLE_RNAseq_rsem_genes_tpm_20180929.txt')
ccle.pcg[, `:=`(transcript_ids = NULL, gene_id = tstrsplit(gene_id, '\\.')[[1]])]
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
pcg.id <- getBM(attributes = 'ensembl_gene_id', filters = "biotype", values = "protein_coding", mart = mart)
setnames(ccle.pcg, 'gene_id', 'ensembl_gene_id')
ccle.pcg <- ccle.pcg[pcg.id, on = 'ensembl_gene_id', nomatch = NULL]
ccle.pcg <- as.matrix(ccle.pcg, rownames = 'ensembl_gene_id')
ccle.pcg <- t(ccle.pcg)

sample.used <- intersect(rownames(ccle.circ), rownames(ccle.pcg))
ccle.circ <- ccle.circ[sample.used, ]
ccle.pcg <- ccle.pcg[sample.used, ]
```

Up to now, `mionco.circ` and `mionco.pcg` could be used as circRNA and
PCG dataset for training or testing, and the `ccle.circ` and `ccle.pcg`
ditto.

## Cross-dataset validation

``` r
# Keep circRNAs expressed in at least 1% sample
# 非常惨QAQ
train.circ <- ICE:::filter_circ(mionco.circ, 1)
new.circ <- ICE:::filter_circ(ccle.circ, 1)
train.pcg <- ICE:::pre_process(mionco.pcg)
new.pcg <- ICE:::pre_process(ccle.pcg)
# 取完交集就剩一点点了
used.circ.name <- intersect(colnames(new.circ), colnames(train.circ))

fwrite(train.circ[, used.circ.name], 'train_circ.csv')
fwrite(new.circ[, used.circ.name], 'new_circ.csv')
fwrite(train.pcg, 'train_pcg.csv')
fwrite(new.pcg, 'new_pcg.csv')

sfInit(parallel = TRUE, cpu = 56)
cl <- sfGetCluster()
sfExport('train.pcg', 'train.circ', 'new.pcg')
sfExport('ICE', namespace = 'ICE')
pboptions(use_lb = TRUE)

predicted <- pblapply(used.circ.name, function(x) ICE(train.pcg, train.circ, new.pcg, x, method = 'RF', filter = FALSE), cl = cl)
names(predicted) <- used.circ.name
cd.cc <- unlist(lapply(used.circ.name, function(x) cor(predicted[[x]], ccle.circ[, x], method = 'spearman')))
```

``` r
meta <- fread('meta_update.csv')
meta[, count := .N, by = 'Analysis Cohort']
head(unique(meta, by = 'Analysis Cohort')[order(-count)])

fwrite(ICE:::filter_circ(mionco.circ), 'mionco_circ.csv')
fwrite(ICE:::pre_process(mionco.pcg), 'mionco_pcg.csv')

cohort_cor <- function(x){
  cohort.sample.used <- intersect(sample.used, meta[x, on = 'Analysis Cohort', ID])
  if(length(cohort.sample.used) == 0) return()
  circ.cohort <- mionco.circ[cohort.sample.used, ]
  pcg.cohort <- mionco.pcg[cohort.sample.used, ]
  
  cv <- list(
  ICE_cv_entire(train.pcg = pcg.cohort, train.circ = circ.cohort, cor.method = 'pearson', folds = 3, method = 'KNN'),
  ICE_cv_entire(train.pcg = pcg.cohort, train.circ = circ.cohort, cor.method = 'pearson', folds = 3, method = 'RF'))
  names(cv) <- c('KNN', 'RF')
  cv <- lapply(names(cv), function(x) data.table(cv[[x]], method = x))
  cv <- do.call(rbind, cv)
  (data.table(cv, cohort = x))
}

cv.full <- lapply(head(unique(meta, by = 'Analysis Cohort')[order(-count)])[['Analysis Cohort']], function(x) cohort_cor(x))

names(cv.full) <- head(unique(meta, by = 'Analysis Cohort')[order(-count)])[['Analysis Cohort']]

cv.full <- Filter(function(x) !is.null(x), cv.full)
cv.full <- do.call(rbind, cv.full)
ggplot() + geom_density(data = cv.full[method == 'RF' & cohort == 'PRAD'], aes(CC))
d_fun <- ecdf (cv.full[method == 'RF' & cohort == 'PRAD', CC])
ggplot(data = cv.full, aes(x = cohort, y = CC, fill = method)) + geom_boxplot() + geom_jitter(aes(fill = method), shape=16, position=position_jitterdodge(0.2))
```

## m6A

``` r
library(data.table)
ip <- fread('rpkm_ip_norm_m.txt')
anno <- fread('fix.anno.txt')
```

## Tuning for static ML parameters

``` r
library(mlr3)
library(data.table)
# install.packages("mlr3learners.randomforest", repos = "https://mlr3learners.github.io/mlr3learners.drat")
library(mlr3learners.randomforest)
library(paradox)
library(mlr3tuning)

train.circ <- scale(ICE:::filter_circ(mionco.circ))
train.pcg <- ICE:::pre_process(mionco.pcg)
train.pcg <- scale(ICE:::corf(train.pcg, train.circ, "chr7:23611170-23611553", 50))

task.circ <- TaskRegr$new(id = 'circ', backend = data.table(target = train.circ[, "chr7:23611170-23611553"], train.pcg), target = "target")
```

### Tuning for random forest

``` r
learner.circ <- lrn('regr.randomForest')
# learner.circ$param_set
tune.ps <- ParamSet$new(list(
  ParamInt$new('ntree', lower = 1, upper = 500),
  ParamInt$new('mtry', lower = 1, upper = 33),
  ParamInt$new('nodesize', lower = 1, upper = 50)
))

instance <- TuningInstance$new(
  task = task.circ,
  learner = learner.circ,
  resampling = rsmp('cv'),
  measures = msr('regr.mse'),
  param_set = tune.ps,
  terminator = term("stagnation", threshold = 1e-5)
)

tuner <- tnr("grid_search", resolution = 10)

result <- tuner$tune(instance)
```

### Tuning for SVM

``` r
learner.circ <- lrn('regr.svm')
# learner.circ$param_set
tune.ps <- ParamSet$new(list(
  ParamDbl$new('cost', lower = 0.0001, upper = 10000),
  ParamDbl$new('gamma', lower = 0.01, upper = 1),
  ParamFct$new('type', 'eps-regression'),
  ParamFct$new('kernel', 'radial')
))

instance <- TuningInstance$new(
  task = task.circ,
  learner = learner.circ,
  resampling = rsmp('cv'),
  measures = msr('regr.mse'),
  param_set = tune.ps,
  terminator = term("stagnation", threshold = 1e-5)
)

tuner <- tnr("grid_search", resolution = 10)

result <- tuner$tune(instance)
```

### Tuning for KNN

``` r
learner.circ <- lrn('regr.fnn')
# learner.circ$param_set
tune.ps <- ParamSet$new(list(
  ParamInt$new('k', lower = 1, upper = task.circ$nrow - 100)
))

instance <- TuningInstance$new(
  task = task.circ,
  learner = learner.circ,
  resampling = rsmp('cv'),
  measures = msr('regr.mse'),
  param_set = tune.ps,
  terminator = term("stagnation", threshold = 1e-5)
)

tuner <- tnr("grid_search", resolution = 10)

result <- tuner$tune(instance)
```

### Tuning for lasso regression

``` r
library(glmnet)
cv.lasso <- cv.glmnet(train.pcg, train.circ[, "chr7:23611170-23611553"], alpha = 1, family = 'gaussian')
cv.lasso$lambda.min
```

### Tuning for Elastic Net

``` r
library(caret)
en.model <- suppressWarnings(train(
        train.pcg, train.circ[, "chr7:23611170-23611553"], method = 'glmnet',
        trControl = trainControl("cv", number = 10),
        tuneLength = 100))
en.model$bestTune
```
