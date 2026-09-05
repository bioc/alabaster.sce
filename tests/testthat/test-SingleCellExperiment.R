# library(testthat); library(alabaster.sce); source("test-SingleCellExperiment.R")

# Making an SE and annotating it.
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

se <- SingleCellExperiment(list(counts=mat, cpm=mat/10))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))
rownames(se) <- sprintf("GENE_%i", seq_len(nrow(se)))

test_that("saveObject works as expected for SCE objects", {
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_s4_class(out2, "SingleCellExperiment")
    expect_identical(rowData(out2), rowData(se))
    expect_identical(colData(out2), colData(se))
})

test_that("saveObject works with some non-trivial rowRanges", {
    all.ranges <- GRanges("chrX", IRanges(seq_len(nrow(mat) * 2), width=1))
    rowRanges(se) <- splitAsList(all.ranges, rep(seq_len(nrow(mat)), length.out=length(all.ranges)))

    tmp <- tempfile()
    saveObject(se, tmp, "rnaseq")
    out2 <- readObject(tmp)

    expect_identical(rowData(out2), rowData(se))
    expect_identical(rowRanges(out2), rowRanges(se))
    expect_identical(rownames(out2), rownames(se))
})

test_that("saveObject works as expected with reduced dims inside", {
    reducedDims(se) <- list(PCA=matrix(rnorm(ncol(mat)*50), ncol=50), TSNE=cbind(TSNE1=runif(ncol(mat)), TSNE2=runif(ncol(mat))))

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(reducedDimNames(out2), reducedDimNames(se))
    expect_identical(as.matrix(reducedDim(out2, "PCA")), reducedDim(se, "PCA"))
    expect_identical(as.matrix(reducedDim(out2, "TSNE")), reducedDim(se, "TSNE"))
})

test_that("saveObject works as expected with alternative experiments inside", {
    copy <- se
    altExps(se) <- list(spikes=copy[1:2,], protein=copy[3:5,])

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(altExpNames(out2), altExpNames(se))
    expect_identical(rownames(altExp(out2, 1)), rownames(altExp(se, 1)))
    expect_identical(rownames(altExp(out2, 2)), rownames(altExp(se, 2)))
})

test_that("saveObject works as expected when we slap in a main name", {
    mainExpName(se) <- "FOO"

    tmp <- tempfile()
    saveObject(se, tmp)
    out <- readObject(tmp)
    expect_identical(mainExpName(out), "FOO")

    # Fails if there's an alternative experiment with the same name.
    altExp(se, "FOO") <- se[1:10,]
    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "expected 'main_experiment_name'")
})
