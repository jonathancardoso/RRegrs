library(RRegrs)
context("simple")

test_that("it fails when only one modeling method is selected", {
  dsData = system.file("extdata", "ds.House.csv", package = "RRegrs")
  expect_error(
    RRegrs(
      DataFileName=dsData,
      fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
      fRFRFE="F",fSVMRFE="F",fENET="F",
      trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="LOOCV"
    ),
    "You must select at least two modelling methods to compare."
  )
})

test_that("the chemical classics data can be loaded", {
  chemClassicsData = system.file("extdata", "chemClassic.tsv", package = "RRegrs")
  data = read.table(chemClassicsData, sep="\t", header=TRUE)
  expect_equal(9, nrow(data))
})

test_that("the method doesn't fail for just MLR", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=49,CVtypes="LOOCV"
  )
})

