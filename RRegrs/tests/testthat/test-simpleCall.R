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

test_that("the method doesn't fail for just MLR and PLS", {
  dsData = system.file("extdata", "ds.House.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="T",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=5,CVtypes="repeatedcv"
  )
})

test_that("the method works without Y-randomization", {
  dsData = system.file("extdata", "ds.House.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="T",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv"
  )
})

test_that("the method works for ds.gajewicz.csv", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="T",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=5,CVtypes="repeatedcv"
  )
})

test_that("the method works with both repeatedcv and LOOCV", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="T",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})

#test_that("the method works with GLM", {
#  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
#  results = RRegrs(
#    DataFileName=dsData,
#    fLM="T",fGLM="T",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
#    fRFRFE="F",fSVMRFE="F",fENET="F",
#    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
#  )
#})

test_that("the method works with LASSO", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="T",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})


test_that("the method works with SVRM", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="T",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})

test_that("the method works with NN", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="T",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})

test_that("the method works with RF", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="T",
    fRFRFE="F",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})

test_that("the method works with RFRFE", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="T",fSVMRFE="F",fENET="F",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})

#test_that("the method works with SVRMRFE", {
#  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
#  results = RRegrs(
#    DataFileName=dsData,
#    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
#    fRFRFE="F",fSVMRFE="T",fENET="F",
#    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
#  )
#})

test_that("the method works with ENET", {
  dsData = system.file("extdata", "ds.gajewicz.csv", package = "RRegrs")
  results = RRegrs(
    DataFileName=dsData,
    fLM="T",fGLM="F",fPLS="F",fLASSO="F",fSVRM="F",fNN="F",fRF="F",
    fRFRFE="F",fSVMRFE="F",fENET="T",
    trainFrac=0.75,iSplitTimes=3,noYrand=0,CVtypes="repeatedcv;LOOCV"
  )
})
