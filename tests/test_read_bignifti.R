 library(niftir)
library(testthat)
library(stringr)

context("Read 4D File as a 2D Big Matrix")

#infile <- system.file("data/test_func.nii.gz", package="niftir")
#hdr <- read.nifti.header(infile)
#ref.dim <- c(hdr$dim[4], prod(hdr$dim[1:3]))
#ref.arr <- as.integer(array(1:prod(hdr$dim), hdr$dim))
#outfile <- file.path(getwd(), "pkg/niftir/data/test_func_inds.nii.gz")
#write.nifti(ref.arr, hdr, outfile=outfile, odt="int")
#    
## Save nifti data into the big matrix
#rows <- as.double(1:idim[1]); cols <- as.double(1:idim[2])
#.Call("read_bignifti_data2", x$address, bigmat@address, rows, cols)

test_that("Create appropriately sized big matrix", {
    # Setup
    funcfile <- system.file("data/test_func_inds.nii.gz", package="niftir")
    nifti_obj <- read.big.nifti_obj(funcfile)
    
    # Create big matrix without setting anything
    bigmat <- read.big.nifti_gen.bigmat(nifti_obj)
    expect_that(dim(bigmat), equals(c(5, 24)))
    
    # Create big matrix with default sizes
    tpts <- as.double(1:nifti_obj$header$dim[4])
    voxs <- as.double(1:prod(nifti_obj$header$dim[-4]))
    bigmat <- read.big.nifti_gen.bigmat(nifti_obj, tpts, voxs)
    expect_that(dim(bigmat), equals(c(5, 24)))
    
    # Create big matrix as parts
    tpts <- as.double(1:2)
    voxs <- as.double(1:4)
    bigmat <- read.big.nifti_gen.bigmat(nifti_obj, tpts, voxs)
    expect_that(dim(bigmat), equals(c(2, 4)))
})

test_that("Reads in nifti file", {
    # Setup
    funcfile <- system.file("data/test_func_inds.nii.gz", package="niftir")
    nifti_obj <- read.big.nifti_obj(funcfile)
    hdr <- nifti_obj$header
    
    ## Read in full nifti file
    # Reference
    ref.dim <- c(hdr$dim[4], prod(hdr$dim[1:3]))
    ref.mat <- read.nifti.image(funcfile)
    dim(ref.mat) <- rev(ref.dim)
    ref.mat <- t(ref.mat)
    # Comparison
    ref.bigmat <- read.big.nifti_gen.bigmat(nifti_obj)
    ref.bigmat <- read.big.nifti_read(nifti_obj, ref.bigmat)
    # Compare
    expect_that(ref.mat[,], is_equivalent_to(ref.bigmat[,]))
    
    ## Read in partial nifti file
    tpts <- as.double(1:2)
    voxs <- as.double(1:4)
    # Reference
    ref.dim <- c(hdr$dim[4], prod(hdr$dim[1:3]))
    ref.mat <- read.nifti.image(funcfile)
    dim(ref.mat) <- rev(ref.dim)
    ref.mat <- t(ref.mat[voxs,tpts])
    # Comparison
    ref.bigmat <- read.big.nifti_gen.bigmat(nifti_obj, tpts, voxs)
    ref.bigmat <- read.big.nifti_read(nifti_obj, ref.bigmat)
    # Compare
    expect_that(ref.mat[,], is_equivalent_to(ref.bigmat[,]))
    
    ## Read in partial nifti file
    tpts <- as.double(2:3)
    voxs <- as.double(2:5)
    # Reference
    ref.dim <- c(hdr$dim[4], prod(hdr$dim[1:3]))
    ref.mat <- read.nifti.image(funcfile)
    dim(ref.mat) <- rev(ref.dim)
    ref.mat <- t(ref.mat[voxs,tpts])
    # Comparison
    ref.bigmat <- read.big.nifti_gen.bigmat(nifti_obj, tpts, voxs)
    ref.bigmat <- read.big.nifti_read(nifti_obj, ref.bigmat)
    # Compare
    expect_that(ref.mat[,], is_equivalent_to(ref.bigmat[,]))
})

test_that("Read in nifti as big nifti", {
    funcfile <- system.file("data/test_func_inds.nii.gz", package="niftir")
    
    ###
    # Full File
    ###
    
    # Reference
    ref.mat <- array(1:prod(2:5), 2:5)
    dim(ref.mat) <- c(24,5)
    ref.mat <- t(ref.mat)
    ref <- as.big.matrix(ref.mat)
    
    # Comparison
    comp <- read.big.nifti(funcfile)
    
    expect_that(ref[,], equals(comp[,]))

    ###
    # Partial File
    ###
    
    # Reference
    ref.mat <- array(1:prod(2:5), 2:5)
    dim(ref.mat) <- c(24,5)
    ref.mat <- t(ref.mat[2:5,2:3])
    ref <- as.big.matrix(ref.mat)
    
    # Comparison
    comp <- read.big.nifti(funcfile, tpts=2:3, voxs=2:5)
    
    expect_that(ref[,], equals(comp[,]))
})
