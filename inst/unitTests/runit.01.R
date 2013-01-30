test.read.nifti <- function() {
    datadir <- system.file("data", package="niftir")
    testfile1 <- file.path(datadir, "test.nii.gz")
    testfile2 <- file.path(datadir, "test.rda")
    
    tryCatch(
        nim1 <- read.nifti(testfile1),
        error=function(e) checkTrue(FALSE, "error reading nifti file")
    )
    load(testfile2)
    
    checkIdentical(nim1, nim2, "nifti file read incorrectly")   
}

test.write.nifti <- function() {
    datadir <- system.file("data", package="niftir")
    testfile1 <- file.path(datadir, "test.rda")
    testfile2 <- file.path(datadir, "test2.nii.gz")
    
    load(testfile1)
    
    tryCatch(
        write.nifti(nim2, testfile2),
        error=function(e) checkTrue(FALSE, "error writing nifti file")
    )
    
    nim1 <- read.nifti(testfile2)
    
    checkIdentical(nim2, nim1, "nifti file not written correctly")
}
