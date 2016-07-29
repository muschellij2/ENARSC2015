
## ----knit-setup, echo=FALSE, results='hide', eval=TRUE-------------------
rm(list=ls())
library(knitr)
opts_chunk$set(cache = TRUE)


## ----, eval = FALSE, echo = TRUE-----------------------------------------
## library(devtools)
## install_github('stnava/ANTsR')
## install_github('muschellij2/ENARSC2015')


## ------------------------------------------------------------------------
library(ANTsR)
library(ENARSC2015)


## ----set_enar_dir--------------------------------------------------------
rootdir = "~/ENARSC2015"
rootdir = Sys.readlink(rootdir)
outdir <- file.path(rootdir, 'inst', 'ANTsR')
stopifnot(file.exists(outdir))


## ----get_img_names-------------------------------------------------------
t1.path <- system.file("NIfTI/T1.nii.gz", package = 'ENARSC2015')
t2.path <- system.file("NIfTI/T2.nii.gz", package = 'ENARSC2015')
flair.path <- system.file("NIfTI/FLAIR.nii.gz", package = 'ENARSC2015')


## ----read_imgs-----------------------------------------------------------
t1 <- antsImageRead(t1.path, 3)
t2 <- antsImageRead(t2.path, 3)
flair <- antsImageRead(flair.path, 3)


## ----plot_imgs, fig.width= 10, fig.height = 3----------------------------
par(mfrow = c(1,3))
image(as.array(t1)[,,10], col = gray(0:64/64))
image(as.array(t2)[,,10], col = gray(0:64/64))
image(as.array(flair)[,,10], col = gray(0:64/64))


## ----clone_imgs----------------------------------------------------------
t1N3 <- antsImageClone(t1)
t2N3 <- antsImageClone(t2)
flairN3 <- antsImageClone(flair)


## ----n3_correct----------------------------------------------------------
N3BiasFieldCorrection(t1@dimension, t1, t1N3, "4")
N3BiasFieldCorrection(t2@dimension, t2, t2N3, "4")
N3BiasFieldCorrection(flair@dimension, flair, flairN3, "4")


## ----antsWrite-----------------------------------------------------------
antsImageWrite(t1N3, file.path(outdir, 'T1_ANTsR_N3Correct.nii.gz'))
antsImageWrite(t2N3, file.path(outdir,'T2_ANTsR_N3Correct.nii.gz'))
antsImageWrite(flairN3, file.path(outdir,'FLAIR_ANTsR_N3Correct.nii.gz'))


## ----plot_n3_imgs, fig.width= 10, fig.height = 3-------------------------
par(mfrow = c(1,3))
image(as.array(t1N3)[,,10], col = gray(0:64/64))
image(as.array(t2N3)[,,10], col = gray(0:64/64))
image(as.array(flairN3)[,,10], col = gray(0:64/64))


## ----read_template-------------------------------------------------------
template.path <- system.file("MNI152_T1_1mm_brain.nii.gz", package = 'ENARSC2015')
template.skull.path <- system.file("MNI152_T1_1mm.nii.gz", package = 'ENARSC2015')
template <- antsImageRead(template.path, 3)
template.skull <- antsImageRead(template.skull.path, 3)


## ----plot_template,  fig.width= 6.5, fig.height = 3----------------------
par(mfrow = c(1,2))
image(as.array(template)[,,90], col = gray(0:64/64))
image(as.array(template.skull)[,,90], col = gray(0:64/64))



## ----register_rigid------------------------------------------------------
antsRegOut <- antsRegistration(fixed = flairN3, moving = t1N3 , typeofTransform = "Rigid",  outprefix = "./test")
t1.to.flair <-antsImageClone(antsRegOut$warpedmovout)


## ----write_rigid---------------------------------------------------------
antsImageWrite(t1.to.flair, file.path(outdir, 'T1_ANTsR_N3Correct_regFLAIR.nii.gz'))


## ----plot_rigid, fig.width= 6.5, fig.height = 10-------------------------
par(mfrow = c(3,2))
slices <- c(5, 10, 15)
for(i in slices){
image(as.array(flairN3)[,,i], col = gray(0:64/64))
image(as.array(t1.to.flair)[,,i], col = gray(0:64/64))
}


## ----subtract_image------------------------------------------------------
subtraction.image <- antsImageClone(flairN3)
ImageMath( 3 , subtraction.image  , "-", flairN3 , t1.to.flair )
image(as.array(subtraction.image)[,,10], col = gray(0:64/128))


## ----roi_img-------------------------------------------------------------
roi.path <- system.file("NIfTI/ROI.nii.gz", package = 'ENARSC2015')
roi <- antsImageRead(roi.path, 3)


## ----plot_roi------------------------------------------------------------
par(mfrow = c(1,2))
image(as.array(flairN3)[,,15], col = gray(0:64/64))
image(as.array(roi)[,,15], col = gray(0:64/64))


## ----register_syn--------------------------------------------------------
outprefix = file.path(rootdir, "inst", 'ANTsR', "ants")
antsRegOut.nonlin <- antsRegistration(fixed = template.skull, moving = t1.to.flair, typeofTransform = "SyN",  outprefix = outprefix)


## ----clone_t1------------------------------------------------------------
t1.to.flair.to.template <-antsImageClone(antsRegOut.nonlin$warpedmovout)


## ----warp_flair----------------------------------------------------------
flair.to.template <- antsApplyTransforms(fixed=template.skull , moving=flairN3 , transformlist=antsRegOut.nonlin$fwdtransforms , interpolator="Linear")
roi.to.template <-antsApplyTransforms(fixed=template.skull , moving=roi , transformlist=antsRegOut.nonlin$fwdtransforms , interpolator="Linear")


## ----warp_roi------------------------------------------------------------
roi.to.template.threshold <-antsImageClone(roi.to.template)
ThresholdImage(3, roi.to.template, roi.to.template.threshold, .5, 1)


## ----plot_warped_imgs, fig.width= 10, fig.height = 3---------------------
par(mfrow = c(1,3))
image(as.array(flair.to.template)[,,110], col = gray(0:64/64))
image(as.array(roi.to.template.threshold)[,,110], col = gray(0:64/64))
image(as.array(template.skull)[,,110], col = gray(0:64/64))


## ----write_warped_imgs---------------------------------------------------
antsImageWrite(flair.to.template, file.path(outdir, 'FLAIR_ANTsR_N3Correct_regMNI_SyN.nii.gz'))
antsImageWrite(roi.to.template , file.path(outdir, 'ROI_regMNI_SyN.nii.gz'))
antsImageWrite(roi.to.template.threshold, file.path(outdir, 'ROI_threshold_regMNI_SyN.nii.gz'))


## ----remove_warp_img-----------------------------------------------------
file.remove(paste0(outprefix, "1InverseWarp.nii.gz" ) )
file.remove(paste0(outprefix, "1Warp.nii.gz" ) )

