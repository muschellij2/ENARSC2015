## ----knit-setup, echo=FALSE, results='hide', eval=TRUE-------------------
rm(list=ls())
library(knitr)
enar_dir = '~/ENARSC2015'
enar_dir = Sys.readlink(enar_dir)
opts_chunk$set(comment="")
# cache.path <- path.expand(file.path(enar_dir, "/vignettes/tutorial_cache/"))
# opts_chunk$set(cache.path=cache.path)

## ------------------------------------------------------------------------
## library(fslr)
library(rgl)
library(misc3d)
library(brainR)
library(scales)
library(ENARSC2015)
library(brainR)
library(xtable)
library(plyr)

## ----kki_readin----------------------------------------------------------

##Baseline data##
enar_dir = system.file(package="ENARSC2015")
mridir = file.path(enar_dir, "Kirby21")
MPRAGE_base <- readNIfTI(
  file.path(mridir, 'SUBJ0001-01-MPRAGE.nii.gz'), 
  reorient=FALSE)

## ------------------------------------------------------------------------
mask <- readNIfTI(
  file.path(mridir, 'SUBJ0001_mask.nii.gz'), 
  reorient = FALSE)

## ------------------------------------------------------------------------
MPRAGE_follow <- readNIfTI(
  file.path(mridir, 'SUBJ0001-02-MPRAGE.nii.gz'),
  reorient=FALSE)

## ------------------------------------------------------------------------
dim(MPRAGE_base)
dim(MPRAGE_follow)

## ----kki_plotimage, cache=TRUE-------------------------------------------
##axial slice##
image(MPRAGE_base[,,128])
image(MPRAGE_base, z = 128, plot.type = "single")
image(mask[,,128])

##coronal slice##
image(MPRAGE_base[,128,], col = rainbow(12))
image(MPRAGE_base[,128,], col = gray(0:64/64))

##sagittal slice##
image(MPRAGE_base[85,,], col = topo.colors(12))
image(MPRAGE_base[85,,], col = gray(0:64/64))

## ------------------------------------------------------------------------
orthographic(MPRAGE_base)
orthographic(MPRAGE_base, xyz = c(90, 100, 15))

## ------------------------------------------------------------------------
mean(MPRAGE_base[mask == 1])
sd(MPRAGE_base[mask ==1])

## ------------------------------------------------------------------------
MPRAGE_diff <- MPRAGE_base - MPRAGE_follow
image(MPRAGE_diff[,,128])

## ------------------------------------------------------------------------
MPRAGE_base_n3 <- readNIfTI(
  file.path(mridir, 'SUBJ0001-01-MPRAGE_N3.nii.gz'),
  reorient=FALSE)

## ------------------------------------------------------------------------
MPRAGE_follow_n3_reg_base <- readNIfTI(
  file.path(mridir, 'SUBJ0001-02-MPRAGE_N3_REG.nii.gz'),
  reorient=FALSE)

## ----kki_imginfo---------------------------------------------------------
MPRAGE_base
slotNames(MPRAGE_base)

## ----kki_imgdiff, cache=TRUE---------------------------------------------
difference <- MPRAGE_follow_n3_reg_base - MPRAGE_base 
class(difference)
difference = niftiarr(MPRAGE_follow_n3_reg_base, difference)
image(difference[,,128])
image(difference, z= 128, plot.type="single")

## ----kki_n3vsOrig, cache=TRUE--------------------------------------------
sample <- sample(1:sum(mask), 4000, replace = FALSE)

plot(MPRAGE_base[mask ==1][sample], MPRAGE_base_n3[mask == 1][sample], 
     main = "N3 Corrected Intensities", xlab = "Raw MPRAGE Intensities", 
     ylab = "N3 Corrected Intensities", pch = 20)
abline(0,1,lty = 2, col = "red") 

## ------------------------------------------------------------------------
plot(MPRAGE_base_n3[mask ==1][sample], MPRAGE_follow_n3_reg_base[mask == 1]
     [sample], main = "Registered Image", xlab = "Baseline", 
     ylab = "Follow Up", pch = 20)
abline(0,1,lty = 2, col = "red") 

## ----kki_rglexample, cache=TRUE------------------------------------------
template <- readNIfTI(
  file.path(enar_dir, "Template", "MNI152_T1_2mm_brain.nii.gz"), 
  reorient=FALSE)
orthographic(template)

## ------------------------------------------------------------------------
cut <- 4500
dtemp <- dim(template)
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = cut,
          alpha = 0.1, draw = TRUE, color="black")

text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*0.98, text="Top")
text3d(x=-1, y=dtemp[2]/2, z = dtemp[3]/2, text="Right")

roi = cal_img(template > 8100)
contour3d(roi, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = .99,
          alpha = 0.5, draw = TRUE, add=TRUE, color = "red")
index <- writeWebGL_split(
  dir=file.path(getwd(), "WebGL"), 
  width=700, height=500, 
  template= system.file("my_template.html", package="brainR"))


## ----, kki_blur, cache=TRUE----------------------------------------------

library(AnalyzeFMRI)
sigma.smooth<-diag(3,3)
k.size<- 7
MPRAGE_base_smoothed <-GaussSmoothArray(MPRAGE_base,sigma=sigma.smooth,
                                        ksize=k.size,mask=mask)

## ------------------------------------------------------------------------
orthographic(MPRAGE_base_smoothed)
orthographic(MPRAGE_base * mask)

## ------------------------------------------------------------------------
writeNIfTI(MPRAGE_base_smoothed, filename = 
             "MPRAGE_base_smoothed", verbose = TRUE, gzipped = TRUE)

## ----n3_compare, cache=FALSE, results = 'hide'---------------------------
ext = ".nii.gz"
mods = c("T1", "T2", "FLAIR")
pipes = c("FSL", "ANTsR")

filedir = system.file("FSL", package="ENARSC2015")

# 
# n3.files = file.path(filedir, 
#                      paste0(mods, "_FSL_N3Correct", ext))
# ex = file.exists(n3.files)
# stopifnot(all(ex))  
# names(n3.files) = mods
# 
# n3.imgs = lapply(n3.files, readNIfTI, reorient=FALSE)


## ----reg_compare, cache=TRUE, results='hide'-----------------------------
reg.files = c(FLAIR = file.path(filedir, paste0("FLAIR_FSL_N3Correct_regMNI_FNIRT", ext)), 
          ROI = file.path(filedir, paste0("ROI_threshold_regMNI_FNIRT", ext)))
          

reg.imgs = lapply(reg.files, readNIfTI, reorient=FALSE)


## ----plot_reg, cache=TRUE------------------------------------------------

fslROI = reg.imgs$ROI
fslROI[ fslROI == 0 ] = NA
orthographic(reg.imgs$FLAIR, fslROI,
             col.y=  alpha("red", .25),             
             text = "Registered FLAIR \n FSL FNIRT")

## ----atlas_list, cache=TRUE, results='hide'------------------------------
### Getting Harvard-Oxford Subcortical Atlas
data(Atlas_Labels)

### returns a list of indices from the template which match a certiain label
get.ind = function(img, df){
  df = df[order(df$index),]
  ind.list = llply(df$index, function(x){
  	return(which(img %in% x))
	}, .progress = "text")
	names(ind.list) = df$Label
	return(ind.list)
}

hoxsubcort.list = get.ind(hoxsubcort.img, hoxsubcort.df)

## ----atlas_breakdown, cache=TRUE, results='hide'-------------------------
### in each area, takes the sum of an image (in this case the ROI)
tab.area = function(img, ind.list, keepall = FALSE) {
  ## get overlap of indices
  raw.mat = sapply(ind.list, function(x) sum(img[x]))
  names(raw.mat) = names(ind.list)
  ## cs is sum of indices of overlap
  cs.raw = data.frame(N_Voxels=raw.mat) 
  rownames(cs.raw) = names(ind.list)
  if (!keepall) cs.raw = cs.raw[cs.raw != 0, , drop=FALSE]
  return(cs.raw)
}

fslROI = reg.imgs$ROI
fslArea = tab.area(fslROI, hoxsubcort.list)
fslArea$cm3 = fslArea$N_Voxels * prod(pixdim(reg.imgs$ROI)[2:4])/1000

## ----, results='asis', echo=FALSE----------------------------------------

cat("### FSL ROI Breakdown\n")
xtab = xtable(fslArea)
print.xtable( xtab, type="html")


## ----plot_atlas, cache = TRUE--------------------------------------------
# xyz = ceiling(cog(fslROI))
xyz = c(134, 120, 85)
res = hoxsubcort.img
res[ fslROI == 0 ] = NA
res[ res %in% -99] = 0
res = cal_img(res)
n = res@cal_max
bcols <- brewer_pal("div", pal = "RdBu")(5)
cols = gradient_n_pal(bcols)(seq(0, 1, length = n))
cols = alpha(cols, 0.25)
orthographic(reg.imgs$FLAIR, res, 
             xyz=xyz, col.y= cols)

