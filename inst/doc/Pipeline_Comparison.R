## ----knit-setup, echo=FALSE, results='hide', eval=TRUE-------------------
rm(list=ls())
library(knitr)
enar_dir = '~/ENARSC2015'
enar_dir = Sys.readlink(enar_dir)
opts_chunk$set(echo=TRUE, prompt=FALSE, message=TRUE, warning=TRUE, comment="")
# cache.path <- path.expand(file.path(enar_dir, "/vignettes/compare_cache/"))
# opts_chunk$set(cache.path=cache.path)

## ------------------------------------------------------------------------
library(fslr)
library(plyr)
library(scales)
library(xtable)
# outdir = "~/Desktop/results"
# rootdir = '~/ENARSC2015'
# rootdir = Sys.readlink(rootdir)
# outdir = file.path(rootdir, "inst", "Comparison")
rootdir = system.file(package="ENARSC2015")
# outdir = tempdir()
# outdir = path.expand(outdir)
# stopifnot(file.exists(outdir))


# files = system.file(file.path("NIfTI", niis), package="ENARSC2015")
# roifile = system.file(file.path("NIfTI", "ROI.nii.gz"), package="ENARSC2015")

ext = ".nii.gz"
mods = c("T1", "T2", "FLAIR")
pipes = c("FSL", "ANTsR")

out.dirs = sapply( pipes,  function(x) {
  file.path(rootdir,  x)
})

antsout = out.dirs['ANTsR']
fslout = out.dirs['FSL']

## ----n3_compare, cache=FALSE, results = 'hide'---------------------------
n3.files = lapply(pipes, function(proc){
  n3 = file.path(out.dirs[proc], paste0(mods, "_", proc, "_N3Correct", ext))
  ex = file.exists(n3)
  names(n3) = mods
#   print(ex)
  stopifnot(all(ex))  
  n3
})
names(n3.files) = pipes

n3.imgs = llply(n3.files, function(x){
  im = lapply(x, readNIfTI, reorient=FALSE)
  names(im) = names(x)
  im
}, .progress = "text")

n3.diff = list(T1 = niftiarr(n3.imgs$ANTsR$T1, 
                             n3.imgs$ANTsR$T1 - n3.imgs$FSL$T1),
               T2 = niftiarr(n3.imgs$ANTsR$T2, 
                             n3.imgs$ANTsR$T2 - n3.imgs$FSL$T2),
               FLAIR = niftiarr(n3.imgs$ANTsR$FLAIR, 
                             n3.imgs$ANTsR$FLAIR - n3.imgs$FSL$FLAIR)
               )


## ----ortho_n3, cache=TRUE------------------------------------------------
orthographic(n3.diff$T1, text = "N3 Corr Diff T1 \n FSL vs. ANTsR")
orthographic(n3.diff$T2, text = "N3 Corr Diff T2 \n FSL vs. ANTsR")
orthographic(n3.diff$FLAIR, text = "N3 Corr Diff FLAIR \n FSL vs. ANTsR")

## ----reg_compare, cache=TRUE, results='hide'-----------------------------
reg.files = list(
  FSL = c(FLAIR = file.path(fslout, paste0("FLAIR_FSL_N3Correct_regMNI_FNIRT", ext)), 
          ROI = file.path(fslout, paste0("ROI_threshold_regMNI_FNIRT", ext))),
  ANTsR = c(FLAIR = file.path(antsout, paste0("FLAIR_ANTsR_N3Correct_regMNI_SyN", ext)), 
          ROI = file.path(antsout, paste0("ROI_threshold_regMNI_SyN", ext)))
                            )

reg.imgs = llply(reg.files, function(x){
  
  im = lapply(x, readNIfTI, reorient=FALSE)
  names(im) = names(x)
  im
}, .progress = "text")


## ----plot_reg, cache=TRUE------------------------------------------------

antsROI = reg.imgs$ANTsR$ROI
fslROI = reg.imgs$FSL$ROI

antsROI[ antsROI == 0 ] = NA
fslROI[ fslROI == 0 ] = NA

orthographic(reg.imgs$ANTsR$FLAIR, antsROI,
             col.y=  alpha("red", .25),
             text = "Registered FLAIR \n ANTsR SyN")
orthographic(reg.imgs$FSL$FLAIR, fslROI,
             col.y=  alpha("red", .25),             
             text = "Registered FLAIR \n FSL FNIRT")

reg.diff = list(FLAIR = niftiarr(reg.imgs$ANTsR$FLAIR, 
                             reg.imgs$ANTsR$FLAIR - reg.imgs$FSL$FLAIR),
                ROI = niftiarr(reg.imgs$ANTsR$ROI, 
                             reg.imgs$ANTsR$ROI - reg.imgs$FSL$ROI)
               )

## ----atlas_breakdown, cache=TRUE, results='hide'-------------------------

data(Atlas_Labels)
# hoxsubcort

get.ind = function(img, df){
  df = df[order(df$index),]
	ind.list = llply(df$index, function(x){
		return(which(img %in% x))
	}, .progress = "text")
	names(ind.list) = df$Label
	return(ind.list)
}

hoxsubcort.list = get.ind(hoxsubcort.img, hoxsubcort.df)
  
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

antsROI = reg.imgs$ANTsR$ROI
fslROI = reg.imgs$FSL$ROI

antsArea = tab.area(antsROI, hoxsubcort.list)
fslArea = tab.area(fslROI, hoxsubcort.list)


## ----, results='asis', echo=FALSE----------------------------------------
cat("### ANTsR ROI Breakdown\n")
xtab = xtable(antsArea)
print.xtable( xtab, type="html")

cat("### FSL ROI Breakdown\n")
xtab = xtable(fslArea)
print.xtable( xtab, type="html")


## ----plot_atlas, cache = TRUE--------------------------------------------
xyz = ceiling(cog(antsROI))
res = hoxsubcort.img
res[ antsROI == 0 ] = NA
res[ res %in% -99] = 0
res = cal_img(res)
orthographic(reg.imgs$ANTsR$FLAIR, res, 
             xyz=xyz, col.y= alpha(rainbow(n = res@cal_max), 0.25))


xyz = ceiling(cog(fslROI))
res = hoxsubcort.img
res[ fslROI == 0 ] = NA
res[ res %in% -99] = 0
res = cal_img(res)
orthographic(reg.imgs$FSL$FLAIR, res, 
             xyz=xyz, col.y= alpha(rainbow(n = res@cal_max), 0.25))

## ----clean-up, include=FALSE---------------------------------------------
# R compiles all vignettes in the same session, which can be bad
rm(list = ls(all = TRUE))

