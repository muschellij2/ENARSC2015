
## ----knit-setup, echo=FALSE, results='hide', eval=TRUE-------------------
rm(list=ls())
library(knitr)
enar_dir = '~/Dropbox/Packages/ENAR_SC_2015'
# opts_chunk$set(echo=TRUE, prompt=FALSE, message=TRUE, warning=TRUE, comment="", cache=FALSE)
# library(knitcitations)


## ----readin--------------------------------------------------------------
library(XML)
library(stringr)
library(oro.nifti)

fsldir = Sys.getenv("FSLDIR")
if (fsldir == "") {
	fsldir = '/usr/local/fsl'
}
fsltemp = file.path(fsldir, "data", "standard")
atlas_dir = file.path(fsldir, "data", "atlases")

bmask = readNIfTI(
  file.path(fsltemp, "MNI152_T1_1mm_brain_mask"))
reg.bmask = bmask > 0


## ------------------------------------------------------------------------
#### extracting Harvard-Oxford Cortical Image an labels###########
atlas = "HarvardOxford"
tatlas_dir = file.path(atlas_dir, atlas)
xmlfile = file.path(atlas_dir, paste0(atlas, "-Subcortical.xml"))

xx = xmlParse(xmlfile)
indices = xpathSApply(xx, "/atlas/data/label", xmlGetAttr, "index")
labels = xpathSApply(xx, "/atlas/data/label", xmlValue)
labs = str_trim(labels)


## ------------------------------------------------------------------------
df = data.frame(labs, stringsAsFactors=FALSE)
colnames(df) = c("Label")
df = rbind(rep("Outside Brain Mask", ncol(df)), df)
df = rbind(rep("Uncategorized", ncol(df)), df)
df$index = c(0, -99, as.numeric(indices) + 1)

hoxsubcort.df = df
hoxsubcort.df$Label = gsub("Ventrical", "Ventricle", 
	hoxsubcort.df$Label)
head(hoxsubcort.df)


## ----imgread-------------------------------------------------------------
img = readNIfTI(
	file.path(tatlas_dir, 
		"HarvardOxford-sub-maxprob-thr0-1mm.nii.gz"))
img[ !reg.bmask ] = -99
uimg = sort(unique(c(img)))
all.ind = sort(unique(c(0, df$index)))
stopifnot(all(uimg %in% all.ind))
hoxsubcort.img = img


## ----save, echo=FALSE----------------------------------------------------
outfile = file.path(enar_dir, "data", "Atlas_Labels.rda")
save(hoxsubcort.img, hoxsubcort.df, file=outfile )


## ----biblio, results='asis'----------------------------------------------
bibliography(bib.style = "authortitle") 


