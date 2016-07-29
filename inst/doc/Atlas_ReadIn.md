---
title: "Atlas to R Object"
author: "John Muschelli"
date: "August 25, 2014"
output: html_document
---

# Making Atlas R Object




## Introduction

We would like to use the Harvard-Oxford atlas, included with FSL for seeing where our region interest (ROI) for the brain tumor is located in the brain.  

We must first read in the MNI atlas to get a brain mask.


```r
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
```

The labels for the brain image are located in an XML file which we will read with the XML package (Temple Lang, 2013).



```r
#### extracting Harvard-Oxford Cortical Image an labels###########
atlas = "HarvardOxford"
tatlas_dir = file.path(atlas_dir, atlas)
xmlfile = file.path(atlas_dir, paste0(atlas, "-Subcortical.xml"))

xx = xmlParse(xmlfile)
indices = xpathSApply(xx, "/atlas/data/label", xmlGetAttr, "index")
labels = xpathSApply(xx, "/atlas/data/label", xmlValue)
labs = str_trim(labels)
```

We will create a `data.frame`, add a label for outside the brain mask, and also add a label for uncategorized information (some small number of voxels is not labeled ).


```r
df = data.frame(labs, stringsAsFactors=FALSE)
colnames(df) = c("Label")
df = rbind(rep("Outside Brain Mask", ncol(df)), df)
df = rbind(rep("Uncategorized", ncol(df)), df)
df$index = c(0, -99, as.numeric(indices) + 1)

hoxsubcort.df = df
hoxsubcort.df$Label = gsub("Ventrical", "Ventricle", 
	hoxsubcort.df$Label)
head(hoxsubcort.df)
```

```
##                        Label index
## 1              Uncategorized     0
## 2         Outside Brain Mask   -99
## 3 Left Cerebral White Matter     1
## 4       Left Cerebral Cortex     2
## 5     Left Lateral Ventricle     3
## 6              Left Thalamus     4
```

We see the values for the image and the corresponding brain structure labels in our `data.frame` and can use this to map areas we find for our ROI onto real cortical structure.

## Read in Image
Let's read in the image, call any values outside of the brain mask to be -99 and then check to see that all values in the image are in the indices from the `data.frame`.


```r
img = readNIfTI(
	file.path(tatlas_dir, 
		"HarvardOxford-sub-maxprob-thr0-1mm.nii.gz"))
img[ !reg.bmask ] = -99
uimg = sort(unique(c(img)))
all.ind = sort(unique(c(0, df$index)))
stopifnot(all(uimg %in% all.ind))
hoxsubcort.img = img
```






### References


```r
bibliography(bib.style = "authortitle") 
```

Temple Lang, D. _XML: Tools for parsing and generating XML within
R and S-Plus._ R package version 3.98-1.1. 2013. <URL:
http://CRAN.R-project.org/package=XML>.

Makris N, Goldstein JM, Kennedy D, Hodge SM, Caviness VS, Faraone SV, Tsuang MT, Seidman LJ. Decreased volume of left and total anterior insular lobule in schizophrenia. Schizophr Res. 2006 Apr;83(2-3):155-71

Frazier JA, Chiu S, Breeze JL, Makris N, Lange N, Kennedy DN, Herbert MR, Bent EK, Koneru VK, Dieterich ME, Hodge SM, Rauch SL, Grant PE, Cohen BM, Seidman LJ, Caviness VS, Biederman J. Structural brain magnetic resonance imaging of limbic and thalamic volumes in pediatric bipolar disorder. Am J Psychiatry. 2005 Jul;162(7):1256-65

Desikan RS, Ségonne F, Fischl B, Quinn BT, Dickerson BC, Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS, Killiany RJ. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. Neuroimage. 2006 Jul 1;31(3):968-80.

Goldstein JM, Seidman LJ, Makris N, Ahern T, O'Brien LM, Caviness VS Jr, Kennedy DN, Faraone SV, Tsuang MT. Hypothalamic abnormalities in schizophrenia: sex effects and genetic vulnerability. Biol Psychiatry. 2007 Apr 15;61(8):935-45
