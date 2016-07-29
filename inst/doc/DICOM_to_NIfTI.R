## ----knit-setup, echo=FALSE, results='hide', eval=TRUE-------------------
library(knitr)
opts_chunk$set(echo=TRUE, prompt=FALSE, message=TRUE, warning=TRUE, comment="", cache=TRUE)
# cache.path <- path.expand("~/Dropbox/Packages/ENAR_SC_2015/vignettes/DCM2NII_cache/")
# opts_chunk$set(cache.path=cache.path)

## ----, eval=FALSE--------------------------------------------------------
## install_github('muschellij2/ENARSC2015')
## # install.packages(c("oro.dicom", "oro.nifti"))

## ------------------------------------------------------------------------
dcm2nii <- function(path, ## path to DICOM files
                    outfile = NULL, ## output filename
                    retimg = TRUE, # write file out to disk?
                    ... # additional options to \code{\link{readDICOM}}
                    ){
  dicom <- readDICOM(path, flipud = FALSE, ...)
  nifti <- dicom2nifti(dicom)
  write = !is.null(outfile)
  if (retimg){
    if (is.null(outfile)) {
      outfile = tempfile()
      outfile = nii.stub(outfile)      
    }
  } else {
    stopifnot(!is.null(outfile))
  }
  
  if (write){
    writeNIfTI(nifti, filename = outfile, 
               verbose = TRUE, gzipped = TRUE)
  } 
  return(nifti)
}

## ----readin,  warning=FALSE, cache=TRUE----------------------------------
library(oro.dicom)
library(oro.nifti)

# output directory for nifti files
outdir = tempdir()
outdir = path.expand(outdir)

if (!file.exists(outdir)){
  dir.create(outdir, showWarnings =  FALSE)
}

# Data taken from http://www.osirix-viewer.com/datasets/
mods = c("T1", "T2", "FLAIR", "ROI")
hdrs = imgs = vector(mode="list", length= length(mods))
names(imgs) = names(hdrs) = mods
imod = 1
for (imod in seq_along(mods)){
  mod = mods[imod]
  dicom_path = system.file(file.path("BRAINIX/DICOM", mod), package="ENARSC2015")
  fname = file.path(outdir, mod)
  headers <- readDICOM(dicom_path, pixelData=FALSE)$hdr
  hdrs[[imod]] = dicomTable(headers)
  img = dcm2nii(dicom_path, retimg = TRUE, outfile = fname)
  imgs[[imod]] = img
}

## ----headers-------------------------------------------------------------
class(hdrs[["T1"]])
rn = rownames(hdrs[["T1"]])
head(rn, 2)
rownames(hdrs[["T1"]]) = basename(rn)
hdrs[["T1"]][1:5, 1:5]

## ----header2-------------------------------------------------------------
colnames(hdrs[["T1"]])[1:5]

## ----cn------------------------------------------------------------------
head(gsub(".*-(.*)", "\\1", colnames(hdrs[["T1"]])), 10)

## ----nif_class-----------------------------------------------------------
sapply(imgs, class)
imgs[["FLAIR"]]

## ----img_info------------------------------------------------------------
head(slotNames(imgs[["FLAIR"]]))
imgs[["FLAIR"]]@dim
dim(imgs[["FLAIR"]])
imgs[["FLAIR"]]@pixdim

## ----img_op--------------------------------------------------------------
i = imgs[["FLAIR"]]
mean(i > 0)
i[i < 0 ]= NA

## ----ex_readin-----------------------------------------------------------
fname = file.path(outdir, "FLAIR")
img = readNIfTI(fname)

## ----ortho---------------------------------------------------------------
orthographic(img, text="FLAIR image")
image(img)
image(img, z = 10, plot.type="single")

## ----makenii-------------------------------------------------------------
arr = img > 300
class(arr)
nim = nifti(arr, dim = dim(img), pixdim=pixdim(img))
orthographic(nim)

## ----over, dev='png'-----------------------------------------------------
library(scales)
nim[nim == 0] = NA
nim[nim == TRUE] = 1
orthographic(img, nim, col.y=alpha("red", 0.25))

