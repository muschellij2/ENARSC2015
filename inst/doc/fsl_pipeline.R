
## ----knit-setup, echo=FALSE, results='hide', eval=TRUE-------------------
rm(list=ls())
library(knitr)
enar_dir = '~/ENARSC2015'
enar_dir = Sys.readlink(enar_dir)
opts_chunk$set(echo=TRUE, prompt=FALSE, message=TRUE, warning=TRUE, comment="")
cache.path <- path.expand(file.path(enar_dir, "/Presentation/fsl_cache/"))
opts_chunk$set(cache.path=cache.path)


## ----, eval=FALSE--------------------------------------------------------
## library(devtools)
## install_github("muschellij2/fslr")


## ------------------------------------------------------------------------
library(fslr)
library(ENARSC2015)
# outdir = "~/Desktop/results"
rootdir = '~/ENARSC2015'
rootdir = Sys.readlink(rootdir)
outdir = file.path(rootdir, "inst", "FSL")
outdir = path.expand(outdir)
stopifnot(file.exists(outdir))

# if (!file.exists(outdir)){
#   dir.create(outdir, showWarnings =  FALSE)
# }

## If you are running interactively, you need to specify fsl.path
options(fsl.path="/usr/local/fsl")

# copy_data(outdir, get_mipav = FALSE)
mods = c("T1", "T2", "FLAIR")
niis = paste0(mods, ".nii.gz")

files = system.file(file.path("NIfTI", niis), package="ENARSC2015")
roifile = system.file(file.path("NIfTI", "ROI.nii.gz"), package="ENARSC2015")

file.copy(roifile, outdir)
outroifile = file.path(outdir, "ROI")
# files = file.path(outdir, "NIfTI", niis)
# roifile = file.path(outdir, "NIfTI", "ROI.nii.gz")

#### remove .nii.gz extension and use basename
stub = nii.stub(files, bn=TRUE)

bias_files = file.path(outdir, paste0(mods, "_FSL_N3Correct"))
names(bias_files) = names(files) = mods
files
bias_files


## ----run_fast, eval= TRUE, cache=TRUE------------------------------------
for (ifile in seq_along(files)){
  bias_file = bias_files[ifile]
  file = files[ifile]
  ext = get.imgext()
  bfile = paste0(bias_file, ext)
  if (!file.exists(bfile)){
    fast(file, opts = "-B --nopve -v", outfile=bias_file)
  
    ### remove extra files from fast
    seg_file = paste0(bias_file, "_seg", ext)
    file.remove(seg_file)

    output = paste0(bias_file, "_restore", ext)
    file.rename(output, paste0(bias_file, ".nii.gz"))
  }
}


## ----bet, eval=TRUE, cache=TRUE------------------------------------------
bet_files = paste0(bias_files, "_BET")
names(bet_files) = mods
ext = get.imgext()
bfiles = paste0(bet_files, ext)
if (!all(file.exists(bfiles))){
  fslbet(infile = bias_files["FLAIR"], outfile = bet_files["FLAIR"], 
         opts = "-v", reorient=TRUE)
  fslbet(infile = bias_files["T1"], outfile = bet_files["T1"], 
         opts = "-v", reorient=TRUE)
  fslbet(infile = bias_files["T2"], outfile = bet_files["T2"], 
         opts = "-v", reorient=TRUE)
}


## ----bet_flirt, eval=TRUE, cache=TRUE------------------------------------
reffile = bet_files["FLAIR"]
infile= bet_files['T1']
ofile = paste0(bet_files['T1'], "_regFLAIR")
omat = paste0(ofile, ".mat")
regt1 = flirt(infile = infile, reffile= reffile, dof = 6, 
      outfile = ofile, omat = omat, opts="-v", retimg=FALSE, 
      reorient=TRUE)

infile= bet_files['T2']
ofile = paste0(bet_files['T2'], "_regFLAIR")
omat = paste0(ofile, ".mat")
regt2 = flirt(infile = infile, reffile= reffile, dof = 6, 
      outfile = ofile, omat = omat, opts="-v", 
      retimg=TRUE, reorient=FALSE)



## ----flirt, eval=TRUE, cache=TRUE----------------------------------------
reffile = bias_files["FLAIR"]
infile= bias_files['T1']
ofile = paste0(bias_files['T1'], "_regFLAIR")
omat = paste0(ofile, ".mat")
regt1 = flirt(infile = infile, reffile= reffile, dof = 6, 
      outfile = ofile, omat = omat, opts="-v", retimg=TRUE,
      reorient = FALSE)

infile= bias_files['T2']
ofile = paste0(bias_files['T2'], "_regFLAIR")
omat = paste0(ofile, ".mat")
regt2 = flirt(infile = infile, reffile= reffile, dof = 6, 
      outfile = ofile, omat = omat, opts="-v", retimg=TRUE, 
      reorient = FALSE)



## ----plotortho, cache=TRUE-----------------------------------------------
bet_files = paste0(bias_files, "_BET")
names(bet_files) = mods
flair = readNIfTI(bet_files['FLAIR'], reorient=FALSE)
roi = readNIfTI(roifile, reorient=FALSE)
mask = flair > 0
regt1 = readNIfTI(paste0(bet_files['T1'], "_regFLAIR"), reorient=FALSE)
regt2 = readNIfTI(paste0(bet_files['T2'], "_regFLAIR"), reorient=FALSE)
orthographic(flair)
orthographic(regt1)
orthographic(regt2)

diff.t1 = flair
diff.t1@.Data = flair - regt1
diff.t1[mask == 0] = min(diff.t1[mask==1])
### to set cal_min and @cal_max for plotting
diff.t1 = cal_img(diff.t1)
orthographic(diff.t1)

diff.t2 = flair
diff.t2@.Data = flair - regt2
diff.t2[mask == 0] = min(diff.t2[mask==1])
### to set cal_min and @cal_max for plotting
diff.t2 = cal_img(diff.t2)
orthographic(diff.t2)

rat.t1 = flair
rat.t1@.Data = flair / regt1
rat.t1[mask == 0] = min(rat.t1[mask==1])
rat.t1[is.infinite(rat.t1)]= 0
rat.t1[rat.t1 > 1000]= 1000
rat.t1 = log(rat.t1 + 1)
### to set cal_min and @cal_max for plotting
rat.t1 = cal_img(rat.t1)
orthographic(rat.t1)

rat.t2 = flair
rat.t2@.Data = flair / regt2
rat.t2[mask == 0] = min(rat.t2[mask==1])
rat.t2[is.infinite(rat.t2)]= 0
rat.t2[rat.t2 > 1000]= 1000
rat.t2 = log(rat.t2 + 1)
### to set cal_min and @cal_max for plotting
rat.t2 = cal_img(rat.t2)
orthographic(rat.t2)



## ----df, cache=TRUE------------------------------------------------------
df = data.frame(FLAIR=c(flair), T1= c(regt1), T2 = c(regt2), ROI=c(roi))
df = df[ c(mask) == 1, ]
samp = df[ sample(nrow(df), size=1e5), ]
plot(FLAIR ~ T1, data = samp, pch='.')
smoothScatter(samp$T1, samp$FLAIR)

plot(FLAIR ~ T2, data = samp, pch='.')
smoothScatter(samp$T2, samp$FLAIR)

plot(T1 ~ T2, data = samp, pch='.')
smoothScatter(samp$T2, samp$T1)


## ----flirt_temp, cache=TRUE----------------------------------------------
atlas = file.path(fsldir(), "data", "standard", "MNI152_T1_1mm.nii.gz")
t1 = paste0(bias_files['T1'], "_regFLAIR")
reg.t1 = paste0(t1, '_regMNI')
omat = paste0(t1, ".mat")
res = flirt(infile = t1, reffile = atlas, 
      omat = omat,
      outfile = reg.t1,
      dof=12, opts= "-v"
      )


## ----fnirt_temp, cache = TRUE--------------------------------------------
fnirt.t1 = paste0(t1, '_regMNI_FNIRT')
cimg = paste0(fnirt.t1, "_coeff")

opts = '--verbose'
opts = paste(opts, sprintf(" --cout=%s --subsamp=8,8,4,2", cimg))
res = fnirt(infile = reg.t1, reffile = atlas, 
      outfile = fnirt.t1,
      opts= opts, 
      intern = FALSE
      )



## ----applywarp, cache = TRUE---------------------------------------------
t2 = paste0(bias_files['T2'], "_regFLAIR")
reg.t2 = paste0(t2, '_regMNI')
fnirt.t2 = paste0(t2, '_regMNI_FNIRT')

flirt_apply(infile = t2, 
            reffile = atlas, 
            initmat = omat,
            outfile = reg.t2, 
            opts = "-v")

fsl_applywarp(infile = reg.t2, 
                 reffile = atlas, 
                 warpfile = cimg,                  
                 outfile = fnirt.t2,               
                 intern=FALSE,
                 opts= "-v")


flair = bias_files['FLAIR']
reg.flair = paste0(flair, '_regMNI')
fnirt.flair = paste0(flair, '_regMNI_FNIRT')

flirt_apply(infile = flair, 
            reffile = atlas, 
            initmat = omat,
            outfile = reg.flair, 
            opts = "-v")

fsl_applywarp(infile = reg.flair, 
                 reffile = atlas, 
                 warpfile = cimg,                  
                 outfile = fnirt.flair,               
                 intern=FALSE,
                 opts= "-v")


roi = nii.stub(outroifile)
reg.roi = paste0(roi, '_regMNI')
fnirt.roi = paste0(roi, '_regMNI_FNIRT')

flirt_apply(infile = roi, 
            reffile = atlas, 
            initmat = omat,
            outfile = reg.roi, 
            opts = "-v")

fsl_applywarp(infile = reg.roi, 
                 reffile = atlas, 
                 warpfile = cimg,                  
                 outfile = fnirt.roi,               
                 intern=FALSE,
                 opts= "-v")

fnirt.troi = paste0(roi, '_threshold_regMNI_FNIRT')

fslthresh(fnirt.roi, outfile = fnirt.troi, thresh = 0.5, opts = "-bin")


# fsl_applywarp()
# cimg

