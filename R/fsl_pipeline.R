#' @title FSL Analysis Pipeline
#'
#' @description Run my analysis pipeline
#' @param t1 filename of T1 Image
#' @param t2 filename of T2 Image
#' @param flair filename of FLAIR Image
#' @param pd filename of PD Image
#' @param bias_correction Should Bias correction be done?
#' @param ... pass to other stuff
#' @export
#' @import fslr
#' @import plyr
#' @import oro.nifti
#' @return Stuff
fsl_pipeline = function(t1, # filename of T1 Image
	t2, # filename of T2 Image
	flair, # filename of FLAIR Image
	pd, # filename of PD Image
  bias_correction = TRUE,
  ...
	){
  
  outdir = "~/Desktop/results"
  outdir = path.expand(outdir)
  
  options(fsl.path="/usr/local/fsl")
  
  # Data taken from http://www.osirix-viewer.com/datasets/
#   img.names = c("T1", "T2", "FLAIR", "PD")
#   fnames = paste0(img.names, "norm", ".nii.gz")
#   
#   imgs = system.file(paste0("Baseline_MRI/", fnames), package="ENARSC2015")
#   ### read in Data
#   t1 = imgs[1]
#   t2 = imgs[2]
#   flair = imgs[3]
#   pd = imgs[4]
  
#   outfile = file.path(outdir, paste0(nii.stub(imgs, bn=TRUE), "_N3Correct"))
  dicom_path = system.file("FLAIR_DICOM/Brainix", package="ENARSC2015")
  outfile = file.path(outdir, "FLAIR_NIfTI")
  img = dcm2nii(dicom_path, retimg = FALSE, outfile = outfile)
  
  bias_file = paste0(outfile, "_FSL_N3Correct")
  
  fast(outfile, opts = "-B --nopve -v", outfile=bias_file)

  ext = get.imgext()
  seg_file = paste0(bias_file, "_seg", ext)
  file.remove(seg_file)

  bias_file = paste0(bias_file, "_restore")

}



