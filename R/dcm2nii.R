#' @title Dicom to NIfTI wrapper
#'
#' @description Converts a directory of DICOMs to NIfTI
#' @param path path to DICOM files
#' @param outfile output filename
#' @param retimg Return Image
#' @param ... additional options to \code{\link{readDICOM}}
#' @export
#' @importFrom oro.dicom readDICOM
#' @examples \dontrun{
#' dicom_path = system.file("FLAIR_DICOM/Brainix", package="ENARSC2015")
#' img = dcm2nii(dicom_path, retimg = TRUE)
#'}
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
  } else{ 
    return(nifti)
  }
}


##