#' @title Copy Data for ENAR SC 2015
#'
#' @description Copies files from system.file for you
#' @param outdir Output directory for data
#' @param get_mipav Should it open URL for mipav?
#' @export
copy_data = function(outdir, get_mipav=TRUE){
  niis = list.files(path = system.file(package="ENARSC2015"),
                    pattern="*.nii.gz", recursive=FALSE, full.names=TRUE)
  niis = basename(niis)
  folders = c("BRAINIX/NIfTI", "Kirby21")
#   folders = system.file(folders, package="ENARSC2015")
  files = c(niis, folders)
  if (!file.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }
  for (ifile in seq_along(files)){
    file = files[ifile]
    file.copy(system.file(files[ifile], package="ENARSC2015"), 
              to = outdir, recursive=TRUE)
  }
  if (get_mipav){
    browseURL("http://mipav.cit.nih.gov/download.php")
  }
}