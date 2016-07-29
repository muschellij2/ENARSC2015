vignettes/Pipeline_Comparison.html vignettes/ENAR_Tutorial.html vignettes/DICOM_to_NIfTI.html: \
	vignettes/Pipeline_Comparison.Rmd \
	vignettes/ENAR_Tutorial.Rmd vignettes/DICOM_to_NIfTI.Rmd
	echo "building_vignettes"
	Rscript -e "library(methods); library(devtools); devtools::build_vignettes()"
inst/doc/ANTsR_pipeline.html:
	cd inst/doc/; Rscript -e "library(knitr); knit2hthml('ANTsR_pipeline.Rmd')"
inst/doc/ANTsR_pipeline.R:
	cd inst/doc/; Rscript -e "library(knitr); purl('ANTsR_pipeline.Rmd')"
inst/doc/fsl_pipeline.html:
	cd inst/doc/; Rscript -e "library(knitr); knit2hthml('fsl_pipeline.Rmd')"
inst/doc/fsl_pipeline.R:
	cd inst/doc/; Rscript -e "library(knitr); purl('fsl_pipeline.Rmd')"
inst/doc/Atlas_ReadIn.html:
	cd inst/doc/; Rscript -e "library(knitr); knit2hthml('Atlas_ReadIn.Rmd')"
inst/doc/Atlas_ReadIn.R:
	cd inst/doc/; Rscript -e "library(knitr); purl('Atlas_ReadIn.Rmd')"
