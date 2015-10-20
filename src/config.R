#-----------------Alignment--------------------------------------------#
usearch <- ""
if (Sys.info()['sysname'] == "Darwin") { # OS-X
	if(system(paste("if [ -f ", getwd(), "/bin/usearch* ] ; then echo \"yes\" ; else echo \"no\" ; fi", sep=''), intern=TRUE) == "yes") {
		usearch <- "bin/./usearch8*"
	} else {
		stop("USEARCH does not exist in <metamp home>/bin/>.")
	}
}
if (Sys.info()['sysname'] == "Linux") { # Linux
  if(system("if [ -f ~/Projects/metamp/bin/usearch* ] ; then echo \"yes\" ; else echo \"no\" ; fi") == "yes") {
		usearch <- "bin/./usearch8*"
  } else {
		stop("USEARCH does not exist in <metamp home>/bin/>.")
  }
}
#--------------Miscellaneous------------------------------------------#
# Path to installed BLASTParser library:
R_LIBS <- "R_Lib"
# Output file with final clusters:
clust_filename <- "clusters.clstr"
# Output OTU table:
otu_table_filename <- "otu_table.txt"
# Output file with coordinates:
coord_filename <- "coordinates.crd"
# Chimeric reference database:
chime_ref <- "data/gold/gold.fa"
# A directory that contains temporary files:
tmp_dir <- paste(analysis_dir, "/tmp", sep='')
# Keep or not temporary files:
keep_tmp_files <- T