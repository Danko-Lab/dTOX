install.packages("rphast", repos = "http://cran.us.r-project.org")
install.packages("rtfbs", repos = "http://cran.us.r-project.org")
install.packages("cluster", repos = "http://cran.us.r-project.org")
install.packages("apcluster", repos = "http://cran.us.r-project.org")
install.packages("vioplot", repos = "http://cran.us.r-project.org")
install.packages("tools", repos = "http://cran.us.r-project.org")
install.packages("latticeExtra", repos = "http://cran.us.r-project.org")
install.packages("lattice", repos = "http://cran.us.r-project.org")
install.packages("Rcurl", repos = "http://cran.us.r-project.org")

library(devtools)
install_github("/Danko-Lab/rtfbs_db/rtfbsdb");

if ( any(Sys.which(c('bedtools', 'bedops', 'twoBitToFa', 'awk', 'sort-bed'))=="") )
{
	Stop("Need to confirm all required tools are installed (awk, bedops. bedtools, twoBitToFa)");
}



