## NOTE: In process of making installation easier, as follows:
##       The bigWig package provides all required Kent source dependencies,
##       and should be appearing in CRAN soon.
#R_LIBS := ~/bin/dREG
#export R_LIBS

.PHONY: dreg

R_dependencies:
	@echo "Installing R dependencies" # to:" ${R_LIBS}
	#mkdir -p ${R_LIBS}
	R --no-save < rDeps.Rtfbsdb.R

dreg:
	@echo "Installing dREG" # to:" ${R_LIBS}
	R --no-save < rDeps.dREG.R
	@echo "Installing dREG" # to:" ${R_LIBS}
	#mkdir -p ${R_LIBS}
	#make topLibs -C kent/src
	#R CMD INSTALL bigWig --clean
	#make -C dREG/src/lib
	R CMD INSTALL dREG --clean

Rgtsvm:
	@echo "Installing Rgtsvm" # to:" ${R_LIBS}
	R --no-save < rDeps.Rgtsvm.R
	R CMD INSTALL Rgtsvm --clean

dTOX:
        git clone https://github.com/Danko-Lab/dTOX
        wget -r -nH --cut-dirs=3 ftp://cbsuftp.tc.cornell.edu/danko/hub/dTOX/*

download:
      wget -r -nH --cut-dirs=3 ftp://cbsuftp.tc.cornell.edu/danko/hub/dTOX/*
