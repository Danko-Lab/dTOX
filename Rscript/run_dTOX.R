library(bigWig)
library(dREG)
library(Rgtsvm)
library(parallel)
library(randomForest);
library(snowfall);

options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE);

## Read arguments from main thread
species            <- args[1];
seq.type           <- toupper(args[2]);
file.bw.plus       <- args[3];
file.bw.minus      <- args[4];
file.fullpath.ret  <- args[5];
file.home.dtox     <- args[6];
gpu.total          <- as.numeric(args[7])
cpu.total          <- as.numeric(args[8])


cat("bigWig plus=", file.bw.plus, "\n")
cat("bigWig minus=", file.bw.minus , "\n")
cat("Seq type=", seq.type, "\n")
cat("Species=", species, "\n")
cat("Result path=", file.fullpath.ret, "\n")
cat("Home path=", file.home.dtox,  "\n")
cat("GPU total=", gpu.total, "\n")


write_boundlike_bed<-function( file.full.motif, file.filter.model, file.bw.plus, file.bw.minus, file.pred.bed, ncores=1 )
{
	library(parallel)
	library(bigWig)
	library(randomForest);
	library(dREG)
	
    source(paste0(file.home.dtox,"/Rscript/pre-reduce.R"));

	tb <- try(read.table(file.full.motif));
	if(class(tb)=="try-error")
	{   
		cat("Error in ", file.full.motif, "\n");
		return("Error");
	}	 

    colnames(tb) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V4.1", "V5.1", "V6.1");
	
    #tb.mat0 <- cbind( get_read_mat_slow(tb[,c(1:3)], file.bw.plus, ncores),  get_read_mat_slow(tb[,c(1:3)], file.bw.minus, ncores ))/ get_total_number(file.bw.plus, file.bw.minus) * 10^9;
    #tb.mat0 <- cbind(tb[,-c(1:4, 6)], tb.mat0);
    tb.mat <- get_read_mat_dreg(tb[,c(1:3)], file.bw.plus, file.bw.minus, ncores=ncores)
    tb.mat <- cbind(tb[,-c(1:4, 6)], tb.mat);
    
	load(file.filter.model);

    tb.pred <- predict(model, newdata=tb.mat);
    tb <- tb[tb.pred>0.5,, drop=T ]
    tmp.file <- "None";
    if(NROW(tb)>0)
    {
        write.bed( tb, file=file.pred.bed, compress=TRUE);
        tmp.file <- file.pred.bed;
    }
    
	rm(tb);
	rm(tb.mat);
	gc(reset=TRUE);
    return(tmp.file); 	
} 

get_ps_count<-function(ps.key)
{
   ret <- read.table(pipe(paste("ps -ef | grep ", ps.key, " | grep -v 'grep' | wc -l ")));
   return(ret$V1[1]);
}

predict_dTOX <- function(file.bw.plus, file.bw.minus, file.filter.model, file.dtox.model, gpu.idx=0, cpu.cores=16 )
{
    source(paste0(file.home.dtox,"/Rscript/pre-reduce.R"));

	path.result <- tempfile( pattern="tmp.pred", tmpdir =".", fileext = "")
	motif.list.table <- tempfile(pattern="tmp.motif.list", tmpdir =".", fileext = ".tab")
 	dir.create(path.result);

cat("motif.list.table=", motif.list.table, "\n");	
cat("path of result=", path.result, "\n");	


	# invoke multiple thread running on GPU cores
	for(gpu in gpu.idx)
	{
	   system(paste("R --vanilla --slave --args dTOX-svm-predict-thread ", gpu, which(gpu==gpu.idx), NROW(gpu.idx), file.dtox.model, file.bw.plus, file.bw.minus, motif.list.table, path.result, " < ", paste0(file.home.dtox, "/Rscript/motif-pred-gpu-thred.R ")), wait=FALSE, intern=FALSE)
	}

    cpu.fun<-function(file.motif)
    {
	    file.full.motif <- paste0(path.motif, "/", file.motif)
        file.pred.bed <- paste0( path.result, "/pred-",file.motif);
	    file.tmp <- write_boundlike_bed(file.full.motif, file.filter.model, file.bw.plus, file.bw.minus, file.pred.bed , ncores=1);
	    return(file.tmp)
    }

    sfInit(parallel = TRUE, cpus = cpu.cores, type = "SOCK" )
    sfExport("file.filter.model", "file.bw.plus", "file.bw.minus", "path.motif", "path.result", "file.home.dtox", "write_boundlike_bed" );
    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();
    
	file.beds <- list.files(pattern = "\\.bed.gz$", path=path.motif);
	L0 <- lapply (1:ceiling(NROW(file.beds)/80), function(i)
	{
        file.round.beds <- file.beds[ (i-1)*80+c(1:80) ];
        file.round.beds <- file.round.beds[!is.na(file.round.beds)];
        file.pred.beds <- unlist(sfLapply(file.round.beds, fun));
show(file.pred.beds);

        ret.code <- 0;
	    for (fs.pred in file.pred.beds) 
	        if( fs.pred != "Error" && fs.pred != "None" )
	        {
	            cat( fs.pred, "\n", file=motif.list.table, append=TRUE);
	            ret.code <- ret.code+1;
	        } 
	    return(ret.code)    
	});

    sfStop();

	cat( "END\n", file=motif.list.table, append=TRUE);

	motif.file.count <- sum(unlist(L0));
cat("Motif file count=", motif.file.count, "\n");

	repeat
	{
	    file.beds <- list.files(pattern = "\\.bound.gz$", path=path.result);
	    if (NROW(file.beds)==motif.file.count)
	        break;
	    
	    if( get_ps_count("dTOX-svm-predict-thread") == 0 )
	       break;

	    Sys.sleep(60);
	}

	Sys.sleep(300);
	
	L0 <- mclapply (list.files( pattern = "\\.bound.gz$", path=path.result), function(file.bound) {
	    tb<-try(read.table(paste0(path.result, "/", file.bound)));
	    if(class(tb)=="try-error")
	       return(NULL)
	    else
	       return(tb)
	 
	}, mc.cores=10);
  
	
	bed.bound <- rbindlist(L0);
	write.bed( as.data.frame(bed.bound), file.fullpath.ret, compress=TRUE, mkindex=TRUE );

	cat("Bound motif count=", NROW(bed.bound), "\n");
	
	unlink( path.result, recursive=TRUE );
	rm(L0);
	rm(bed.bound);
}	

if (species!="hg19")
	stop("Currently only hg19 can be predicted.")

path.motif <- paste0(file.home.dtox, "/TFBS/motif_rtfbsdb2_hg19_ext")
	
if (seq.type=="ATAC-SEQ")
{
	file.filter.model = paste0(file.home.dtox, "/models/atac-reduce2-model-400K-chr1.rdata");
	file.dtox.model   = paste0(file.home.dtox,"/models/kg_atac_model_v1_20181205.rdata");
}

if (seq.type=="DNASE-SEQ")
{
	file.filter.model = paste0(file.home.dtox, "/models/dnase-reduce2-model-400K-chr1.rdata");
	file.dtox.model   = paste0(file.home.dtox, "/models/kg_dnase_model_v1_20181205.rdata");
}

if (seq.type=="PRO-SEQ")
{
	file.filter.model = paste0(file.home.dtox, "/models/proseq-reduce2-model-400K-chr1.rdata");
	file.dtox.model   = paste0(file.home.dtox, "/models/kg_proseq_model_v1_20181205.rdata");
}

ti <- system.time( predict_dTOX (file.bw.plus, file.bw.minus, file.filter.model, file.dtox.model, gpu.idx=c(1:gpu.total)-1, cpu.cores=cpu.total ) ); 
print(ti);


