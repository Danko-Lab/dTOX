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
str.out.prefix     <- args[5];
file.home.dtox     <- args[6];
using.filter       <- as.logical(args[7]);
gpu.total          <- as.numeric(args[8])
cpu.total          <- as.numeric(args[9])

file.fullpath.ret=paste( str.out.prefix, ".dTOX.bound.bed.gz", sep="")
file.infp.ret=paste( str.out.prefix, ".infp.bed.gz", sep="")

cat("bigWig plus=", file.bw.plus, "\n")
cat("bigWig minus=", file.bw.minus , "\n")
cat("Seq type=", seq.type, "\n")
cat("Species=", species, "\n")
cat("Result path=", file.fullpath.ret, "\n")
cat("Home path=", file.home.dtox,  "\n")
cat("Filter=", if(using.filter)"YES" else "NO", "\n" );
cat("GPU total=", gpu.total, "\n")
cat("CPU total=", cpu.total, "\n")



get_info_pos <- function(bigwig_plus, bigwig_minus)
{
    # source("get_info_pos_k562_atacNEW.r");

    infoPos <- get_informative_positions(bw_path = bigwig_plus, bw_minus_path = bigwig_minus, depth= 0, window= 400, step=50, use_OR=TRUE, use_ANDOR=TRUE, debug= TRUE);
    file.temp <- tempfile(fileext=".bed");
    write.table(infoPos, file=file.temp, quote=F ,row.names=F, col.names=F, sep="\t");

    # intersect with
    tb <- read.table(pipe(paste0("bedtools merge -i ", file.temp)));
    unlink(file.temp);
    return(tb);
}

write_boundlike_bed<-function( file.full.motif, file.filter.model, file.bw.plus, file.bw.minus, file.pred.bed, filter.thres=0.5, ncores=1 )
{
    library(parallel)
    library(bigWig)
    library(randomForest);
    library(dREG)
    
    source(paste0(file.home.dtox,"/Rscript/pre-reduce.R"));
  
    if( !file.exists(file.full.motif) )
        return("Error");

    ## generate the bed file including chr, start,end, score in motif bed and flag_inclduing_infp
    tb <- try( read.table( pipe(paste( "zcat ", file.full.motif, " | awk 'BEGIN{OFS=\"\t\"}{ print $0 }' | bedmap --echo --delim '\t' --indicator - ", file.infp.bed, " | sort-bed - " ) ) ) );
    if(class(tb)=="try-error")
    {   
        cat("Error in ", file.full.motif, "\n");
        return("Error");
    }     
    
    if(NROW(tb)==0)
       return("None");
  
    # select the TFBS with infp, remove the 10th column(indicator 0 or 1) 
    tb <- tb[tb$V9==1, -9 ];
    if(using.filter) 
    {
        colnames(tb) <- c("V1", "V2", "V3", "V5", "V4.1", "V5.1", "V6.1", "V7");
        tb.mat <- get_read_mat_dreg(tb[,c(1:3)], file.bw.plus, file.bw.minus, ncores=ncores)
        tb.mat <- cbind(tb[,-c(1:3, 8)], tb.mat);
    
        load(file.filter.model);

        tb.pred <- predict(model, newdata=tb.mat);

        #tb <- tb[tb.pred>filter.thres, , drop=T ]
        tb <- tb[as.numeric(levels(tb.pred))[tb.pred]==1, , drop=T ]
        tmp.file <- "None";
        if(NROW(tb)>0)
        {
            write.bed( data.frame(tb[, c(1:3, 8, 4)], "."), file=file.pred.bed, compress=TRUE);
            tmp.file <- file.pred.bed;
        }
           rm(tb.mat);
       }    
    else
    {
        write.bed( data.frame(tb[, c(1:3, 8, 4)], "."), file=file.pred.bed, compress=TRUE);
        tmp.file <- file.pred.bed;
    }
    
    rm(tb);
    gc(reset=TRUE);
    return(tmp.file);     
} 

get_ps_count<-function(ps.key)
{
   ret <- read.table(pipe(paste("ps -ef | grep ", ps.key, " | grep -v 'grep' | wc -l ")));
   return(ret$V1[1]);
}

predict_dTOX <- function(file.bw.plus, file.bw.minus, file.dtox.model, file.filter.model, filter.threshold, using.filter=TRUE, gpu.idx=0, cpu.cores=16 )
{
    source(paste0(file.home.dtox,"/Rscript/pre-reduce.R"));

    file.infp.bed <- tempfile( pattern="tmp.infp.", tmpdir =".", fileext = ".bed")
    infp.tb <- get_info_pos( file.bw.plus, file.bw.minus );
    infp.tb[,3] <- infp.tb[,3] + 49;
    write.bed( infp.tb, file.infp.bed);
    write.bed( infp.tb, file.infp.ret, compress=TRUE);

    path.result <- tempfile( pattern="tmp.pred.", tmpdir =".", fileext = "")
    motif.list.table <- tempfile(pattern="tmp.motif.", tmpdir =".", fileext = ".tab")
     dir.create(path.result);

cat("motif.list.table=", motif.list.table, "\n");    
cat("infp file=", file.infp.bed, "\n");    
cat("path of result=", path.result, "\n");    
cat("File Infp=", file.infp.ret, "\n");    

    # invoke multiple thread running on GPU cores
    for(gpu in gpu.idx)
    {
       system(paste("R --vanilla --slave --args dTOX-svm-predict-thread ", gpu, which(gpu==gpu.idx), NROW(gpu.idx), file.dtox.model, file.bw.plus, file.bw.minus, motif.list.table, path.result, " < ", paste0(file.home.dtox, "/Rscript/motif-pred-gpu-thred.R ")), wait=FALSE, intern=FALSE)
    }

    cpu.fun<-function(file.motif)
    {
        file.full.motif <- paste0(path.motif, "/", file.motif)
        file.pred.bed <- paste0( path.result, "/pred-",file.motif);
        file.tmp <- write_boundlike_bed(file.full.motif, file.filter.model, file.bw.plus, file.bw.minus, file.pred.bed, filter.thres=filter.threshold, ncores=1);
        return(file.tmp)
    }

    sfInit(parallel = TRUE, cpus = cpu.cores, type = "SOCK" )
    sfExport("file.filter.model", "file.bw.plus", "file.bw.minus", "path.motif", "path.result", "file.home.dtox", "write_boundlike_bed", "file.infp.bed", "using.filter", "filter.threshold" );
    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();
    
    motif.TFBS.count <<- 0;
    file.beds <- list.files(pattern = "\\.bed.gz$", path=path.motif);
    L0 <- lapply (1:ceiling(NROW(file.beds)/48), function(i)
    {
        file.round.beds <- file.beds[ (i-1)*48+c(1:48) ];
        file.round.beds <- file.round.beds[!is.na(file.round.beds)];
        file.pred.beds <- unlist(sfLapply(file.round.beds, fun));

        ret.code <- 0;
        for (fs.pred in file.pred.beds) 
            if( fs.pred != "Error" && fs.pred != "None" )
            {
                motif.TFBS.count <<- motif.TFBS.count + NROW(read.table(fs.pred));
                cat( fs.pred, "\n", file=motif.list.table, append=TRUE);
                ret.code <- ret.code+1;
            } 
        
        gc(reset=TRUE);
        return(ret.code)    
    });

    sfStop();

    cat( "END\n", file=motif.list.table, append=TRUE);

    motif.file.count <- sum(unlist(L0));

    repeat
    {
        file.beds <- list.files(pattern = "\\.bound.gz$", path=path.result);
        if (NROW(file.beds)==motif.file.count)
            break;
        
        if( get_ps_count("dTOX-svm-predict-thread") == 0 )
           break;

        Sys.sleep(60);
    }

    Sys.sleep(120);
    
    L0 <- mclapply (list.files( pattern = "\\.bound.gz$", path=path.result), function(file.bound) {
        tb<-try(read.table(paste0(path.result, "/", file.bound)));
        if(class(tb)=="try-error")
           return(NULL)
        else
           return(tb)
     
    }, mc.cores=round((cpu.cores+1)/2));
  
    
    bed.bound <- rbindlist(L0);
    write.bed( as.data.frame(bed.bound), file.fullpath.ret, compress=TRUE, mkindex=TRUE );

cat("Motif file count=", motif.file.count, "\n");
cat("Motif TFBS count=", motif.TFBS.count, "\n");
cat("Bound motif count=", NROW(bed.bound), "\n");
    
    unlink( motif.list.table );
    unlink( file.infp.bed );
    unlink( path.result, recursive=TRUE );
    rm(L0);
    rm(bed.bound);
}    

if (species!="hg19" && species!="mm10")
   stop("Currently only hg19 and mm10 can be predicted.")

if (species=="hg19")
   path.motif <- paste0(file.home.dtox, "/TFBS/motif_rtfbsdb_hg19_ext2")

if (species=="mm10")
   path.motif <- paste0(file.home.dtox, "/TFBS/motif_rtfbsdb_mm10_ext2")
    
if (seq.type=="ATAC-SEQ")
{
    file.dtox.model   = paste0(file.home.dtox,"/models/kg_atac_model.rdata");
    file.filter.model = paste0(file.home.dtox, "/models/atac-reduce2-model-400K-chr1.rdata");
    filter.threshold = 0.5;
}

if (seq.type=="DNASE-1-SEQ")
{
    #file.filter.model = paste0(file.home.dtox, "/models/dnase-reduce2-model-400K-chr1.rdata");
    #file.dtox.model   = paste0(file.home.dtox, "/models/kg_dnase_model.rdata");
    file.dtox.model   = paste0(file.home.dtox, "/models/kg_dnase_model_6.7.19.rdata");
    file.filter.model = paste0(file.home.dtox, "/models/dnase-rf-model-chr1.rdata");
    filter.threshold = 0.01;
}

if (seq.type=="PRO-SEQ")
{
    file.dtox.model   = paste0(file.home.dtox, "/models/kg_proseq_model.rdata");
    file.filter.model = paste0(file.home.dtox, "/models/proseq-reduce2-model-400K-chr1.rdata");
    filter.threshold = 0.5;
}

ti <- system.time( predict_dTOX (file.bw.plus, file.bw.minus,  file.dtox.model, file.filter.model, filter.threshold, using.filter=using.filter, gpu.idx=c(1:gpu.total)-1, cpu.cores=cpu.total ) ); 
print(ti);
