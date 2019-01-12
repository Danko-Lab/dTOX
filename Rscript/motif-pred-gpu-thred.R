# Prototype:
# R --vanilla --slave --args dTOX-svm-predict-thread $GPU.idx $GPU.total $model $bw.plus $bw.minus $motif.table $path.tmp < motif-pred-gpu-thred.R
#
# example:
# R --vanilla --slave --args dTOX-svm-predict-thread 0 1 2 /home/zw355/src/dTOX/models/kg_dnase_model_v1_20181205.rdata /fs/cbsudanko/storage/projects/TFBindingSVR/dTOX/hela/dnase/wgEncodeUwDnaseHelaAln_plus.bw /fs/cbsudanko/storage/projects/TFBindingSVR/dTOX/hela/dnase/wgEncodeUwDnaseHelaAln_minus.bw ./tmp.motif.listb375b331e4e.tab ./tmp.predb371b9fedec

library(bigWig)
library(dREG)
library(Rgtsvm)

options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE);

## Read arguments from main thread
ps.key     <- args[1]
gpu.devno  <- as.numeric(args[2])
gpu.idx    <- as.numeric(args[3])
gpu.total  <- as.numeric(args[4])
file.dtox.model  <- args[5]
file.bw.plus     <- args[6];
file.bw.minus     <- args[7];
file.motif.table  <- args[8]
file.path.result  <- args[9]

write.bed<-function ( df.bed, file.bed, compress=FALSE, mkindex=FALSE )
{
	options("scipen"=100, "digits"=4);
	temp <- tempfile(fileext=".bed");
	write.table( df.bed, file=temp, quote=F, row.names=F, col.names=F, sep="\t");
	if(compress) 
	   system(paste0("sort-bed ", temp,  " | bgzip  > ",  file.bed ))
	else
	   system(paste0("sort-bed ", temp,  " > ",  file.bed ))

    if (mkindex)
        system(paste("tabix -p bed",file.bed));
    
	invisible(unlink(temp));
}

get_ps_count<-function(ps.key)
{
   ret <- read.table(pipe(paste("ps -ef | grep ", ps.key, " | grep -v 'grep' | wc -l ")));
   return(ret$V1[1]);
}

PredictOnRtfbsdb <- function( def.ncores=4 )
{
	Motif_predict <- function(file.motif, ncores)
	{
		test.bed <- file.motif;
		file.ret.bed <- paste0(  substring(file.motif, 1, nchar(file.motif)-7), ".bound.gz" );

		if(file.exists( file.ret.bed ))
			return( -1 );

	    tbr <- read.table(test.bed, header=F)

  	    mat = read_genomic_data(gdm, as.data.frame(tbr[,1:3]), file.bw.plus, file.bw.minus, ncores=ncores)
	    mat = cbind( mat, tbr[,4])

	    gt.predict <- predict.run( gt.model, newdata=mat, decision.values=TRUE)
	    scores  <- attr( gt.predict, "decision.values" )

		tb_pred <- cbind(tbr[,c(1:6)], p_score=round(scores,4), p_status=gt.predict);
        tb_pred <- tb_pred[ tb_pred$p_status==1,,drop=F ]
        if(NROW(tb_pred)>0)
		   write.bed( tb_pred, file=file.ret.bed, compress=TRUE);

		gc(reset=TRUE);
	    return( 0 );
	}


	load(file.dtox.model)
	gt.model <- predict.load(gt.model);
    
    task.last.row <- 0;
    repeat
    {
        if(!file.exists(file.motif.table) )
            {Sys.sleep(60);next;}
        
        motif.files <- as.vector(read.table(file.motif.table)$V1);
        if(NROW(motif.files)<1)
            {Sys.sleep(60);next;}


        last.motif.item <- motif.files[NROW(motif.files)];
        motif.files.ends <- rep( if(last.motif.item=="END") "END" else "NA", gpu.total - NROW(motif.files)%%gpu.total );
		motif.file.matrix <- matrix(c(motif.files, motif.files.ends), ncol=gpu.total, byrow=TRUE);
		
		task.motif.file <- motif.file.matrix[, gpu.idx];
		if (NROW(task.motif.file) <= task.last.row )
	    {
	       if(last.motif.item=="END" )
	          break
	        else
	          {Sys.sleep(60);next;}
		}
		
		if( task.motif.file[ task.last.row + 1 ] == "END" )
		    break;
		
		if( task.motif.file[ task.last.row + 1 ] == "NA" )
	        {Sys.sleep(60);next;}
 		
		file.motif <- task.motif.file[ task.last.row + 1 ];
		ncores <- def.ncores;
        if( get_ps_count("RSOCKnode.R") < 8 )		
           ncores <- def.ncores*2;   
		Motif_predict(file.motif, ncores);
		
		task.last.row <- task.last.row + 1
    }
    
	predict.unload(gt.model);
	return();
}

selectGPUdevice( gpu.devno );
PredictOnRtfbsdb( def.ncores=4 );
