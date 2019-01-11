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
        system(paste("tabix -p bed",file.bed) );
	   
	invisible(unlink(temp));
}

get_reads<-function(infoReads, file.bw.plus, file.bw.minus, range=50)
{
    infoReads[,2] <- infoReads[,2] - range;
    infoReads[,3] <- infoReads[,3] + range;
    infoReads[infoReads[,2] < 0, 2] <- 0;
    
	bw.plus <- load.bigWig( file.bw.plus );
	bw.minus <- load.bigWig( file.bw.minus );
	readPlus  <- bed.region.bpQuery.bigWig(bw.plus, infoReads[, c(1:3)], op = "sum", abs.value = TRUE);
	readMinus <- bed.region.bpQuery.bigWig(bw.minus, infoReads[, c(1:3)], op = "sum", abs.value = TRUE);
	unload.bigWig( bw.plus );
	unload.bigWig( bw.minus );
	reads= readPlus + readMinus
	return(reads)
}

get_read_onebigwig<-function(infoReads, file.bw, range=50)
{
    inf_c <- floor((infoReads[,2]+infoReads[,3])/2);
    infoReads[,2] <- inf_c - range;
    infoReads[,3] <- inf_c + range;
    infoReads[infoReads[,2] < 0, 2] <- 0;
    
	bw <- load.bigWig( file.bw );
	reads  <- bed.region.bpQuery.bigWig(bw, infoReads[, c(1:3)], op = "sum", abs.value = TRUE);
	unload.bigWig( bw );
	return(reads)
}

get_read_mat_slow <- function(bed3, file.bw,  ncores=10)
{
    block= 100000;
    
	df <- do.call("rbind", mclapply(1:ceiling(NROW(bed3)/block), function(i){
		i.start <- (i-1)* block + 1
		i.stop <- (i-1)* block + block
		if(i.stop > NROW(bed3)) i.stop <- NROW(bed3); 

		read50  <- get_read_onebigwig( bed3[i.start:i.stop, c(1:3)], file.bw, 50);
		read250 <- get_read_onebigwig( bed3[i.start:i.stop, c(1:3)], file.bw, 250);
		read1k  <- get_read_onebigwig( bed3[i.start:i.stop, c(1:3)], file.bw, 1000);
		read5k  <- get_read_onebigwig( bed3[i.start:i.stop, c(1:3)], file.bw, 5000);
	
        df0 <- as.data.frame(cbind(read50, read250, read1k, read5k));
		return(df0);
	}, mc.cores=ncores));
    
    return(as.data.frame(df));
}

get_read_mat_dreg <- function(bed3, file.bw.plus, file.bw.minus, ncores=10)
{
    gdm <- genomic_data_model(window_sizes = c(50, 250, 1000, 5000), half_nWindows= c(1, 1, 1, 1 ) );
    mat <- read_genomic_data( gdm, bed3, file.bw.plus, file.bw.minus, as_matrix = TRUE, scale.method ="linear", batch_size = 100000, ncores = ncores )* 10^9
    mat0 <- cbind( mat[,1]+mat[, 2], mat[, 3]+mat[, 4], mat[, 5]+mat[, 6], mat[, 7]+mat[, 8],
                   mat[,9]+mat[,10], mat[,11]+mat[,12], mat[,13]+mat[,14], mat[,15]+mat[,16] );
                   
    colnames(mat0) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8");
    return(mat0);
}

get_read_mat_delt <- function(bed3, file.bw, ncores=10)
{
    r <- get_read_mat(bed3, file.bw,  ncores); 	
    r$diff250 <- r[,2] - r[,1];
    r$diff1K  <- r[,3] - r[,2];
    r$diff5K  <- r[,4] - r[,3];
    return(r);
}

merge_bed <- function(bed3, distance=50)
{
    tmp.bed <- tempfile(fileext=".bed");
	write.bed( bed3, tmp.bed );
	
    tbr <- read.table(pipe(paste("bedtools merge -d", distance, " -i",  tmp.bed)))
    unlink(tmp.bed);
    return(tbr);
}

substract_bed <- function(bedA, bedB)
{
    tmpA.bed <- tempfile(fileext=".bed");
	write.bed( bedA, tmpA.bed );
    tmpB.bed <- tempfile(fileext=".bed");
	write.bed( bedB[,c(1:3)], tmpB.bed );
	
    tbr <- read.table(pipe(paste("bedtools subtract  -a", tmpA.bed, " -b",  tmpB.bed)))
    unlink(tmpA.bed);
    unlink(tmpB.bed);
    return(tbr);
}

get_total_number <- function(file.bw.plus, file.bw.minus)
{
	bw1 <- load.bigWig(file.bw.plus);
	bw2 <- load.bigWig(file.bw.minus);

    total <- abs(bw1$mean*bw1$basesCovered) + abs(bw2$mean*bw2$basesCovered)

	unload.bigWig(bw1);
	unload.bigWig(bw2);
	
	return(total);
}

coverage_bed <- function(bedA, bed.map, range=50 )
{
    tmpA.bed <- tempfile(fileext=".bed");
    write.table( bedA, file=tmpA.bed, quote=F, row.names=F, col.names=F, sep="\t" );
    
    if(is.data.frame( bed.map))
    {
       tmp.map <- tempfile(fileext=".bed");
	   write.table( bed.map, file=tmp.map, quote=F, sep="\t", row.names=F, col.names=F );
    }
    else
       tmp.map <- bed.map;
    
cat( paste("bedmap --echo --min --max --mean --range ",  range, " --delim '\t'", tmpA.bed, tmp.map, "\n" ) );

    tbr <- read.table(pipe(paste("bedmap --echo --min --max --mean --range ",  range, " --delim '\t'", tmpA.bed, tmp.map ) ) );
    unlink(tmpA.bed);
	
	##important
    if (is.data.frame(bed.map)) unlink(tmp.map);

    return(tbr);
}


get_peak_from_MotifInfo <- function(path.motif, chrs=c("chr1", "chr22"), ncores=15 )
{    
	file.beds <- list.files(pattern = "\\.bed.gz$", path=path.motif);
	L0 <- mclapply (file.beds, function(motif)
	{
cat(paste0(path.motif, "/", motif), "\n")
		tb <- try(read.table(paste0(path.motif, "/", motif)));
		if(class(tb)=="try-error")
		{   
			cat("Error in ", paste0(path.motif, "/", motif), "\n");
			return(NULL);
		}	
		tb <- tb [tb$V1 %in% chrs, ]
		gc(reset=TRUE);
		return(tb);
	}, mc.cores=ncores);

	tb  <- do.call("rbind.data.frame", L0);
	tb  <- tb[order(tb$V1, tb$V2),]; 
	tb1 <- tb[tb$V5==1,];
	tb0 <- tb[tb$V5==0,];
	rm(tb);
	rm(L0);
	gc(reset=TRUE);
	
	return(list(tb1, tb0));
}

get_peak_from_PredictRdata <- function(path.rdata, ncores=15 )
{
	file.rdata <- list.files(pattern = "\\.rdata$", path=path.rdata);

	L0 <- mclapply (file.rdata, function(motif)
	{
		load(paste0(path.rdata, "/", motif))
		return(df_pred)
	},  mc.cores=ncores );

	tb  <- do.call("rbind.data.frame", L0)
	tb  <- tb[order(tb$V1, tb$V2),]; 
	tb1 <- tb[tb$pred_status==1,];
	tb0 <- tb[tb$pred_status==0,];
	rm(tb);
	rm(L0);

	tb1 <- tb1[tb1$pred_score > quantile(tb1$pred_score, 0.1),-c(9:10)];
	tb0 <- tb0[tb0$pred_score < quantile(tb0$pred_score, 0.9),-c(9:10)];
	tb1 <- tb1[order(tb1$V1, tb1$V2),]
	tb0 <- tb0[order(tb0$V1, tb0$V2),]

	return(list(tb1, tb0));
}
