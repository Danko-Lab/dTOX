library(parallel)
library(sqldf)

write.bed<-function ( df.bed, file.bed, compress=FALSE )
{
    options("scipen"=100, "digits"=4);
    temp <- tempfile(fileext=".bed");
    write.table( df.bed, file=temp, quote=F, row.names=F, col.names=F, sep="\t");
    if(compress)
        system(paste0("sort-bed ", temp,  " | bgzip > ",  file.bed ))
    else
        system(paste0("sort-bed ", temp,  " > ",  file.bed ));
        
    invisible(unlink(temp));
}


get_center_from_MotifInfo <- function(path.motif, chr=c("chr1"), path.output, ncores=12, file.prefix.bed="tfbs.occupy" )
{
cat("chr=", chr, "\n");
    
    file.beds <- list.files(pattern = "\\.bed.gz$", path=path.motif);
    L0 <- mclapply (file.beds, function(motif)
    {
cat(paste0(path.motif, "/", motif), "\n")
        file.motif = paste0(path.motif, "/", motif); 
        str.cmd = paste0("zcat ", file.motif, " | awk '{if ($1==\"", chr, "\") print $0 }' -");
        tb <- try(read.table(pipe(str.cmd)));
        if(class(tb)=="try-error")
        {   
            cat("Error in ", paste0(path.motif, "/", motif), "\n");
            return(NULL);
        }    
        #tb <- tb [tb$V1 %in% chr, ]
        gc(reset=TRUE);
        return(tb);
    }, mc.cores=ncores);

    tb  <- do.call("rbind.data.frame", L0);
    rm(L0);
    gc(reset=TRUE);

#browser();
    tb$center<- round((tb$V2+tb$V3)/2);
    tb2 <- sqldf("select V1, center, group_concat(V4) as V4, max(V5) as V5, min(V7) as V7, max(V8) as V8, avg(V9) as V9 from tb group by V1, center order by V1, center", method = "raw") 
    colnames(tb2) <- c("chr", "start", "motifs", "scores", "min_ext25", "max_ext25", "avg_ext25" ); 
    tb2$stop <- tb2$start+1;
    
    tb2 <- tb2[,c("chr", "start", "stop", "scores", "min_ext25", "max_ext25", "avg_ext25", "motifs")]
    
cat("NROW(tb)", NROW(tb), "==>", "NROW(tb2)", NROW(tb2), "\n");
    rm(tb);
    gc(reset=TRUE);

    data.idx <- unique( c( seq(1, NROW(tb2), 1000*1000),NROW(tb2)+1 ) );
    for(i in 1:(NROW(data.idx)-1))
    {
        file.bed.gz <- paste( paste0(path.output, "/", file.prefix.bed), chr, i, "bed.gz", sep=".");
        write.bed( tb2[data.idx[i]:(data.idx[i+1]-1),], file.bed.gz, compress=TRUE);
    }    
    
    rm(tb2);
        gc(reset=TRUE); 
    return(1);
}


#for(chr in paste("chr", 1:21, sep=""))
#    {get_center_from_MotifInfo("/local/ftp/pub/hub/dTOX/TFBS/motif_rtfbsdb_hg19_ext", chr=chr, file.prefix.bed="hg19.tfbs.occupy");gc()}



#tb <- read.table("Spur_3.1.LinearScaffold.fa.fai", stringsAsFactors =FALSE)
#for(i in 1:NROW(tb))
#   get_center_from_MotifInfo("/workdir/zw355/spur_rtfbsdb_ext", chr=tb[i,1], "/workdir/zw355/spur_rtfbsdb_ext2", file.prefix.bed="spur.tfbs.occupy")

if(1)
{
path.ext="./tfbs_ext"
path.output="./tfbs_ext2"
ncores=12
file.prefix.bed="mm10.tfbs.occupy"

for(chr in paste("chr", 19:1, sep=""))
    {get_center_from_MotifInfo(path.ext, chr=chr, path.output, ncores=ncores, file.prefix.bed);gc()}
}


