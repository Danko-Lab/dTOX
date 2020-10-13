library(parallel)
library(sqldf)

options("scipen"=100, "digits"=4);

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


get_motifinfo_on_chr <- function(path.motif, chr=c("chr1"), path.output, ncores=12 )
{
cat("chr=", chr, "\n");
    
    file.beds <- list.files(pattern = "\\.bed.gz$", path=path.motif);
    L0 <- mclapply (file.beds, function(motif)
    {
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
    
    file.tfbs.chrx  = tempfile(tmpdir = "./", fileext =".bed")
    write.bed(tb, file.tfbs.chrx, compress=FALSE )

    rm(tb);
    rm(L0);
    gc(reset=TRUE);

    file.beds <- list.files(pattern = "\\.bed.gz$", path=path.motif);
    L0 <- mclapply (file.beds, function(motif)
    {
cat(paste0(path.motif, "/", motif), "\n")
        file.motif = paste0(path.motif, "/", motif); 
        str.cmd = paste0("zcat ", file.motif, " | awk '{if ($1==\"", chr, "\") print $0 }' -");
        tb <- try(read.table(pipe(str.cmd)));


        file.bed <- paste0(path.output, "/ext-", basename(motif), ".bed");
browser();

        tb_ext <- coverage_bed( tb, file.tfbs.chrx, 25);
        
        write.table( tb_ext, file=file.bed, quote=F, row.names=F, col.names=F, sep="\t", append = TRUE);
    }

    unlink(file.tfbs.chrx);  

    gc(reset=TRUE); 
    return(1);
}

if(0)
{
path.motif="./hg19-tfbs";
path.output="./hg19-ext";

system(paste("rm", path.output, "/ext-*.bed")
system(paste("rm", path.output, "/ext-*.bed.gz")


for(chr in paste0("chr", 22:1))
   get_motifinfo_on_chr( path.motif, chr, path.output, ncores=12 )

system("bgzip *.bed")

}

if(1)
{
path.motif="./mm10-tfbs";
path.output="./mm10-ext";

system(paste("rm", path.output, "/ext-*.bed")
system(paste("rm", path.output, "/ext-*.bed.gz")

for(chr in paste0("chr", 19:1))
   get_motifinfo_on_chr( path.motif, chr, path.output, ncores=12 )

system("bgzip *.bed")
}