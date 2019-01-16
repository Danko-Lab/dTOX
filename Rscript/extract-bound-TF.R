options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE);

## Read arguments from thwe web page.
arg_species       <- args[1]
arg_bound_file    <- args[2]
arg_TF            <- args[3]
arg_output        <- args[4] 

cat("------------ Parameters ------------- \n");
cat("Input bound regions:", arg_bound_file, "\n");
cat("Species", arg_species, "\n");
cat("TF:", arg_TF, "\n");
cat("Output file:", arg_output, "\n");
cat("-------------------------------------\n ");

if( arg_species != "hg19" && arg_species != "mm10" )
	stop( "Only hg19 and mm10 are available now.");
	
if(!file.exists(arg_bound_file))
	stop( paste("Can't find the bound file.(", arg_bound_file, ")"));

if(file.exists(arg_output))
	stop( paste("Output file is existing(", arg_output, ")"));

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

extract_TF_bound<-function(arg_species, arg_bound_file, arg_TF, arg_output)
{
	require(rtfbsdb);
	require(parallel);
	require(data.table);

	if(arg_species=="hg19")
		db <- CisBP.extdata("Homo_sapiens")

	if(arg_species=="mm10")
		db <- CisBP.extdata("Mus_musculus");

	tfs <- tfbs.createFromCisBP(db, tf_name=arg_TF);
    if(is.null(tfs@tf_info) || is.null(tfs@tf_info) || NROW(tfs@tf_info)==0 )
  	   stop( paste("NO motif is belonging to TF(", arg_TF, ")"));
    
	motif_list <- unique(tfs@tf_info[, c("Motif_ID")]);
	
	motif_bound <- rbindlist( mclapply( as.character(motif_list), function(x){ 
		   str.cmd <- paste0( "zcat ", arg_bound_file, " | awk 'BEGIN{OFS=\"\\t\"}{ if ($4==\"", x, "\") print $0 }' - ");
		   tb <- read.table(pipe(str.cmd), stringsAsFactors=FALSE);
		   show(head(tb));
		   return(as.data.frame(tb))
			}, mc.cores=16));

	write.bed(as.data.frame(motif_bound), arg_output, compress=T);

	cat("Done!\n");
}

extract_TF_bound( arg_species, arg_bound_file, arg_TF, arg_output )
