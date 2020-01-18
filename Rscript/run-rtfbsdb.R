options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE);

## Read arguments from thwe web page.
arg_bsh_path    <- args[1]
arg_species     <- args[2]
arg_cisbp       <- args[3]
arg_genome_file <- args[4] 
arg_output_dir  <- args[5]

cpu_cores <- as.integer(args[6])
if (is.na(cpu_cores)) cpu_cores <- 1;

cat("------------ Parameters ------------- \n");
cat("Script Path:", arg_bsh_path, "\n");
cat("Species:", arg_species, "\n");
cat("CIS-BP data file:", arg_cisbp, "\n");
cat("Genome sequence:", arg_genome_file, "\n");
cat("TFBS output dir:", arg_output_dir, "\n");
cat("CPU cores:", cpu_cores, "\n");
cat("-------------------------------------\n ");

if(!file.exists(arg_cisbp))
	stop( paste("Can't find the CIS-BP database file.(", arg_cisbp, ")"));
if(!file.exists(arg_genome_file))
	stop( paste("Can't find the Genome sequece file(", arg_genome_file, ")"));
if(!file.exists(arg_output_dir))
	stop( paste("Can't find the output folder(", arg_output_dir, ")"));

make_bed<-function( df_bed, out_file)
{
    options(scipen=99, digits=4);
    file.tmp <- tempfile(fileext=".bed");
    write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
    system(paste( "sort-bed ", file.tmp, " | bgzip -f > ", out_file, ".gz", sep="") );
    unlink(file.tmp)
}

require(rtfbsdb);

file.twoBit_path <- arg_genome_file;

db <- CisBP.zipload(arg_cisbp, arg_species);
tfs <- tfbs.createFromCisBP(db);
df <- as.data.frame(table(tfs@tf_info$Motif_Type))
df <- df[df[,2]!=0,]

tf_info <- tfs@tf_info[, c("Motif_ID", "TF_Name")];
tf_info <- do.call("rbind", lapply(unique(tf_info$Motif_ID), function(x){ 
          idx <- which( as.character(tf_info$Motif_ID) == as.character(x) ); 
          if(NROW(idx)>4)
            return(data.frame(Motif_ID=x, TF_Name=paste(c(as.character(tf_info[idx, "TF_Name"])[1:4], "..."), collapse=",")))
          else if( NROW(idx)>1 )
             return(data.frame(Motif_ID=x, TF_Name=paste(as.character(tf_info[idx, "TF_Name"]), collapse=",")))
          else 
            return(data.frame(Motif_ID=x, TF_Name=as.character(tf_info[idx, "TF_Name"])))
        }));

write.table(tf_info, file=paste0(arg_output_dir, "/", arg_species, ".tf.info"), quote=F, row.names=F, col.names=F, sep="\t");

cat("Total Motif_Type", NROW(tf_info), "\n");

for (Motif_Type in df$Var1)
{
  cat("Scan TFBS site for Motif_Type=", Motif_Type, "\n");
  
  if (!file.exists(paste0(arg_output_dir, "/", arg_species, "_Motif_Type_", Motif_Type, ".rdata")))
  {
    tfs1 <- tfbs.createFromCisBP(db, motif_type=Motif_Type);
    dScan <- tfbs.scanTFsite( tfs1, arg_genome_file, ncores=as.numeric(cpu_cores) );
    save(dScan, file=paste0(arg_output_dir, "/", arg_species, "_Motif_Type_", Motif_Type, ".rdata") );
    rm(dScan);
    gc();
  }
}

bed.count = 0;
for (Motif_Type in df$Var1)
{
  if (file.exists(paste0(arg_output_dir, "/", arg_species, "_Motif_Type_", Motif_Type, ".rdata")))
  {
    load( paste0(arg_output_dir, "/", arg_species, "_Motif_Type_", Motif_Type, ".rdata") ) ;

    for (i in 1:NROW(dScan$summary))
    {
       cat("Output bed filr for Motif_Type=", dScan$summary$Motif_ID[i], "\n");

	   dScan.bed <- dScan$result[[i]];
	   out_file <- paste0(arg_output_dir, "/", as.character(dScan$summary$Motif_ID[i]), ".bed");
	   if(!file.exists(paste0(out_file, ".gz") ) )
           make_bed( dScan.bed, out_file);
        
       bed.count = bed.count + 1;
    }
  }
}

cat(bed.count, " bed files are generated.\n");
cat("Done!\n");
