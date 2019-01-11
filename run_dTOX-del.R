options("scipen"=100, "digits"=4)

#R --vanilla --slave --args $bsh_path $genomeId $seqType $bwPL $bwMN $out $ncores $gpu $bwCheck /run_TOX.R

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Read arguments from thwe web page.

arg_path_bash  <- args[1]
arg_model      <- args[2]
arg_tfbs_dir   <- args[3]
arg_plus_path  <- args[4]
arg_minus_path <- args[5]
arg_out_prefix <- args[6]
arg_ncores     <- as.integer(args[7])
arg_gpu_used   <- as.logical(args[8])

out.log  <- paste(arg_out_prefix, "dTOX.log", sep=".")

if (is.na(arg_ncores)) 
    arg_ncores <- 1;
if (is.na(arg_gpu_used)) 
    arg_gpu_used <- TRUE;

check_bw <- TRUE;
if(NROW(args)>=9)
   check_bw <- as.logical(args[9]);

cat("Home path:", arg_path_bash, "\n");
cat("Model file:", arg_model, "\n");
cat("TFBS dir:", arg_tfbs_dir, "\n");
cat("Bigwig(plus):", arg_plus_path, "\n");
cat("Bigwig(minus):", arg_minus_path, "\n");
cat("Output:", arg_out_prefix, "\n");
cat("ncores:", arg_ncores, "\n");
cat("GPU:", arg_gpu_used, "\n");
cat("Check BigWigs", check_bw, "\n");

source(paste(arg_path_bash, "Rscript/dtox-routines.R", sep="/"));

cat("1) -------- Checking bigWig files\n");
if(check_bw)
{
  b1 <- check_bigwig(arg_plus_path, strand="+", out.file = out.log );
  b2 <- check_bigwig(arg_minus_path, strand="-", out.file = out.log );
  if( !b1 || !b2 )
  {
    cat("Warning: bigWig files maybe not meet the requirements. See dREG requirement in https://github.com/Danko-Lab/dREG");
    stop("Stop");
  }
}

## Now scan all positions in the genome ...
cat("2) -------- Prediction\n");
t.run <- system.time( r <- predictH0( arg_path_bash, arg_model, arg_plus_path, arg_minus_path,  arg_tfbs_dir, cpu_cores=arg_ncores, use_rgtsvm=arg_gpu_used  ) )
r <- fix_TF_name0( r, arg_tfbs_dir );
r$V5 <- round(r$V5, digits=4);
r$pred <- round(r$pred, digits=4);

cat("3) -------- outputting result\n");

# columns in r
# chr, start, pos, MOTIF_ID, score, strand, TF, predict_score, predict_status

out.chrom.info <- paste(arg_out_prefix, "chrom.info", sep=".")
out.full.bed<- paste(arg_out_prefix, "dTOX.full.bed", sep=".")
out.rtfbsdb.bed<- paste(arg_out_prefix, "dTOX.rtfbsdb.bed", sep=".")
out.bound.bed<- paste(arg_out_prefix, "dTOX.bound.bed", sep=".")

make.chrom.info(arg_plus_path, arg_minus_path, out.chrom.info);
make_index_gz( r, out.full.bed);

make_index_gz( r[,c(1,2,3,5)], out.rtfbsdb.bed);
make_index_gz( r[,c(1,2,3,8)], out.bound.bed);

system( paste("tar -cvzf ", arg_out_prefix, ".tar.gz", " ", arg_out_prefix, ".dTOX.*", sep="") );
cat("Result:", paste(arg_out_prefix, ".tar.gz", sep=""), "\n");

cat("Done!\n");
