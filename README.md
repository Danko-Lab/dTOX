dTOX
===============

Discriminative Transcription factor Occupancy eXtraction

Online computational gateway
-------------------------
We provide a computational gateway to run dTOX on GPU server, the users don't need to install any software, only upload the bigWig files and wait for the results, it is simple and easy. Please click the link to try this site:

https://dreg.dnasequence.org/

In the online service, dTOX only provide the service for 4 mammals: (human/hg19, mouse/mm10, horse/euuCab2, dog/canFam3). If your studies are not limited to these species, you have to install the dTOX pipeline and run your data locally. 

Before you run your data on the dREG gateway, please check the server status [here](https://github.com/Danko-Lab/dREG/blob/master/gateway-update.md).

![Hi](https://github.com/Danko-Lab/dTOX/raw/master/imgs/dreg-gateway.png?v=4&s=200 "dREG gateway")

### Important note for the Exchange email users:

The Exchange email system might quarantine all emails including the word  “password” or other sensitive stuffs in links. (https://technet.microsoft.com/en-us/library/aa997692(v=exchg.160).aspx).

Unfortunately, some emails from dREG gateway are quarantined by this spam policy. Usually these quarantined emails are not delivered to the email box, so they can not be checked in any email folders, including junk, spam or inbox. If you find the emails from dREG gateway are not delivered into your email box, please conect the administrator of your email system. For the Cornell email, please check this link:

https://it.cornell.edu/spam-control/log-quarantine-management-spam-control

Abstract
--------
Identification of the genomic regions that regulate transcription remains an important open problem.  We have recently shown that global run-on and sequencing (GRO-seq) with enrichment for 5-prime-capped RNAs reveals patterns of divergent transcription that accurately mark active transcriptional regulatory elements (TREs), including enhancers and promoters.  Here, we demonstrate that active TREs can be identified with comparable accuracy by applying sensitive machine-learning methods to standard GRO-seq and PRO-seq data, allowing TREs to be assayed together with transcription levels, elongation rates, and other transcriptional features, in a single experiment.  Our method, called discriminative Regulatory Element detection from GRO-seq (dREG), summarizes GRO-seq read counts at multiple scales and uses support vector regression to predict active TREs.  The predicted TREs are strongly enriched for marks associated with functional elements, including H3K27ac, transcription factor binding sites, eQTLs, and GWAS-associated SNPs.  Using dREG, we survey TREs in eight cell types and provide new insights into global patterns of TRE assembly and function. 

Data preparation: 
==========================

dTOX takes bigWig files with double strands for PRO-seq, ATAC-seq and DNase-I-seq, as the input. The bigWig files should follow 3 rules:

1) Each read is mapped at 5’ (GRO-seq) or 3’ (PRO-seq) position (point mode) , not mapped to a continuous region starting from 5’ or 3’.  This is different with the software Tfit.

2) Only positive values or only negative values in each strand, no mixture.

3) No normalization

As for how to generate bigWig files from fastq data, please refer to https://github.com/Danko-Lab/proseq2.0/.


Installation instructions: 
==========================

The source code and models of dTOX will be availiable on GitHub (https://github.com/Danko-Lab/dTOX).  

Supported OS:
-------------
Linux and Mac OSX are currently supported.

Required software
-----------------
* bedops (http://bedops.readthedocs.org/en/latest/index.html)
* R (http://www.r-project.org/)
* bigWig R package (https://github.com/andrelmartins/bigWig; will be public very soon).

This software is already installed on many UNIX systems.  Users can install the most appropriate version of these files for Ubuntu using: 

    sudo apt-get install r-base-core
    sudo apt-get install libssl1.0.0 libssl-dev

Users who are not sure how to install the proper dependencies on their system should check with their system administrator for help.  

dTOX also has several dependencies within R. These include **data.table**, **e1071**, **mvtnorm**, **parallel**, **randomForest**, **rmutil**, **rphast**, and **snowfall**. These packages are all availiable on the CRAN repository.  For convenience, users can install these packages using the makefile:

    make R_dependencies

If users run into any problems they should contact the package author for assistance.

Install dTOX
------------
Users should change to the directory containing this README.md file, and can then install dTOX by typing the following:

(1a) Install R dependencies

    make R_dependencies

(1b) Install dTOX

    make dTOX

(1c) Install Rgtsvm if you have GPU nodes.

    make Rgtsvm


Get the dTOX models
-------------------
Pre-trained models that can be used to predict TF binding status across the genome are availiable in mammals. Download the newest models and TFBS data set from here:
<span style="color:blue"> *ftp://cbsuftp.tc.cornell.edu/danko/hub/dTOX/* </span>


Usage instructions:
===================

dTOX provides one solutions to identify TF bound or unbound status for 3 sequence data,  PRO-seq, ATAC-seq, and DNase-I-seq. 

The first solution implemented in the early package, is to ***predict dREG scores*** and detect the broad dREG peaks with the aid of Perl program. In order to identify narrow peak, these broad peaks need to be refined using [dREG-HD package](https://github.com/Danko-Lab/dREG.HD).

The second solution implements the ***peak calling*** function using the dREG scores based on the imporved SVR model. Compared with the broad peaks in the first solution, this solution generates the narrow peaks with peak score, probability, center position. Although this solution simplies the dREG process, it relies on GPU computing nodes to acceleratethe computational speed. If GPU computing resource is not available for you, please try our online computational gateway (https://dreg.dnasequence.org/). 

In this section, we will introduce new solution following by old one.

## 1) TFBS calling

To use this solution, type: 

    bash run_rtfbsdb.bsh species cisbp_zip genome_2bit TFBS_output_dir [cpu_cores] 

    species	        -- any names for the species.
    cisbp_zip	    -- The zip file including motifs list and PWM file. (downloaded from cis-bp).
    genome_2bit     -- The 2 bit file of genome 
    TFBS_output_dir -- The path to output the TFBS bed files for each motif define in cisbp_zip
    cpu_cores       -- [optional, default=1] indicating how many CPU cores can be used.

How to download motif information from the CIS-BP database (http://cisbp.ccbr.utoronto.ca/)

![Hi](https://raw.githubusercontent.com/Danko-Lab/dTOX/master/imgs/cisbp-download.png?v=4&s=600 "Download from CIS-BP")


## 2) TF binding predict

To use this solution, type: 

    bash run_dTOX.bsh hg19 seq_type plus_strand_bw minus_strand_bw out_prefix gpu_cores cpu_cores

    Species         -- Currently only hg19 
    Sequeue Type    -- ATAC-seq, DNase-1-seq or PRO-seq 
    plus_strand.bw	-- Seqence data (plus strand).  Read counts (not normalized) formatted as a bigWig file.
    minus_strand.bw	-- Seqence data (minus strand). Read counts (not normalized) formatted as a bigWig file.
    out_prefix      -- The prefix of the output file.
    gpu_cores       -- [optional, default=1] indicating how many GPU cores can be used.
    cpu_cores       -- [optional, default=16] indicating how many CPU cores can be used.


For example, to run dTOX on the ATAC-seq data, use:

    bash run_dTOX.bsh hg19 ATAC-seq  atac.plus.bw atac.minus.bw atac.out 1 16

One file below is generated in this solution:  

1. <out_prefix>.dTOX.bound.bed.gz:

    reduced peak information, only including peak position, max score.

**Notice:** 

(1) That command takes 8~20 hours to execute on NVIDA K80 GPU using Rgtsvm package. Due to very long computational time, we don't suggest to run peak calling on CPU nodes, even in parallel mode.


# Document

dTOX is an R package, and that provides some additional flexibility for users familiar with R. Currently you can get details about each R function from the dREG manual (https://github.com/Danko-Lab/dREG/blob/master/dREG-manual.pdf).  We are actively working to document each function in the package.  

How to cite
===================
(1)  Identification of transcription factor binding patterns using genome-wide accessibility and transcription



