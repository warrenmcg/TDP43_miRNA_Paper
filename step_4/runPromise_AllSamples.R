#!/usr/bin/Rscript

args <- commandArgs(TRUE)
options(stringsAsFactors=F)
library("Roleswitch") # needed package to do the PROMiSe analysis

if(length(args) != 4 | args[1] != "file_manifest.txt") {
	stop("This script has four arguments - a file manifest list, a predictions table, the miRNA counts directory, the RNASeq counts directory")
}

print("Loading data files …")
# assumes you've either given the absolute path or the correct relative path to the files
file_manifest <- read.table(file.path(args[1]),header=T,sep="\t")
miRNA_dir <- file.path(args[3])
rna_dir <- file.path(args[4])

miRNA_files <- file_manifest[file_manifest$Platform.Type=="miRNASeq",]
miRNA_isoform_files <- miRNA_files[grep(".isoform.quantification.txt$",miRNA_files$File.Name),]

#load the miRNA-mRNA target pair matrix	
targetPairsDF <- read.table(file.path(args[2]), header=TRUE, sep="\t", stringsAsFactors=FALSE)
targetMat <- data.matrix(targetPairsDF[,2:ncol(targetPairsDF)])
row.names(targetMat) <- targetPairsDF[,1]

dir.create("output")

##############################################################################
# This method normalizes the rows into a doubly stochastic matrix
normlize.p <- function(p) {
	
	p <- apply(p, 1, function(x)x/sum(x))
	p[is.na(p)] <- 0
	p <- apply(p, 2, function(x)x/sum(x))
	p[is.na(p)] <- 0
	p	
}

# This method takes the two files as input, as well as a matrix of predicted miRNA/mRNA interactions, and performs PROMiSe ("roleswitch") after some clean-up
run_promise <- function(miRNAFile,rnaFile,targetMatrix) {

	#load the miRNA expression profile
	miRNA_table <- read.table(miRNAFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
	miRNA_EP <- data.matrix(miRNA_table[,4])
	row.names(miRNA_EP) <- miRNA_table[,1]

	#load the target expression profile
	target_table <- read.table(rnaFile,header=TRUE,sep="\t",stringsAsFactors=FALSE)
	target_EP <- data.matrix(target_table[,2])
	row.names(target_EP) <- target_table[,1]

	# remove miRNAs and mRNAs not covered in the expression profile or data matrix
	if(dim(miRNA_EP)[1] < dim(targetMatrix)[2]) {
		targetMatrix <- targetMatrix[,which(colnames(targetMatrix) %in% row.names(miRNA_EP))]
		miRNA_EP <- as.matrix(miRNA_EP[which(row.names(miRNA_EP) %in% colnames(targetMatrix)),])
	} else { 
		miRNA_EP <- as.matrix(miRNA_EP[which(row.names(miRNA_EP) %in% colnames(targetMatrix)),])
		targetMatrix <- targetMatrix[,which(colnames(targetMatrix) %in% row.names(miRNA_EP))]
	}

	if(dim(target_EP)[1] < dim(targetMatrix)[1]) {
		targetMatrix <- targetMatrix[which(row.names(targetMatrix) %in% row.names(target_EP)),]
		target_EP <- as.matrix(target_EP[which(row.names(target_EP) %in% row.names(targetMatrix)),])
	} else { 
		target_EP <- as.matrix(target_EP[which(row.names(target_EP) %in% row.names(targetMatrix)),])
		targetMatrix <- targetMatrix[which(row.names(targetMatrix) %in% row.names(target_EP)),]
	}

	rs.pred <- roleswitch(target_EP, miRNA_EP, targetMatrix)
	rs.pred$p.x[is.na(rs.pred$p.x)] <- 0
	rs.pred$p.z[is.na(rs.pred$p.z)] <- 0
	rs.pred$p.xz <- normlize.p(rs.pred$p.x * rs.pred$p.z)
	save(rs.pred, file=paste("output/", sample_name,"_ProMISe.RData",sep=""))
}
##############################################################################
	
#rs.pred.list <- vector(mode="list",nrow(miRNA_isoform_files))

print("Processing profiles …")
for (i in 1:nrow(miRNA_isoform_files))
{
	sample_name <- miRNA_isoform_files[i,5]
	miRNA_counts <- file.path(paste(miRNA_dir, sample_name,"_miRNAcounts.txt",sep=""))
	rnaSeq_file <- file.path(paste(rna_dir, sample_name,".rsem.isoforms.normalized_results",sep=""))
	if (file.exists(miRNA_counts) && file.exists(rnaSeq_file)) {
		run_promise(miRNA_counts, rnaSeq_file, targetMat)	
	}
	else{
		sink(paste("output/", sample_name,"_ProMISe.txt",sep=""),append=F,split=F)
		print(c(sample_name, "no matched miRNA and RNA profiles found"))
		sink()
	}
}
#save.image(file=file.path(getwd(),"promise_results.RData"))
q()