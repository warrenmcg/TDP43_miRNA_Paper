#!/usr/bin/Rscript --slave

options(stringsAsFactors=F)
suppressMessages(library(optparse))

option.list <- list(
    make_option(c("-m","--mirna"), type="character", default=NULL, dest="mirna.names",
        metavar="file", help="[REQUIRED] Full path to file with miRNA names and miRBase accessions"),
    make_option(c("-c","--ctrdir"), type="character", default=NULL, dest="ctr.dir",
        metavar="directory", help="[REQUIRED] Full path to directory with control files."),
    make_option(c("-p","--ctrpattern"), type="character", default="-11_", dest="ctr.pattern",
        metavar="pattern", help="Regular expression to find control files [default %default]"),
    make_option(c("-t","--tumordir"), type="character", default=NULL, dest="tumor.dir",
	metavar="directory", help="[REQUIRED] Full path to directory with tumor files"),
    make_option(c("-u","--tumorpattern"), type="character", default="-0[1,2]_", dest="tumor.pattern",
        metavar="pattern", help="Regular expression to find tumor files [default %default]"),
    make_option(c("-o","--outprefix"), type="character", default="./output", dest="out.prefix",
        metavar="prefix", help="Output prefix for results [default %default]")
)

parser <- OptionParser(option_list=option.list)
opt <- parse_args(parser)

#mirna.names.file <- readline("What is the file you want to read mirna names from? ")
#control.dir <- readline("What is the directory for control files? ")
#control.pattern <- readline("What is the pattern to find control files? ")
#tumor.dir <- readline("What is the directory for case files? ")
#tumor.pattern <- readline("What is the pattern to find case files? ")
#output.prefix <- readline("What would you like the output prefix to be? ")

mirna.names.file <- opt$mirna.names
control.dir <- opt$ctr.dir
control.pattern <- opt$ctr.pattern
tumor.dir <- opt$tumor.dir
tumor.pattern <- opt$tumor.pattern
output.prefix <- opt$out.prefix

print("Finishing set-up...")

control.files <- dir(path=control.dir, pattern=control.pattern)
tumor.files <- dir(path=tumor.dir, pattern=tumor.pattern)

if (length(control.files) == 0) {
    stop(paste("Control directory", control.dir, "with pattern", control.pattern, "contains no files."))
} else if (length(tumor.files) == 0) {
    stop(paste("Tumor directory", tumor.dir, "with pattern", tumor.pattern, "contains no files."))
}

mirna.names.table <- read.table(mirna.names.file,header=T,sep="\t")

mirna.counts.df <- data.frame(mirBase.Accession=mirna.names.table[,1],mirna.name=as.factor(mirna.names.table[,2]))

avgAndSD <- data.frame(matrix(0, nrow(mirna.counts.df), ncol=5))
colnames(avgAndSD) <- c("mirna.name","mean.ctrl","sd.ctrl","mean.tumor","sd.tumor")
counter <- 1

updatemirnaCounts <- function(file,dir=NULL) {
	table <- read.table(file.path(dir,file),header=T,sep="\t")
	table <- table[,c("mirBase.Accession","raw.read.count")]
	mirna.counts.df <<- merge(mirna.counts.df, table, by=c("mirBase.Accession"), all=T)
	mirna.counts.df[is.na(mirna.counts.df)] <<- 0
	names(mirna.counts.df)[length(mirna.counts.df)] <<- file
}

calcSNR <- function(x) {
	mirna <- x[1]
	control.vals <- as.numeric(x[3:length(control.files)+2])
	tumor.vals <- as.numeric(x[(length(control.files)+3):length(x)])
	mean.ctrl <- mean(control.vals,na.rm=T)
	mean.tumor <- mean(tumor.vals,na.rm=T)
	sd.ctrl <- sd(control.vals,na.rm=T)
	sd.tumor <- sd(tumor.vals,na.rm=T)
	
	avgAndSD[counter,] <<- c(mirna,mean.ctrl,sd.ctrl,mean.tumor,sd.tumor)
	counter <<- counter + 1 
	
	if (!is.na(sd.ctrl) && !is.na(sd.tumor) && (sd.ctrl+sd.tumor) > 0) {
	    return((mean.tumor - mean.ctrl) / (sd.ctrl + sd.tumor))
	} else return(NaN)
}

print("Processing Control Files ...") 
invisible(sapply(control.files,updatemirnaCounts,dir=control.dir))
print("Processing Case Files ...")
invisible(sapply(tumor.files,updatemirnaCounts,dir=tumor.dir))

print("Calculating signal-to-noise ...")
mirna.snr <- apply(mirna.counts.df,1,calcSNR)

mirna.snr.table <- cbind(as.character(mirna.counts.df[,1]),mirna.snr)

print("Outputting files ...")
write.table(mirna.counts.df,file=paste(output.prefix,"_allmirna_mirnaCounts.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(mirna.snr.table,file=paste(output.prefix,"_mirnaCountSNR_rankedList.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(avgAndSD,file=paste(output.prefix,"_mirnaCountAvgAndSD.txt",sep=""),quote=F,sep="\t",row.names=F)
