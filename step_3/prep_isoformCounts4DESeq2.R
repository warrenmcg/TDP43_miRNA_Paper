suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
suppressMessages(library(optparse))
options(stringsAsFactors=F)

option.list <- list(
    make_option(c("-c","--ctrdir"), type="character", default=NULL, dest="ctr.dir",
        metavar="directory", help="[REQUIRED] Full path to directory with control files."),
    make_option(c("-p","--ctrpattern"), type="character", default="-11\\.", dest="ctr.pattern",
        metavar="pattern", help="Regular expression to find control files [default %default]"),
    make_option(c("-t","--tumordir"), type="character", default=NULL, dest="tumor.dir",
        metavar="directory", help="[REQUIRED] Full path to directory with tumor files"),
    make_option(c("-u","--tumorpattern"), type="character", default="-0[1,2]\\.", dest="tumor.pattern",
        metavar="pattern", help="Regular expression to find tumor files [default %default]"),
    make_option(c("-o","--outprefix"), type="character", default="./output", dest="out.prefix",
        metavar="prefix", help="Output prefix for results [default %default]")
)

parser <- OptionParser(option_list=option.list)
opt <- parse_args(parser)

control.dir <- opt$ctr.dir
control.pattern <- opt$ctr.pattern
tumor.dir <- opt$tumor.dir
tumor.pattern <- opt$tumor.pattern
output.prefix <- opt$out.prefix

print("Finishing set-up...")

control.files <- dir(path=control.dir, pattern=control.pattern, full.names=T)
tumor.files <- dir(path=tumor.dir, pattern=tumor.pattern, full.names=T)

data.files <- c(control.files, tumor.files)

split <- as.data.frame( strsplit( data.files, "/" ) ) 
file.names <- as.character(split[nrow(split),])
prefixes <- as.character( as.data.frame( strsplit ( file.names, ".", fixed=T) )[1,] )

first.file <- data.files[1]

first.data <- read.table(first.file, header=T, sep="\t", stringsAsFactors=F)
first.data <- first.data[,1:2]

count.data <- matrix(nrow=nrow(first.data), ncol=length(prefixes))

prepData <- function(index) {
	file <- data.files[index]
	table <- read.table(file, header=T, sep="\t", stringsAsFactors=F)
	table <- table[,1:2]
	stopifnot(identical(table[,1], first.data[,1]))
	table[,2]
}

count.data.mat <- sapply(seq_along(file.names), prepData)
colnames(count.data.mat) <- prefixes
rownames(count.data.mat) <- first.data[,1]
count.data.mat <- round(count.data.mat)

colData <- data.frame(row.names=prefixes)
colData$condition <- factor(c(rep("control",length(control.files)),rep("tumor",length(tumor.files))), levels = c("control","tumor"))

dds <- DESeqDataSetFromMatrix(count.data.mat, colData, ~ condition)
dds <- DESeq(dds, parallel=T, BPPARAM = MulticoreParam(workers=8))
results <- results(dds)
sig.results <- results[which(results$padj<=0.05),]

full.results <- cbind(isoform_id=row.names(results), as.data.frame(results))
full.sig.results <- full.results[which(full.results$padj<=0.05),]

rank.results <- full.results[,c("isoform_id", "baseMean", "log2FoldChange","padj")]
rank.results$logQval <- -log10(rank.results$padj)
rank.results$direction <- ifelse(rank.results$log2FoldChange >0, "up", ifelse(rank.results$log2FoldChange < 0, "down", "no change"))
rank.results$signedLogQ <- apply(rank.results, 1, function(x) {
					options(digits=10)
					logqval <- as.numeric(x[5])
					direction <- as.character(x[6])
					if (is.na(direction) || is.na(logqval)) signedlogq <- 0 else
					if (direction == "up") signedlogq <- logqval else
					if (direction == "down") signedlogq <- -1 * logqval
					else signedlogq <- 0
					signedlogq
				})
rank.signed.summary <- rank.results[,c("isoform_id","signedLogQ")]
rank.signed.summary <- rank.signed.summary[order(rank.signed.summary$signedLogQ, decreasing=T),]

write.table(full.results, file=paste(output.prefix,"transcripts_allResults.txt",sep="_"),quote=F,row.names=F,sep="\t")
write.table(full.sig.results, file=paste(output.prefix,"transcripts_sigResults.txt",sep="_"),quote=F,row.names=F,sep="\t")
write.table(rank.signed.summary, file=paste(output.prefix, "transcriptSignedLogQ_Rank.txt",sep="_"),quote=F,row.names=F,col.names=F,sep="\t")
save(list=c("dds", "results", "sig.results","full.results","full.sig.results","rank.results","rank.signed.summary"),file=paste(output.prefix,"isoform_DESeq2_data.Rdata",sep="_"))
print(sessionInfo())
