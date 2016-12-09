suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
options(stringsAsFactors=F)

option.list <- list(
    make_option(c("-m","--mirna"), type="character", default=NULL, dest="mirna.counts",
        metavar="file", help="[REQUIRED] Full path to file with miRNA counts matrix."),
    make_option(c("-p","--ctrpattern"), type="character", default="[-.]11_", dest="ctr.pattern",
        metavar="pattern", help="Regular expression to find control files [default %default]"),
    make_option(c("-u","--tumorpattern"), type="character", default="[-.]0[12]_", dest="tumor.pattern",
        metavar="pattern", help="Regular expression to find tumor files [default %default]"),
    make_option(c("-o","--outprefix"), type="character", default=NULL, dest="output.prefix",
        metavar="prefix", help="Prefix for output files [default is same prefix as input files]")
)

parser <- OptionParser(option_list=option.list)
opt <- parse_args(parser)

#file <- readline("What is the miRNA counts you want to work with? ")
file <- opt$mirna.counts
ctr.pattern <- opt$ctr.pattern
tumor.pattern <- opt$tumor.pattern

if (is.null(opt$output.prefix)) {
    output.prefix <- as.character(unlist(strsplit(file, "_")))[1] 
} else output.prefix <- opt$output.prefix
print(paste("Out prefix is", output.prefix))

count.data <- read.table(file, header=T, sep="\t", stringsAsFactors=F)
prefixes <- names(count.data)[c(-1,-2)]
controls <- prefixes[grepl(ctr.pattern, prefixes)]
cases <- prefixes[grepl(tumor.pattern, prefixes)]
row.names(count.data) <- count.data[,1]
count.data <- as.matrix(count.data[,c(-1,-2)])
count.data <- round(count.data)

colData <- data.frame(row.names=prefixes)
colData$condition <- as.factor(c(rep("control",length(controls)),rep("tumor",length(cases))))
levels(colData$condition) <- c("control","tumor")

dds <- DESeqDataSetFromMatrix(count.data, colData, ~ condition)
dds <- DESeq(dds)
results <- results(dds, tidy=T)
names(results)[1] <- "mirBaseAccession"

#results$padj[is.na(results$padj)] <- 1

sig.results <- results[which(results$padj<=0.05),]

save(list=c("dds", "results", "sig.results"),file=paste(output.prefix, "DESeq2_data.Rdata", sep="_"))
write.table(results, file=paste(output.prefix,"miRNA_allResults.txt",sep="_"), quote=F, sep="\t", row.names=F)
write.table(sig.results, file=paste(output.prefix,"miRNA_sigResults.txt",sep="_"),q uote=F, sep="\t", row.names=F)
print(sessionInfo())
