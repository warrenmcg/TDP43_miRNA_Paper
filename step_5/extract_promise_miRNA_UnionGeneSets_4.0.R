suppressMessages(library(optparse))

option.list <- list(
    make_option(c("-p","--probthreshold"), type="numeric", default=0, dest="prob.threshold",
        metavar="num [0,1]", help="Probability threshold above which a target is considered predicted (default x > %default)"),
    make_option(c("-n","--numthreshold"), type="character", default="all", dest="num.threshold",
        metavar="num [1,num samples]", help="Threshold number of samples where a target is predicted for interaction to be included (default is 'all' samples)"),
    make_option(c("-o","--optimizer"), type="character", default="N", dest="optimize",
        metavar="y|n", help="Flag to determine whether to optimize the gene sets or not (default %default)"),
    make_option(c("-d","--datadir"), type="character", default="./output", dest="data.dir",
        metavar="pattern", help="Directory containing the ProMISe RData files [default is %default]"),
    make_option(c("-x","--prefix"), type="character", default=NULL, dest="out.prefix",
        metavar="prefix", help="[REQUIRED] Output prefix for results"),
    make_option(c("-a","--accession"), type="character", default=NULL, dest="accession",
	metavar="file", help="[REQUIRED] Filename that converts miRBase accession to miRNA name")
)

parser <- OptionParser(option_list=option.list)
opt <- parse_args(parser)

#prob.threshold <- as.numeric(readline("What is the probability threshold, above which you want to keep target predictions [0,1]? "))
#num.threshold <- as.integer(readline("What is minimum number of samples you wish to have an interaction to keep it [1,number of samples]? "))
#optimize <- readline("Do you wish to simply pick the maximum number of viable gene sets [y/n]? ")
#data_dir <- readline("Where is the RData that you want to load? ")
#prefix <- readline("What is the prefix you want for the output file? ")
#accession2name.file <- readline("Where is the file that converts miRBase accession to miRNA name (please use full directory)? ")
old_dir <- getwd()

prob.threshold <- opt$prob.threshold
optimize <- opt$optimize
data_dir <- opt$data.dir
filenames = dir(path=data_dir, pattern="*.RData")
firstFile <- filenames[1]
otherFiles <- filenames[-1]
if (opt$num.threshold == "all") num.threshold <- length(filenames) else
    num.threshold <- as.integer(opt$num.threshold)
prefix <- opt$out.prefix
accession2name.file <- opt$accession

acc2name.table <- read.table(accession2name.file, header=T, sep="\t", stringsAsFactors=F)
row.names(acc2name.table) <- acc2name.table[,1]

gene.sets <- vector("list",nrow(acc2name.table))
names(gene.sets) <- row.names(acc2name.table)

load(file.path(data_dir, firstFile))
promiseTarget.JointMatrix <- t(rs.pred$p.xz)
gene.names <- rownames(promiseTarget.JointMatrix)

for(i in 1:nrow(promiseTarget.JointMatrix)){
	miRNA <- gene.names[i]
	promise_scores <- promiseTarget.JointMatrix[i,]
	gene.sets[[miRNA]] <- list(transcripts=names(which(promise_scores > prob.threshold)) )
	gene.sets[[miRNA]]$counts <- rep(1, length(gene.sets[[miRNA]]$transcripts))
	names(gene.sets[[miRNA]]$counts) <- gene.sets[[miRNA]]$transcripts
}

update.scores <- function(fileName) {
	load(file.path(data_dir, fileName))
	promiseTarget.JointMatrix <- t(rs.pred$p.xz)
	gene.names <- rownames(promiseTarget.JointMatrix)
	
	for(i in 1:nrow(promiseTarget.JointMatrix)){
		miRNA <- gene.names[i]
		if(!(miRNA %in% names(gene.sets))) {
			stop(paste(miRNA, "not in gene.sets!"))
		}
		promise_scores <- promiseTarget.JointMatrix[i,]
		new_set <- names(which(promise_scores > prob.threshold))
		
    if(is.null(gene.sets[[miRNA]])){
      gene.sets[[miRNA]]$transcripts <<- new_set
      counts <- rep(1, length(gene.sets[[miRNA]]$transcripts))
      names(counts) <- gene.sets[[miRNA]]$transcripts
      gene.sets[[miRNA]]$counts <<- counts
      next()
    }
		
		old_set <- gene.sets[[miRNA]]$transcripts
		gene.sets[[miRNA]]$transcripts <<- union(old_set, new_set)
    old.counts <- gene.sets[[miRNA]]$counts
		new.counts <- rep(1, length(gene.sets[[miRNA]]$transcripts))
		names(new.counts) <- gene.sets[[miRNA]]$transcripts
		update.set <- intersect(old_set, new_set)
		#browser()
    if(length(update.set) > 0) {
			update.counts <- old.counts[update.set]
			if(length(update.counts) < 1) {
				browser()
				stop(paste(miRNA, "apparently does not have updated counts!"))
			}
			new.counts[update.set] <- update.counts + c(1)
		}
		gene.sets[[miRNA]]$counts <<- new.counts
	}
}

getAllMaxAndMin <- function(list, upperBound) {
	maxMin.df <- data.frame(nummax=integer(length=upperBound), min=integer(length=upperBound), nsets=integer(length=upperBound))
	for(i in 0:(upperBound-1)) {
		targetCounts <- sapply(list, function(x) length(which(x$counts > i)))
		nonzero.targets <- targetCounts[which(targetCounts>0)]
		num.max <- length(nonzero.targets[which(nonzero.targets>1000)])
		min <- min(nonzero.targets)
		length <- length(nonzero.targets[which(nonzero.targets<=1000)])
		maxMin.df[i+1,] <- c(num.max, min, length)
	}
	return(maxMin.df)
}

getFinalGenesets <- function(list, threshold) {
	final.gene.sets <- sapply(list, function(x) x$transcripts[which(x$counts >= threshold)])
	return(final.gene.sets)
}

output.genesets <- function(name) {
	vector <- t(c(name, acc2name.table[name,2], final.gene.sets[[name]]))
	write.table(vector, file=paste(prefix, opt$num.threshold, prob.threshold, "UnionGeneSetsWithCounts.gmt", sep="_"), append=T, sep="\t", quote=F, row.names=F, col.names=F)
}

sapply(otherFiles, update.scores)
maxmin.df <- getAllMaxAndMin(gene.sets, length(filenames))

if(grepl("^[Y/y]",optimize)){
	candidate.rows <- maxmin.df[which(maxmin.df[,3]==max(maxmin.df[,3])),]
	row <- which(candidate.rows[,1]==max(candidate.rows[,1]))
	num.threshold <- as.integer(rownames(candidate.rows)[row])
}

final.gene.sets <- getFinalGenesets(gene.sets, num.threshold)

prob.threshold <- as.character(prob.threshold)
num.threshold <- as.character(num.threshold)
genesets.file.name <- paste(prefix, opt$num.threshold, prob.threshold, "UnionGeneSetsWithCounts.gmt", sep="_")
if (file.exists(genesets.file.name)) {
	file.remove(genesets.file.name)
	sapply(names(final.gene.sets), output.genesets)
} else sapply(names(final.gene.sets), output.genesets)
setwd(old_dir)
