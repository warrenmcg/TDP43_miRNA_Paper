suppressMessages(library(HiveR))
suppressMessages(library(gridExtra))
suppressMessages(library(optparse))

options(error=traceback)

option.list <- list(
    make_option(c("-d","--downDir"), type="character", default=NULL, dest="down_dir", 
                metavar="directory", help="[REQUIRED] Input directory with the downGene HivePlot Dot Files."),
    make_option(c("-u","--upDir"), type="character", default=NULL, dest="up_dir", 
                metavar="directory", help="[REQUIRED] Input directory with the upGene HivePlot Dot Files."),
    make_option(c("-p","--prefix"), type="character", default="result", dest="out_prefix",
                metavar="prefix", help="Prefix for output files")
)

parser <- OptionParser(option_list=option.list)
opt <- parse_args(parser)

down_dot_file <- list.files(opt$down_dir,pattern=".dot", full.name=T)
down_node_file <- list.files(opt$down_dir,pattern="nodeInst.csv", full.name=T)
down_edge_file <- list.files(opt$down_dir,pattern="edgeInst.csv", full.name=T)

up_dot_file <- list.files(opt$up_dir,pattern=".dot", full.name=T)
up_node_file <- list.files(opt$up_dir,pattern="nodeInst.csv", full.name=T)
up_edge_file <- list.files(opt$up_dir,pattern="edgeInst.csv", full.name=T)

message("creating the hive plot objects (this may take several minutes)...")
downUp.hpd <- dot2HPD(file = down_dot_file, node.inst = down_node_file, 
					  edge.inst = down_edge_file, axis.cols = rep("black",3))
upDown.hpd <- dot2HPD(file = up_dot_file, node.inst = up_node_file, 
					  edge.inst = up_edge_file, axis.cols = rep("black",3))

message("creating the up-gene, down-mirna normalized hive plot")
pdf(paste0(opt$out_prefix, "_upDown_normHivePlot.pdf"))
plotHive(upDown.hpd, 0.1, "norm", bkgnd="white", axLabs=c("Gene", "miRNA", "Pathway / GO"), 
		 axLab.pos = c(0.1, 0.2, 0.2), axLab.gpar = grid::gpar(fontfamily="ArialMT"))
dev.off()

message("creating the down-gene, up-mirna normalized hive plot")
pdf(paste0(opt$out_prefix, "_downUp_normHivePlot.pdf"))
plotHive(downUp.hpd, 0.1, "norm", bkgnd="white", axLabs=c("Gene", "miRNA", "Pathway / GO"), 
		 axLab.pos = c(0.1, 0.15, 0.15), axLab.gpar = grid::gpar(fontfamily="ArialMT"))
dev.off()

message("Generating Hive Plots complete")
