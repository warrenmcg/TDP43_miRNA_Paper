#!/usr/bin/Rscript

suppressMessages(library(GO.db))
suppressMessages(library(KEGGREST))

message("Retrieving GO terms")
terms <- Term(GOTERM)
ontologies <- Ontology(GOTERM)
go_terms <- data.frame(category=ontologies, id=names(terms),name=as.character(terms),stringsAsFactors=F)
go_terms$category <- paste0("GO",go_terms$category)

# The following pathways were missing from the GO.db because they're older/obsolete terms used by FatiGO
go_terms <- rbind(go_terms, c("GOBP", "GO:0034984", "cellular response to DNA damage stimulus"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0006916", "anti-apoptosis"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0007243", "intracellular protein kinase cascade"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0045941", "positive regulation of transcription"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0016044", "cellular membrane organization"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0016568", "chromatin modification"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0006944", "membrane fusion"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0008633", "activation of pro-apoptotic gene products"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0042375", "quinone cofactor metabolic process"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0006917", "induction of apoptosis"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0003078", "regulation of natriuresis")
go_terms <- rbind(go_terms, c("GOBP", "GO:0006467", "protein thiol-disulfide exchange"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0010149", "senescence"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0019047", "provirus integration"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0019059", "initiation of viral infection"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0030146", "diuresis"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0031557", "induction of programmed cell death in response to chemical stimulus"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0031575", "G1/S transition checkpoint"))
go_terms <- rbind(go_terms, c("GOBP", "GO:0035121", "tail morphogenesis"))

message("Retrieving KEGG Pathway terms (this will take a few minutes)")
kegg_ids <- as.character(read.table("pathway_input_lists/kegg_pathway_list.txt", stringsAsFactors=F)[,1])
num_ids <- length(kegg_ids)
kegg_list <- sapply(seq(0, num_ids, 10), function(x) lapply(keggGet(kegg_ids[(x+1):(x+10)]), "[[", "NAME"))

kegg_names <- unlist(kegg_list)
# The following pathways were missing from the KEGGREST database (older KEGG pathways?)
kegg_names <- append(kegg_names, "Methane metabolism - Homo sapiens (human)", after=189)
kegg_names <- append(kegg_names, "Cyanoamino acid metabolism - Homo sapiens (human)", after=183)
kegg_names <- append(kegg_names, "Limonene and pinene degradation - Homo sapiens (human)", after=178)

# Use grep to sub out the redundant human part of the pathway names
kegg_names <- gsub(" - Homo sapiens (human)","",kegg_names, fixed=T)

kegg_terms <- data.frame(category="KEGG",id=kegg_ids,name=kegg_names,stringsAsFactors=F)

message("Retrieving Biocarta Terms")
biocarta <- read.table("pathway_input_lists/biocarta_pathway_list.txt",sep="\t", stringsAsFactors=F, quote="")
names(biocarta) <- c("id","name")
biocarta <- data.frame(category="BIOCARTA",biocarta,stringsAsFactors=F)

message("Retrieving Reactome terms")
reactome <- read.table("pathway_input_lists/reactome_pathway_list.txt",sep="\t", stringsAsFactors=F, quote="")
names(reactome) <- c("id","name")
reactome <- data.frame(category="REACTOME",reactome,stringsAsFactors=F)

message("Combining the pathways and writing the resulting table to file")
full_table <- rbind(go_terms, kegg_terms, biocarta, reactome)
write.table(full_table, file="pathway_ids2names.txt", sep="\t", row.names=F, quote=F)