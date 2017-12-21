args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("The gtf file should be provided as the first argument.n", call.=FALSE)
}

gtffile = args[1]
prefix = tools::file_path_sans_ext(gtffile)

library("GenomicFeatures")
library("spliceR")
library("BSgenome.Hsapiens.UCSC.hg38",character.only = TRUE)

ucscCDS <- getCDS(selectedGenome="hg38", repoName="UCSC")

sample.TxDb = makeTxDbFromGFF(gtffile, organism = "Homo sapiens")

chroms.to.keep = paste0("chr", 1:22)


sample.transcripts = GenomicFeatures::transcripts(sample.TxDb)
sample.transcripts = sample.transcripts[seqnames(sample.transcripts) %in% chroms.to.keep]
sample.transcripts$tx_id = as.character(sample.transcripts$tx_id)
seqlevels(sample.transcripts) = chroms.to.keep
sample.transcripts$"spliceR.isoform_id" = sample.transcripts$tx_id

exon.transcripts = exonsBy(sample.TxDb, by = "tx")
sample.exons = do.call(c, lapply(names(exon.transcripts), function(x) {b = exon.transcripts[[x]]; mcols(b) = cbind(mcols(b), tx_id = x); return(b) } ))

# sample.exons = GenomicFeatures::exons(sample.TxDb)
sample.exons = sample.exons[seqnames(sample.exons) %in% chroms.to.keep
seqlevels(sample.exons) = chroms.to.keep



sample.SpliceRList = SpliceRList(sample.transcripts, sample.exons, "hg38", "cufflinks","sample")

sample.PTC = annotatePTC(sample.SpliceRList, cds=ucscCDS, Hsapiens, PTCDistance=50)

save.image(sample.PTC, file = paste0(prefix, ".RData"))
