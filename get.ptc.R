library("GenomicFeatures")
library("spliceR")

library("BSgenome.Hsapiens.UCSC.hg38",character.only = TRUE)

ucscCDS <- getCDS(selectedGenome="hg38", repoName="UCSC")

A549.gtffile = "A549.transcripts_stranded.gtf"
A549.TxDb = makeTxDbFromGFF(A549.gtffile, organism = "Homo sapiens")

chroms.to.keep = paste0("chr", 1:22)


A549.transcripts = GenomicFeatures::transcripts(A549.TxDb)
A549.transcripts = A549.transcripts[seqnames(A549.transcripts) %in% chroms.to.keep]
seqlevels(A549.transcripts) = chroms.to.keep
A549.transcripts$"spliceR.isoform_id" = A549.transcripts$tx_id

exon.transcripts = exonsBy(A549.TxDb, by = "tx")
A549.exons = do.call(c, lapply(names(exon.transcripts), function(x) {b = exon.transcripts[[x]]; mcols(b) = cbind(mcols(b), tx_id = x); return(b) } ))

# A549.exons = GenomicFeatures::exons(A549.TxDb)
A549.exons = A549.exons[seqnames(A549.exons) %in% chroms.to.keep]
seqlevels(A549.exons) = chroms.to.keep



A549.SpliceRList = SpliceRList(A549.transcripts, A549.exons, "hg38", "cufflinks","A549")

A549.PTC = annotatePTC(A549.SpliceRList, cds=ucscCDS, Hsapiens, PTCDistance=50)
