plot.ptc = function(tx.id, SpliceRListObj, genomeObject) {

    exon.granges = SpliceRListObj[["exon_features"]] [SpliceRListObj[["exon_features"]]$"tx_id.X" == tx.id, ]
    
    transcript.granges = SpliceRListObj[["transcript_features"]] [SpliceRListObj[["transcript_features"]]$"tx_id" == tx.id, ]
    exonSeq = getSeq(genomeObject, exon.granges)

    totalSeq = tolower(paste(exonSeq, collapse = "" ))


                                        # Define a vector with the sequences of potential start and stop codons
    codons <- c("atg", "taa", "tag", "tga")
                                        # Find the number of occurrences of each type of potential start or stop codon
    for (i in 1:4)
    {
        codon <- codons[i]
                                        # Find all occurrences of codon "codon" in sequence "sequence"
        occurrences <- matchPattern(codon, totalSeq)

                                        # Find the start positions of all occurrences of "codon" in sequence "sequence"
        codonpositions <- attr(occurrences@ranges,"start")
                                        # Find the total number of potential start and stop codons in sequence "sequence"
        numoccurrences <- length(codonpositions)
        if (i == 1)
        {
                                        # Make a copy of vector "codonpositions" called "positions"
            positions   <- codonpositions
                                        # Make a vector "types" containing "numoccurrences" copies of "codon"
            types       <- rep(codon, numoccurrences)
        }
        else
        {
                                        # Add the vector "codonpositions" to the end of vector "positions":
            positions   <- append(positions, codonpositions, after=length(positions))
                                        # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
            types       <- append(types, rep(codon, numoccurrences), after=length(types))
        }
    }
                                        # Sort the vectors "positions" and "types" in order of position along the input sequence:
    indices <- order(positions)
    positions <- positions[indices]
    types <- types[indices]
                                        # Make a plot showing the positions of the start and stop codons in the input sequence:

                                        # Draw a line at y=0 from 1 to the length of the sequence:

    num.of.exons = length(exon.granges)
    num.of.introns =  num.of.exons - 1 
    
    
    transcript.length  <- nchar(totalSeq) # length of sequence

    # setting intronic segments to be quarter of exon lengths
    intron.seqlength =  transcript.length %/% 4  %/% num.of.introns + 1 

    total.length = transcript.length + intron.seqlength * num.of.introns

    y <- c(0,0)
    plot(c(1, total.length), y, ylim=c(0,3), type="n", axes=FALSE, xaxt = "n", 
         xlab="Nucleotide", ylab="Reading frame",
         main="Predicted start (red) and stop (blue) codons")
    segments(1,1,total.length,1)
    segments(1,2,total.length,2)
                                        # Add the x-axis at y=0:
    # axis(1, pos=0)
                                        # Add the y-axis labels:
    text(0,0.5,"+1", adj = 1)
    text(0,1.5,"+2", adj = 1)
    text(0,2.5,"+3", adj = 1)
                                        # Draw in each predicted start/stop codon:

#### plotting exons
    exon.meta = as.data.frame(exonSeq@ranges)
    b = exon.meta
    rect(b$start[1], 0, b$end[1], 3, border = "black", lwd = 1.4)
    axis(1, at = c(b$start[1], b$end[1]), labels = FALSE)
    text (x = c(b$start[1], b$end[1]), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
          labels = c(exon.meta$start[1], exon.meta$end[1]), srt = 90, xpd = TRUE) 
    for (i in 2:nrow(b)) {

        b$start[i] = exon.meta$start[i] + (i - 1) * intron.seqlength
        b$end[i] = exon.meta$end[i] + (i - 1) * intron.seqlength
        rect(b$start[i], 0, b$end[i], 3, border = "black", lwd = 1.4)
        segments(b$end[i - 1], 1.5, b$end[i - 1] + intron.seqlength, col = "black", lwd = 1.4)

        axis(1, at = c(b$start[i], b$end[i]), labels = FALSE)
        text (x = c(b$start[i], b$end[i]), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
              labels = c(exon.meta$start[i], exon.meta$end[i]), srt = 90, xpd = TRUE)
    }

    
#### updating exon positions

    update.pos = function(pos, exon.metadata) {
        exon.ends = c(0, cumsum(exon.metadata$width))
        exon.index = which(pos %/%exon.ends == 0 )[1] - 1  ## gets in which exon the variant is
        pos.updated = pos - exon.ends[exon.index] + exon.metadata$start[exon.index]
        pos.updated
    }

    
    numcodons <- length(positions)
    for (i in 1:numcodons)
    {
        position <- update.pos(positions[i], b)
        type <- types[i]
        remainder <- (position-1) %% 3
        if    (remainder == 0) # +1 reading frame
        {
            if (type == "atg") { segments(position,0,position,1,lwd=1,col="red") }
            else               { segments(position,0,position,1,lwd=1,col="blue")}
        }
        else if (remainder == 1)
        {
            if (type == "atg") { segments(position,1,position,2,lwd=1,col="red") }
            else               { segments(position,1,position,2,lwd=1,col="blue")}
        }
        else if (remainder == 2)
        {
            if (type == "atg") { segments(position,2,position,3,lwd=1,col="red") }
            else               { segments(position,2,position,3,lwd=1,col="blue")}
        }
    }    
}


