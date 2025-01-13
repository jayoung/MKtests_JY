#### parses BioPerl module for Bio::MolEvol::CodonModel to obtain the number of synonymous and non-synonymous mutations it takes to mutate from one codon to another (not including stop codons). Jason Stajich <jason-at-bioperl-dot-org.

## cp /fh/fast/malik_h/grp/malik_lab_shared/perl/Bio-MolEvol/lib/Bio/MolEvol/CodonModel.pm .

codonPathsFile <- "CodonModel.pm"

dat <- scan(codonPathsFile, what="character", sep="\n")

firstRow <- grep("sub codon_path", dat) + 2
lastRow <- length(dat) - 3

dat <- dat[firstRow:lastRow]

if (length(dat) != 61^2) {
    stop("ERROR - there should be 3712 (61^2) codon paths. but there are", length(dat), "\n" )
}

dat <- gsub("[ \\'\\[\\]]","",dat, perl=TRUE)
dat <- gsub("=>",",",dat)


dat <- strsplit(dat, ",")

paths <- sapply(dat, "[[", 1)

codonPaths <- data.frame( path=paths,
                   codon1=substr(paths, 1,3),
                   codon2=substr(paths, 4,6),
                   ns=as.integer(sapply(dat, "[[", 2)),
                   s=as.integer(sapply(dat, "[[", 3)), 
                   stringsAsFactors=FALSE )

save(codonPaths, file="codonPathsFromBioperl.Rdata")
