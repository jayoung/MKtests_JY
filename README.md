# MKtests_JY

## Janet Young (started Nov 2020)

# Overview

This git repo shares my R functions to perform McDonald-Kreitman tests, which are explained [here](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test).

Can polarize changes if a outgroup(s) are provided.

Can filter out low frequency variants (this sometimes results in increased fixed changes, which is counter-intuitive, but occurs when there are polymorphisms where the ancestral allele is at low frequency and the derived allele is the only one that remains after filtering).

Can specify a smaller region of the alignment to look at using the regionStartAA / regionEndAA options.


# Notes on the McDonald-Kreitman (MK) test

The MK test test was introduced in a [1991 paper by McDonald and Kreitman](https://www.nature.com/articles/351652a0), and there's a nice description on [Wikipedia](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test).

As well as a p-value testing for departure from neutrality, the MK test calculates: 
- *alpha*, which represents the proportion of substitutions driven by positive selection.  
- the *neutrality index*, where >1 indicates negative selection; <1 indicates positive selection

Reminder - best practise for the MK test is to: 
- use sequences from only a single population (e.g. ZI population of D. melanogaster) 
- to remove rare polymorphisms (e.g. remove alleles with frequency <= 5%)

Some people choose a chi-squared test to test for departures from neutrality (e.g. the MKT website), other people prefer a Fishers exact test ("FET").  My script gives both p-values, so the user can choose which to report.


# Instructions

In R, we first set up for analysis by reading in my functions: `source("scripts/MKfunctions.R")` 

Each input file should be a single fasta-format, in-frame, multiple sequence alignment containing sequences from two species, and population data from at least one of those species.  Optionally, it might also contain one or more outgroup sequences.

We run MK tests using the `doMKtest` function.   A basic example is shown here
```
source("scripts/MKfunctions.R")

alnFile <- "test_data/MKTwebsiteExample/MKTwebsite_testAln/MKTwebsite_testAln.fa"

MKresults_websiteExample <- doMKtest(
          myAlnFile=alnFile, 
          pop1seqs=c("pongo1", "pongo2"), 
          pop2seqs=c("trachy1", "trachy2") )
```

This [test script](scripts/example_script_MKtest.md), gives additional code examples, including how to polarize, remove rare variants, and various other options and tools.

# An alternative: the [MKT website](http://mkt.uab.es/mkt/MKT.asp)

I am using results provided by http://mkt.uab.es/mkt/MKT.asp to check my output. As of Jan 2025 I'm not sure that website is working.  The [iMKT site](https://imkt.uab.cat/index.html) might be intended as a replacement, but I haven't explored it.

A note about results from the MKT website: it ignores alignment positions where >=1 sequence contains a gap OR an N (e.g. demonstrate that using the two test alignments in MKTwebsite_testAln_addNseqs).  In contrast, my script includes those positions but does not count the gap or N as a change.  Only in the case where one of the populations has only gap or N at a position, then I cannot count any fixed changes in that position.

For particularly gappy alignments, this can be a problem. This is why we use the utility script `removeSeqsContainingNs.pl` (see [below](https://github.com/jayoung/MKtests_JY/tree/main?tab=readme-ov-file#utility-scripts-in-perl-to-help-use-popfly-data)) to remove any seq containing Ns if we want to perform MK tests using the website.  As an extreme example, imagine you have an alignment with quite some seqs with Ns, and the Ns are spread around the sequences, something like this:
```
    NNNNNNACGTAGCTA
    ACGTNNNNNNNNNNN
    ACGTAGACGTAGCTA
    ACGTAGACGTAGCTA
    ACGTAGACGTAGCTA
    ACGTAGACGTAGCTA
```

In this case the website ignores every single position in the alignment because of those first two gappy sequences, and there is nothing left to analyze.  


# Test datasets (in test_data)

## 1. website example
MKTwebsiteExample/MKTwebsite_testAln.fa
is the example provided on http://mkt.uab.es/mkt/MKT.asp 

also versions with stretches of Ns (MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addNseqs.fa)

and gaps (---) (MKTwebsiteExample/MKTwebsite_testAln_addNseqs/MKTwebsite_testAln_addGapSeqs.fa)


## 2a. test1.fa to test filtering polymorphisms based on frequency

- file: test_rareVariants/test1
- contains 10 seqs from pop1, 1 seq from pop2
- those seqs include the following nucleotide changes: 
    - Dn=1, Ds=1, Pn=1 (derived allele at 30%), Ps=1 (derived allele at 10%)
- If we use a 20% threshold to test filtering, the results should be:
    - Dn=1, Ds=1, Pn=1, Ps=1   without frequency filter
    - Dn=1, Ds=1, Pn=1, Ps=0   with 20% frequency filter


## 2b. test2.fa to test filtering polymorphisms based on frequency - more complex

- file: test_rareVariants/test2
- contains 10 seqs from pop1, 1 seq from pop2
- those seqs include the following nucleotide changes: 
    - Dn=1, Ds=1, Pn=2 (derived alleles both at 30%), Ps=1 (two derived alleles, at 90% and 10%)
    - The rare derived allele is ignored (Ps is reduced by 1), and the rare ancestral allele gets counted as a fixed change (Ps is reduced by 1, and Ds is increased by 1)
- If we use a 20% threshold to test filtering. Results should be:
    - Dn=1, Ds=1,  Pn=2,  Ps=2   without frequency filter
    - Dn=1, Ds=2,  Pn=2,  Ps=0   with 20% frequency filter


## 3. test3.fa very simple alignment that has only a single change, one rare synon SNP

test3.fa
pop1_seq01 has a synon change at 10% (bp 6)
should be 0 changes with a 20% freq filter

test3a.fa - same, but the rare change is in the second seq rather than the first
should be 0 changes with a 20% freq filter

test3b.fa - a single synonymous change at 90%
should be 1 fixed synonymous change with a 20% freq filter


## 4. test_unusualPolymorphism.fa

this test case inspired by a problem I had running code on one of Ching-Ho's alignments:
```
pop1_seq1   GCG
pop1_seq2   TCG
pop1_seq3   GAG
pop2_seq1   GAG
```

(file = BigH1_CDS_with_sim.fasta.  problem codon = bases 853-855, and choosing seqs 1,16,28,1115 demonstrates the issue)

In this case, the ancestral codon of pop1 is most likely to be GAG.

The A->C polymorphism at position 2 is counted OK as a non-synonymous polymorphism

The G->T polymorphism at position 1 created a problem for one version of my code - if it occurred on the background of the ancestral codon it would create a stop codon (GAG->TAG), which I don't know how to count So it must have occurred on the 'intermediate' codon (GCG->TCG), and therefore is also a non-synonymous codon.

I have fixed the code so that it counts this as 2 non-synonymous polymorphisms, just like the MKT website. 

I count changes from before codon -> after codon.  I use the ancestral codon but also any codons created by other SNPs in the same codon (ignoring any stop codons, as they should be impossible).


# Utility scripts (in perl) to help use Popfly data

## script 1 - takes a fasta file, removes any sequence containing one or more Ns

Useful, because the MK website will ignore changes in alignment positions where some sequences have Ns

usage: 
`removeSeqsContainingNs.pl seqfile(s).fasta`

what it does:
- removes all sequences containing 1 or more N bases. 
- output file names will end in .noNseqs.fa
- also creates report file called removeNseq_results.txt to show how many seqs were removed/retained


## script 2 - takes a popfly fasta file, splits it by population (ignores small populations, by default)

Alternatively, you could select individual populations when you download sequences from Popfly.

usage: 
`splitSeqsByPopulation.pl seqfile(s).fasta`

what it does:
- parses Popfly sequence names to determine which population each strain is from. 
- For any population with enough sequences, writes an output file containing only sequences from that population.  Default number of sequences is 50, but you can change that if you like.
- also creates report file called populationCounts.txt to show how many seqs were in each population
- after it finishes, it writes out names of populations it has seen in any of the files. If it doesn't seem to be processing population names quite right, let me know.
    
advanced usage, specifying an alternative name for the report file, and an alternative minimum population size (you can specify just one of those, don't need to specify both):

`splitSeqsByPopulation.pl --report=JYreport.txt --minPop=80 seqfile(s).fasta`


# To do

## definitely

make this an R package, with good documentation

finish methods description below

## maybe

Add alternative input methods eg separate input files for each population (and outgroup) - that would allow command-line scripting more easily.  

Continue script to get D simulans strain sequences for a gene from the vcf file and the reference assembly.  Don't think that's needed any more. Ching-Ho has a way to get ~34 Dsim seqs.

Figure out a more efficient way to extract D. mel alignments from the Popfly server

# Methods

(still a rough outline right now - just collecting notes as I see them)

Assumptions:  
- stop codons cannot occur.

## counting non-synonymous and synonymous changes for any pair of codons

There are 61 possible non-stop codons. We look up NS and S for a codon's evolutionary path (beforeCodon -> afterCodon) using the `codon_path` object from Bioperl's [`Bio::MolEvol::CodonModel` module](https://metacpan.org/release/CJFIELDS/BioPerl-1.007000_005/view/Bio/MolEvol/CodonModel.pm).  That object gives the NS and S counts for every one of the 3721 (=61*61) combinations of those 61 non-stop codons. 
The Bioperl module references "code and work from Alisha Holloway at UC Davis and Corbin Jones at UNC-Chapel Hill. See [Population Genomics: Whole-Genome Analysis of Polymorphism and Divergence in Drosophila simulans] (http://dx.doi.org/http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050310) and (http://www.dpgp.org/)"

That PLoS Biol article states: "For cases in which two alternative codons differed at more than one position, we used the pathway between codons that minimized the number of nonsynonymous substitutions. This is conservative with respect to the alternative hypothesis of adaptive evolution."

I copied `CodonModel.pm` locally (it's in `scripts/codonPaths`), parsed it to obtain a simple data frame using `scripts/codonPaths/getCodonPaths.R` and saved that data.frame as an R object (`scripts/codonPaths/codonPathsFromBioperl.Rdata`).

## dealing with ambiguity: there's not always one 'correct' answer

Sometimes (especially when we polarize changes) there's ambiguity in one or both of the 'beforeCodon' or 'afterCodon' sequences.  Example (where pop1=Dmel, pop2=Dsim and Dyak is the outgroup): 

```
>Dmel_1
ACA
>Dmel_2
ACA
>Dsim_1
ACC
>Dsim_2
ACC
>Dyak
ACT
```
Dmel-Dsim have one synonymous change.   However, in this case, the outgroup has a third allele: it does NOT help us polarize.  So the Dmel-Dsim ancestor could have been ACA or ACC, and the synonymous change could have happened along the Dmel OR the Dsim branch.

We can deal with sites like this with two possible strategies. Let's call them `mean` and `conservative` (and you can choose which to use in the `doMKtest` function using the `combiningApproach` option (the default behavior is `combiningApproach="conservative"`).

or each branch we're looking at, we first tabulate all the possible beforeCodon-afterCodon combinations. 

If we're using the `mean` approach, we simply take the mean of the `ns` counts over all possible combinations, and the mean of the `s` counts.

If we're using the `conservative` approach, we choose the beforeCodon-afterCodon combination that would give the most conservative result.  We first look at combination(s) with the minimum `ns` count, and if there's still >1 combination, we look at the one with the minimum `s` count.  Now we have a unique `ns` and `s` counts. 

(before Dec 7th 2022, my code had an unintended behavior. For situations where >1 codon combination had the minimum `ns` count, I was arbitrarily choosing the counts for the first of those codon combinations. This led to inconsistent results if we flipped the identities of pop1 and pop2, because the codon combinations appeared in a different order. After this change, results SHOULD be symmetrical no matter which population is pop1 and which is pop2.) 

Either way, we do NOT totally ignore codons with ambiguity.


# Some coding details

Dependencies: 
- `Biostrings` package

And for complete functionality:
- `openxlsx` package

The way I have things set up, it might work better if you're using R/Rstudio in "project" mode.  Rprojects are great, but if you're resistant to that for some reason, and you see weird errors about not being able to read files (e.g. `Warning: cannot open file`, `No such file or directory` ), you should edit the paths given in the top of the `scripts/MKfunctions.R` file. 
