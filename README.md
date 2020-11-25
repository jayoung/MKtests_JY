# MKtests_JY

## Janet Young Nov 2020

# Overview
My R functions to perform McDonald-Kreitman tests

Can polarize changes if an outgroup is provided

Can filter out low frequency variants (this sometimes results in increased fixed changes, which is counter-intuitive, but occurs when there are polymorphisms where the ancestral allele is at low frequency and the derived allele is the only one that remains after filtering)

Reminder:  best practise is to use sequences from only a single population (e.g. ZI population of D. melanogaster) and to remove rare polymorphisms (e.g. remove alleles with frequency <= 5%)

See https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test

## An alternative: the MKT website

I am using results provided by http://mkt.uab.es/mkt/MKT.asp to check my output. 

A note about results from that website: it ignores alignment positions where >=1 sequence contains a gap OR an N (e.g. demonstrate that using the two test alignments in MKTwebsite_testAln_addNseqs).  In contrast, my script includes those positions but does not count the gap or N as a change.  Only in the case where one of the populations has only gap or N at a position, then I cannot count any fixed changes in that position.

For particularly gappy alignments, this can be a problem. This is why we use Lisa's script to remove any seq containing Ns if we want to perform MK tests using the website.  As an extreme example, imagine you have an alignment with quite some seqs with Ns, and the Ns are spread around the sequences, something like this:
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

## 2. test1.fa to test filtering polymorphisms based on frequency

file: test_rareVariants/test1

contains 10 seqs from pop1, 1 seq from pop2

contains these changes: 

Dn=1, Ds=1, Pn=1 (derived allele at 30%), Ps=1 (derived allele at 10%)

Use a 20% threshold to test filtering. Results should be:

Dn=1, Ds=1, Pn=1, Ps=1   without frequency filter

Dn=1, Ds=1, Pn=1, Ps=0   with 20% frequency filter


## 2. test2.fa to test filtering polymorphisms based on frequency - complex

file: test_rareVariants/test2

contains 10 seqs from pop1, 1 seq from pop2

contains these changes: 

Dn=1, Ds=1, Pn=2 (derived alleles both at 30%), Ps=1 (derived alleles at 90% and 10%)

The rare derived allele is ignored (Ps reduced by 1), and the rare ancestral allele gets counted as a fixed change (Ps reduced by 1, and Ds increased by 1)

Use a 20% threshold to test filtering. Results should be:

Dn=1, Ds=1,  Pn=2,  Ps=2   without frequency filter

Dn=1, Ds=2,  Pn=2,  Ps=0   with 20% frequency filter


## 3. test3.fa very simple alignment that has only a single change, one rare synon SNP

test3.fa
pop1_seq01 has a synon change at 10% (bp 6)
should be 0 changes with a 20% freq filter

test3a.fa - same, but the rare change is in the second seq rather than the first
should be 0 changes with a 20% freq filter

test3b.fa - a single synonymous change at 90%
should be 1 fixed synonmouse change with a 20% freq filter


# Utility scripts (perl) to help use Popfly data

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

## maybe

To help troubleshooting, make a trim function, to consider only some portions of the alignment (and it might be useful in general)

Add alternative input methods eg separate input files for each population (and outgroup) - that would allow a script

Add methods tab to spreadsheet to capture parameters and warnings, perhaps also a tab to name the seqs used
