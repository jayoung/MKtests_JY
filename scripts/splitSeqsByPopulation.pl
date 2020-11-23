#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

#### JY script to go through fasta file(s) downloaded from popfly and to create new files where sequences are sorted according to which population they come from 
# This version does NOT depend on Bioperl 
# This version from Nov 16, 2020 (no functional changes since Nov 13 version)


###### default options:
# has an option to only write seqs for populations with at least a certain number of members)
my $minimumPopSizeToWriteFasta = 50;

# we'll summarize what we did to each file here:
my $reportfile = "populationCounts.txt";

###### get any command line options

## get any non-default options from commandline
GetOptions("report=s" => \$reportfile,  ## default spcified above
           "minPopSize=i" => \$minimumPopSizeToWriteFasta
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";

print "\noptions:\n    reportfile $reportfile\n    minimumPopSizeToWriteFasta $minimumPopSizeToWriteFasta\n\n";

#############

## split input file using this delimiter:
$/ = "\n>";


open (REPORT, "> $reportfile");
print REPORT "File\tNum seqs full file\tPopulation\tNum seqs\n";
my %allPopulations; # just for troubleshooting - we will report the names of all populations we encountered after going through all files, to make sure we processed seq names right
foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - file $file does not exist\n\n";
    }
    print "\n## processing file $file\n";
    
    ## read in all seqs
    open (IN, "< $file");
    my @seqs = <IN>;
    close IN;
    # first seq doesn't look quite like the rest
    $seqs[0] =~ s/^>//;
    
    ## process the seqs
    my $seqCount = 0;
    my %populationCounts;
    my %populationSeqs; ## hash of arrays
    foreach my $seq (@seqs) {
        # the file split delimiter stays at the end of each chunk, so fix that:
        $seq =~ s/>$//; $seq = ">$seq";
        $seqCount++;
        
        ## split the sequence up into header row (seq name plus description, if any), and the rest
        my $header = $seq; $header =~ m/^(>.+?\n)/; $header = $1; $header =~ s/\n$//;
        
        my $seqWithoutHeader = $seq;
        # the \Q and \E in the next line are there to handle any special characters found in $header, e.g. ( or )
        $seqWithoutHeader =~ s/\Q$header\E\n//;
        $seqWithoutHeader =~ s/\n//g; # remove line breaks in seq and at end
        
        ## split header into seqname and description (if any)
        my $seqname = $header; $seqname =~ s/^>//; 
        if ($seqname =~ m/\s/) { 
            $seqname =~ m/(.+?)\s/;  $seqname = $1; 
        }
        
        ## determine population
        my $pop = $seqname;
        #my $pop = "ED19-1-3_ChrX_(reversed)"; ## test
        #print "    pop before $pop\n";
        $pop =~ s/_Chr.+?$//;
        $pop =~ s/-.*$//; # the * is a greedy pattern match - maximizes length of matched string (whereas the ? after Chr above minimizes, not that it matters for Chr which I choose. Here we are greedy to make sure we match the first hyphen
        $pop =~ s/[ABN]$//;
        $pop =~ s/\d+?$//;
        #$pop =~ s/-$//;
        $pop =~ s/_$//;
        #print "    pop after $pop\n";
        #die;
        $populationCounts{$pop}++;
        push @{$populationSeqs{$pop}}, $seq;
    }
    print "    $seqCount seqs processed\n";
    
    foreach my $population (sort keys %populationSeqs) {
        $allPopulations{$population}=1;
        my @seqsThisPopulation = @{$populationSeqs{$population}};
        my $numSeqsThisPopulation = @seqsThisPopulation;
        print REPORT "$file\t$seqCount\t$population\t$numSeqsThisPopulation\n";
        print "        $numSeqsThisPopulation seqs in population $population\n";
        ## skip small populations
        if ($numSeqsThisPopulation < $minimumPopSizeToWriteFasta) { next; }
        
        ## figure out name of output fasta file, and open it
        my $outfile = $file; $outfile =~ s/\.fasta$//; $outfile =~ s/\.fa$//;
        $outfile .= ".$population.fa";
        # make sure outfiles go in current directory, even if input files are elsewhere:
        if ($outfile =~ m/\//) { $outfile = (split /\//, $outfile)[-1]; }
        open (OUT, "> $outfile");
        foreach my $seq (@seqsThisPopulation) { print OUT $seq; }
        close OUT;
    }
}
close REPORT;
print "\n#### Done. Encountered these populations in any of the files:\n";
foreach my $population (sort keys %allPopulations) { print "    $population\n"; }
print "\n";
