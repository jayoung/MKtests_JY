#!/usr/bin/perl

use warnings;
use strict;

#### JY script to remove any sequences that contain 1 or more N bases (this version does NOT depend on Bioperl). This version from Nov 13, 2020

# we'll summarize what we did to each file here:
my $reportfile = "removeNseq_results.txt";

## split input file using this delimiter:
$/ = "\n>";


open (REPORT, "> $reportfile");
print REPORT "File\tNum seqs\tNum seqs retained\tNum seqs removed\n";
foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - file $file does not exist\n\n";
    }
    print "\n## processing file $file\n";
    
    ## figure out name of output file, and open it
    my $outfile = $file; $outfile =~ s/\.fasta$//; $outfile =~ s/\.fa$//;
    $outfile .= ".noNseqs.fa";
    # make sure outfiles go in current directory, even if input files are elsewhere:
    if ($outfile =~ m/\//) { $outfile = (split /\//, $outfile)[-1]; }
    open (OUT, "> $outfile");
    
    ## read in all seqs
    open (IN, "< $file");
    my @seqs = <IN>;
    close IN;
    # first seq doesn't look quite like the rest
    $seqs[0] =~ s/^>//;
    
    ## process the seqs
    my $seqCount = 0;
    my $seqCountKept = 0;
    my $seqCountRemoved = 0;
    foreach my $seq (@seqs) {
        # the file split delimiter stays at the end of each chunk, so fix that:
        $seq =~ s/>$//; $seq = ">$seq";
        $seqCount++;
        
        ## for troubleshooting:
        #if ($seqCount > 3) { die; } 
        
        ## split the sequence up into header row (seq name plus description, if any), and the rest
        my $header = $seq; $header =~ m/^(>.+?\n)/; $header = $1; $header =~ s/\n$//;
        
        my $seqWithoutHeader = $seq;
        # the \Q and \E in the next line are there to handle any special characters found in $header, e.g. ( or )
        $seqWithoutHeader =~ s/\Q$header\E\n//;
        $seqWithoutHeader =~ s/\n//g; # remove line breaks in seq and at end
        
        # split header into seqname and description (if any)
        # not actually using seqname and description for now, but if I want to write a report on which seqs were removed it might be useful
        #my $seqname = $header; $seqname =~ s/^>//; 
        #my $description = "";
        #if ($seqname =~ m/\s/) { 
        #    $seqname =~ m/(.+?)\s/;  $seqname = $1; 
        #    $description = $header; $description =~ s/^>//;
        #    $description =~ s/\Q$seqname\E\s//;
        #}
        
        ## test whether seq contains any Ns
        if ($seqWithoutHeader =~ m/N/i) { $seqCountRemoved++; next; }
        ## print sequence to output file (with original line wrapping, in this case)
        print OUT $seq;
        $seqCountKept++;
    }
    close OUT;
    
    ## report what happened
    print REPORT "$file\t$seqCount\t$seqCountKept\t$seqCountRemoved\n";
    print "    $seqCount seqs processed\n";
    print "        $seqCountKept N-less seqs retained\n";
    print "        $seqCountRemoved seqs containing N(s) removed\n";
    if ($seqCountKept == 0) {
        print "\n    WARNING!!! there were no seqs left after removing any that contain Ns.  There will be no output file\n\n";
        unlink $outfile;
    }
}
close REPORT;
