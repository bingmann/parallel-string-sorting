#!/usr/bin/perl -w
#
# Simple script to extract only the DNA sequence from the human genome
# published on Project Gutenberg.
#
# Usage:
#
# 1) download documents 11775.txt to 11799.txt from Project Gutenberg. These
# are titled "Human Genome Project, Build 34", Chromosome Number 1 to 22, X, Y
# and Supplementary Material.
#
# 2) run "perl gutenberg-hg34-split.pl 117*.txt"
# to extract all chrXX.txt files, concatenate them into dna.txt in the order
# they appear in the documents, and output junk.txt which contains license text
# added by Project Gutenberg.
#
# 2013-03-28 Timo Bingmann <tb@panthema.net>
#

use strict;
use warnings;

my $chrfile;

open(DNA,"> dna.txt");
open(JUNK,"> junk.txt");

while ( my $ln = <> )
{
    $ln =~ s/\r\n$//;

    if ($ln =~ /^>(.+)$/)
    {
        if ($chrfile) {
            print("closed chr $chrfile\n");
            close(CHR);
            $chrfile = undef;
        }

        $chrfile = $1;
        print "new chr file $chrfile\n";
        open(CHR, "> $chrfile.txt") or die;
    }
    elsif ($ln =~ /^[AGCTNagctnMR]+$/) {
        print(CHR $ln) or die;
        print(DNA $ln) or die;
    }
    else {
        if ($chrfile) {
            print("closed chr $chrfile\n");
            close(CHR);
            $chrfile = undef;
        }
        print(JUNK "$ln\n");
    }
}
