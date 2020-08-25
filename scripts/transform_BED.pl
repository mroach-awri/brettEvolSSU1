#!/usr/bin/env perl

# Copyright (c) 2020 Michael Roach (Australian Wine Research Institute)

use strict;
use warnings;

# script to translate the contig bed coords to the scaffold coords using the translation BED

my $usage = "
perl transform_BED.pl translation.bed annot.bed

";

open my $TRN, shift or die $usage;
open my $BED, shift or die $usage;

my %trans;

while(<$TRN>){
    chomp($_);
    my@l=split/\s+/;
    $trans{$l[3]}{"$l[4]:$l[5]"}{C} = $l[0];                # chr
    $trans{$l[3]}{"$l[4]:$l[5]"}{O} = $l[1];                # chr start
    @{$trans{$l[3]}{"$l[4]:$l[5]"}{S}} = ($l[4], $l[5]);    # ctg span
    $trans{$l[3]}{"$l[4]:$l[5]"}{D} = $l[6];                # direction - not used
    $trans{$l[3]}{"$l[4]:$l[5]"}{L} = $l[7];                # ctg length
}
close $TRN;

while(<$BED>){
    chomp($_);
    my @l=split(/\t/,$_);
    next if (!($trans{$l[0]}));
    my $outPrint;
    # find which section it belongs to
    T: for my $k (keys %{$trans{$l[0]}}){
        if (($l[1] >= @{$trans{$l[0]}{$k}{S}}[0]) and ($l[2] <= @{$trans{$l[0]}{$k}{S}}[1])){
            # subtract the range offset
            $l[1] -= @{$trans{$l[0]}{$k}{S}}[0];
            $l[2] -= @{$trans{$l[0]}{$k}{S}}[0];
            # transform the coords
            $l[1] += $trans{$l[0]}{$k}{O};
            $l[2] += $trans{$l[0]}{$k}{O};
            # set the new chr name
            $l[0] = $trans{$l[0]}{$k}{C};
            $outPrint=1;
            last T;
        }
    }
    if ($outPrint){
        my $o;
        ($o .= "$_\t") for (@l);
        $o =~ s/\t$/\n/;
        print STDOUT $o;
    }
}
