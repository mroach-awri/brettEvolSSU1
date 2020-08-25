#!/usr/bin/env perl

# Copyright (c) 2020 Michael Roach (Australian Wine Research Institute)

use strict;
use warnings;
use Data::Dumper;


# script to take a fasta .fai and bed annotations and generate table files 
# of ggplot2 geom_path and geom_area for generating a chromosome ideogram

my $usage = "Usage:
bed2geom_ideo.pl  <fasta.fai>  <annot1.bed>  [  <annot2.bed>  ...  ]

ideograms written out in their input order
";


# parse args
my $fai = shift or die $usage;
(@ARGV) or die $usage;
my @bed;
(push @bed, $_) for (@ARGV);

# variables for plotting
my $chr_padding = 1;
my $chr_ideo_width = 1;
my $tracks = $#bed;

# global vars
my @chr;
my @len;
my %chrx;

# read in the fai index
open my $FAI, '<', $fai or die "failed to open $fai for reading";
my$n=0;
while(<$FAI>){
  my@l=split /\s+/;
  push @chr, $l[0];
  push @len, $l[1];
  $chrx{$l[0]}=$n;
  $n++;
}
close $FAI or die 'failed to close file handle';

# write out the geom_path chromosome track outlines
  # open file for writing
open my $OL, '>', "$fai.geom_path" or die "failed to open $fai.geom_path for writing";
my $ciel = (@chr * $chr_padding) + (@chr * @bed * $chr_ideo_width);
my $bedheight = (@bed * $chr_ideo_width);
for my$c(0..$#chr){
  my $bedheight = 0;
  for my$b(0..$#bed){
    my $yoff = ($ciel - $bedheight);
    my $yofflow = $yoff-$chr_ideo_width;
    print $OL "0\t$yoff\n$len[$c]\t$yoff\n$len[$c]\t$yofflow\n0\t$yofflow\n0\t$yoff\nNA\tNA\n";
    $bedheight += $chr_ideo_width;
  }
  $ciel -= ($chr_padding + $bedheight);
}
close $OL or die 'failed to close file handle';

# iterate the bed annotation files
my$b=0;
for(@bed){
  open my $BED, '<', $_ or die "failed to open $_ for reading";
  open my $OA, '>', "$_.geom_area" or die "failed to open $_.geom_area for writing";
  while(<$BED>){
    my@l=split/\s+/;
    # we'll use the fai as a sort of filter and ignore entries not present in the fai
    next unless defined($chrx{$l[0]});
    my $yoh = (@chr * $chr_padding) + (@chr * @bed * $chr_ideo_width);
    $yoh -= ($chrx{$l[0]} * $chr_padding) + ($chrx{$l[0]} * @bed * $chr_ideo_width) + ($b * $chr_ideo_width);
    my $yol = $yoh - $chr_ideo_width;
    print $OA "$l[1]\t$yoh\n$l[2]\t$yoh\n$l[2]\t$yol\n$l[1]\t$yol\nNA\tNA\n";
  }
  close $OA or die 'failed to close file handle';
  $b++;
}














