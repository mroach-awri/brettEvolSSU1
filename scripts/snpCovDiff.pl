#!/usr/bin/env perl

# Copyright (c) 2020 Michael Roach (Australian Wine Research Institute)

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my $usage = "
perl snpCovDiff.pl  -rd ref.seqCov  -qd query.seqCov  -rs ref.snp  -rc ref.cov  -qs evolv.snp  -qc evolv.cov > out.tsv
";

my $rSeq;   # ref total seq coverage
my $eSeq;   # query total seq coverage

my $rSnp;   # ref snp hist
my $rCov;   # ref depth hist
my $eSnp;   # query snp hist
my $eCov;   # query depth hist 

GetOptions (
    "rd=s" => \$rSeq,
    "qd=s" => \$eSeq,
    "rs=s" => \$rSnp,
    "rc=s" => \$rCov,
    "qs=s" => \$eSnp,
    "qc=s" => \$eCov
);

(-s $rSeq) and (-s $eSeq) and (-s $rSnp) and (-s $rCov) and (-s $eSnp) and (-s $eCov) or die $usage;


my $covDiff = int(`cat $rSeq`) / int(`cat $eSeq`);


my %snp;    # $snp{ctg}{start stop} = ref snps
my %cov;    # $cov{ctg}{start stop} = ref coverage
my %out;    # $out{ctg}{start stop}{s/c} = % change


# slurp the reference strain snp and coverage
slurpBed($rSnp, \%snp);
slurpBed($rCov, \%cov);

# parse the evolved strains and build the output
parseEvo($eSnp, \%snp, \%out, 's');
parseEvo($eCov, \%cov, \%out, 'c');

# sort and print the output
for my $c (sort keys %out){
    for my $r (sort keys %{$out{$c}}){
        print STDOUT "$c\t$r\t$out{$c}{$r}{s}\t$out{$c}{$r}{c}\n";
    }
}

# done
exit;




sub slurpBed {
    my $f = $_[0];
    my $h = $_[1];
    
    open my $TMP, '<', $f or die;
    
    while(<$TMP>){
        my@l=split/\s+/;
        $$h{$l[0]}{"$l[1]\t$l[2]"} = $l[3];
    }
    
    return;
}


sub parseEvo {
    my $f = $_[0];  # file
    my $r = $_[1];  # ref hash
    my $o = $_[2];  # out hash
    my $t = $_[3];  # tag ('s' or 'c')
    
    open my $TMP, '<', $f or die;
    
    while(<$TMP>){
        my@l=split/\s+/;
        
        # scale coverage
        if ($t eq 'c'){
            $l[3] *= $covDiff;
        }
        
        if ($$r{$l[0]}{"$l[1]\t$l[2]"} == 0){
            $$o{$l[0]}{"$l[1]\t$l[2]"}{$t} = 0;
        } else {
            $$o{$l[0]}{"$l[1]\t$l[2]"}{$t} = ( abs($l[3]) - abs($$r{$l[0]}{"$l[1]\t$l[2]"}) ) ;
        }
    }
    
    return;
}






