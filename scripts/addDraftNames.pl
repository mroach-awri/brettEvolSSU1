use strict;
use warnings;

my $usage = "
perl addDraftNames.pl  outfmt6.gz  nameTable.tsv  proteins.fasta

";


my $blast = shift or die $usage;
my $table = shift or die $usage;
my $prot = shift or die $usage;

my %nom;
my %noom;

open my $BLS, "zcat $blast |" or die;
while(<$BLS>){
    my@l=split/\s+/;
    $nom{$l[0]}=$l[1];
}
close $BLS;

open my $TBL, '<', $table or die;
while(<$TBL>){
    chomp;
    my@l=split/\t/;
    $noom{$l[0]} = $l[1];
}
close $TBL;

open my $PRO, '<', $prot or die;
while(<$PRO>){
    if($_=~/^>(\S+)/){
        print STDOUT ">$1 ", ($nom{$1}) ? $noom{$nom{$1}} : "hypothetical protein", "\n";
    } else {
        print STDOUT $_;
    }
}
close $PRO;

exit(0);
