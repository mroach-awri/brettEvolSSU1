use strict;
use warnings;


my $usage = "
perl  reads.fastq.gz  bin1.id.gz  bin2.id.gz etc... \
    | gzip - > newfile.fq.gz
";

my $fq = shift or die $usage;
my @bins = @ARGV or die $usage;

my %get;

for my $bin (@bins){
    open my $TMP, '-|', "gunzip -c $bin" or die;
    while(<$TMP>){
        chomp;
        $get{$_}=1;
    }
    close$TMP or die;
}

open my $FQ, '-|', "gunzip -c $fq" or die;

while(1){
    my $h = <$FQ> or last;
    chomp($h);
    my $s = <$FQ>;
    my $b = <$FQ>;
    my $q = <$FQ>;
    
    $h =~ s/\s.*//;
    $h =~ s/^@//;
    if($get{$h}){
        print '@',$h,"\n",$s,$b,$q;
    }
}

close $FQ or die;




