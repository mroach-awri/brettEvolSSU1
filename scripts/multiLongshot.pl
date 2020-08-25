use strict;
use warnings;
use threads;
use Thread::Semaphore;
use Thread::Queue;


my $usage = "
perl multiLongshot.pl  fasta.fai  aligned.bam  threads
";

my $fai = shift or die $usage;
my $fasta = $fai;
$fasta =~ s/\.fai$//;
my $bam = shift or die $usage;
my $threads = shift or die $usage;


my $queue = Thread::Queue->new();
my $available_threads = Thread::Semaphore->new($threads);


open my $FAI, '<', $fai or die;
while(<$FAI>){
    my@l=split/\s+/;
    $queue->enqueue($l[0]);
}
close $FAI;

$queue->end();


for (1..$threads){
    $available_threads->down(1);
    threads->create(\&longshotjob);
}

# wait for workers to finish
$available_threads->down($threads);
$available_threads->up($threads);

# join workers
for my $thr (threads->list()){
    $thr->join();
}

print STDERR "done!\n";

exit(0);



sub longshotjob {
    while (defined(my $ctg = $queue->dequeue())) {
        system("longshot -f $fasta -b $bam -o $ctg.vcf -r $ctg -p $ctg") == 0 or die;
    }
    # exit thread
    $available_threads->up(1);

    return;
}
    








