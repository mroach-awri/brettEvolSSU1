

## Overview

 - Create haploid assembly of AWRI1499
 - Phase and bin reads into the A1 + A2 haplomes and B haplome
 - Separately assembled the A1 + A2 haplomes and the B haplome


## Initial assembly

### Canu

```bash
canu -p 1499 -d CANU-NP-PB-2 \
    -genomeSize=24m -stopOnReadQuality=false \
    -nanopore-raw 1499.fastq.gz -nanopore-raw 1499.failed.fastq.gz \
    -pacbio-raw 1499.pb.fastq.gz \
    -minReadLength=3500 -minOverlapLength=2500 \
    -corOutCoverage=100 -correctedErrorRate=0.11
```

### Purge Haplotigs

We originally tried deduplicating the contigs but the haplotype switching made
completely deduplicating the assembly impossible. We therefore mapped nanpore 
reads to the canu unitigs to allow for the most deduplicated assembly possible.

```bash
minimap2 -t 12 1499.unitigs.fasta 1499.fastq.gz | samtools sort -@ 4 -m 1G -o aligned.bam
purge_haplotigs hist -g 1499.unitigs.fasta -b aligned.bam -t 16
purge_haplotigs cov -l 5 -m 95 -h 190
purge_haplotigs purge -g 1499.unitigs.fasta -c coverage_stats.csv
purge_haplotigs clip -p curated.fasta -h curated.haplotigs.fasta
purge_haplotigs place -p clip.fasta -h clip.haplotigs.fasta -f
```

### Longshot

Longshot had a propensity to phase the divergent haplotype while keeping the two
more similar haplotypes together--which was ideal. Several contigs still had 
hapltype switching, and these contigs needed to be split to properly bin the 
reads. Contigs > ~50kb were manually eyeballed and sorted. contigs < ~50kb with 
roughly 2:1 size ratio were sorted into their anticipated bins (double size = 
the two similar haplotyes, half size = divergent). 

```bash
minimap2 -t 12 clip.FALC.fasta 1499.fastq.gz | samtools sort -@ 4 -m 1G -o clip.FALC.bam
# This script just runs Longshot contig-by-contig over multiple threads
perl multiLongShot.pl  clip.FALC.fasta.fai  clip.FALC.bam  16
```

### Manually fix haplotype switches

Visualise contigs in IGV and identify haplotype switching.
This is a list of haplotype switches and which of hap1/hap2 is the A1+A2 bin
e.g. contig [hap 1/2] breakpoint [1/2] breakpoint ...
Contigs were split at these breakpoints to allow the best possible separation of
A1 + A2 haplomes from B haplome. Manually sorted the BAM files.

```text
000013F [1] 130k [2]
000006F [2] 124k [1]
000005F [1] 73.5k [2]
000009F [1] 82k [2] 127k [1]
000004F [1] 55k [2]
000003F [1] 31k [2] 103k [1] 129k
000014F [2] 82k [1]
000000F [1] 49k [2] 130k [1]
000011F [1] 64k [2]
000010F [2] 93k [1] 134k [2]
000024F [2] 15k [1]
000025F [2] 52k [1]
000027F [1] 74k [2]
000026F [2] 44k [1] 86k [2]
000029F [2] 70k [1]
000031F [1] 59k [2] 69k [1]
000036F [1] 59k [2] 69k [1]
000041F [1] 5k [2] 10k [1]
000045F [1] 80k [2]
000046F [1] 8k [2] 27k [1]
000053F [1] 66k [2]
000055F [2] 3k [1]
000057F [1] 28k [2]
000059F [1] 56k [2]
000060F [1] 43k [2]
000064F [1] 58k [2]
000068F [1] 36k [2] 45k [1]
000074F [1] 46k [2]
000085k [1] 30k [2]
000167F [2] 10k [1] 27k [2]
```

### Read binning

Note: 'haploid' = B haplome and 'diploid' = A1A2 haplomes

```bash
# pull ids for each bin
for i in `ls haploid`; do 
    samtools view haploid/$i | awk '{print $1}'; 
done | gzip - > haploid.reads.gz

for i in `ls diploid`; do 
    samtools view diploid/$i | awk '{print $1}'; 
done | gzip - > diploid.reads.gz

# haploid read bin
perl getBinnedReads.pl 1499.fastq.gz  haploid.reads.gz  unassigned.reads.gz \
    | gzip - > 1499.B.fastq.gz
perl getBinnedReads.pl 1499.fail.fastq.gz  haploid.reads.gz  unassigned.reads.gz \
    | gzip - > 1499.B.fastq.gz
perl getBinnedReads.pl 1499.pb.fastq.gz  haploid.reads.gz  unassigned.reads.gz \
    | gzip - > 1499.B.fastq.gz

# diploid read bin
perl getBinnedReads.pl 1499.fastq.gz  diploid.reads.gz  unassigned.reads.gz \
    | gzip - > 1499.A1A2.fastq.gz
perl getBinnedReads.pl 1499.fail.fastq.gz  diploid.reads.gz  unassigned.reads.gz \
    | gzip - > 1499.A1A2.fastq.gz
perl getBinnedReads.pl 1499.pb.fastq.gz  diploid.reads.gz  unassigned.reads.gz \
    | gzip - > 1499.A1A2.fastq.gz
```

## Reassembly

```bash
# haploid canu assembly
canu -p 1499.B -d bin-B \
    -genomeSize=13m \
    -nanopore-raw 1499.B.fastq.gz \
    -nanopore-raw 1499.B.fail.fastq.gz \
    -pacbio-raw 1499.B.pb.fastq.gz \
    -stopOnReadQuality=false \
    -minReadLength=2500 -minOverlapLength=2000 \
    -correctedErrorRate=0.11

# diploid canu assembly
canu -p 1499.B -d bin-A1A2 \
    -genomeSize=13m \
    -nanopore-raw 1499.A1A2.fastq.gz \
    -nanopore-raw 1499.A1A2.fail.fastq.gz \
    -pacbio-raw 1499.A1A2.pb.fastq.gz \
    -stopOnReadQuality=false \
    -minReadLength=2500 -minOverlapLength=2000 \
    -correctedErrorRate=0.11
```
Purge Haplotigs was run on the read-binned assemblies following the same 
protocol as above. Both assemblies were concatenated to `1499.diploid.fasta`
with the prefixes 'A1A2_' or 'B_' prepended to contig names to denote haplomes.

## Polishing

### Racon

```bash
# map reads
minimap2 -t 12 1499.diploid.fasta ../../1499.fastq.gz -x map-ont \
    | gzip - > NP.paf.gz

# run racon
racon -t 16 ../../1499.fastq.gz NP.paf.gz 1499.diploid.fasta \
    > 1499.racon.fasta
```

### Pilon

```bash
# index fasta
bwa index 1499.racon.fasta

# map reads
bwa mem -t 12 1499.racon.fasta \
    <( zcat ../../1499_S4_L001_R1_001.fastq.gz ../../1499_S29_L001_R1_001.fastq.gz) \
    <(zcat ../../1499_S4_L001_R2_001.fastq.gz ../../1499_S29_L001_R2_001.fastq.gz) \
    | samtools sort -@ 4 -m 1G -o pe.bam

# run pilon
java -jar ~/pilon.jar --genome 1499.racon.fasta --frags pe.bam --threads 16

# rename file
mv pilon.fasta 1499.racon.pilon.fasta

# trim '_pilon' from contig names
cat 1499.racon.pilon.fasta | sed 's/_pilon//' > 1499.diploid.newPolished.fasta
```

## Annotate

We annotated genes with Augustus and used these for initially investigating
the key genes but the results in the manuscript relating to genes affected by
SVs all relate to the AWRI2804 reference annotations. 

```bash
augustus --species=saccharomyces_cerevisiae_S288C 1499.diploid.newPolished.fasta > 1499.gff

# extract prot seqs
cat 1499.gff | perl -e \
    '$/ = "end gene"; 
    while(<>){
        my$gen;
        if(m/start gene (\w+)/){
            $gen=$1;
        }else{
            die;
        }
        s/#//g;
        if(m/protein sequence = \[([\w\s]+)\]/){
            $s=$1;
            $s=~s/\s//g;
            print ">$gen\n$s\n";
        }
    }' \
    > 1499.faa

# prot names
cat 1499.faa \
    | parallel -j 8 --recstart '>' --blocksize 10k --pipe \
        blastp -query - -db /home/mike/blastdb/uniprot_sprot.fasta \
        -evalue 1e-5 -num_alignments 1 -outfmt 6 \
        | gzip - > 1499.uniprotkb.outfmt6.gz

# add the draft names as descriptions
perl addDraftNames.pl 1499.diploid.uniprotkb.outfmt6.gz uniprot_sprot.names.tsv 1499.faa > tmp
mv tmp 1499.faa
```

## SV calling: Sniffles

Most of the Sniffles SV calls overlapped with the below SV calling pipeline, so 
they were left out of the final manuscript.

```bash
# mapping
ngmlr -x ont -t 8 -r 1499.diploid.fasta -q ../1499.fastq.gz | samtools sort -@ 4 -m 1G -o 1499.bam -T tmp
ngmlr -x ont -t 8 -r 1499.diploid.fasta -q ../1499_11A.fastq.gz | samtools sort -@ 4 -m 1G -o 1499_11A.bam -T tmp
ngmlr -x ont -t 8 -r 1499.diploid.fasta -q ../1499_17A.fastq.gz | samtools sort -@ 4 -m 1G -o 1499_17A.bam -T tmp
ngmlr -x ont -t 8 -r 1499.diploid.fasta -q ../1499_52C.fastq.gz | samtools sort -@ 4 -m 1G -o 1499_52C.bam -T tmp

# call SVs with sniffles
for i in 1499 1499_11A 1499_17A 1499_52C; do
    sniffles -m $i.bam -b $i.bedpe -t 8
    sniffles -m $i.bam -v $i.vcf -t 8
    done

# quick filter
for i in 1499 1499_11A 1499_17A 1499_52C; do
    head -1 $i.bedpe > $i.filt.bedpe
    cat $i.bedpe | awk '$12>30{print}' >> $i.filt.bedpe
    done

# manually curate with genome-ribbon: doi.org/10.1101/082123 
```

## SV calling: depth- and SNP-based

```bash
# align 'diploid' queries to haploid reference
nucmer -b 500 -c 40 -d 0.5 -g 200 -l 12 -t 16 \
    AWRI2804_scaffold.fasta 1499.diploid.fasta -p 2804.1499
nucmer -b 500 -c 40 -d 0.5 -g 200 -l 12 -t 16 \
    AWRI2804_scaffold.fasta 1613.H1.fasta -p 2804.1613

# filter all query to best ref
delta-filter -q 2804.1499.delta > 2804.1499.qdelta
delta-filter -q 2804.1613.delta > 2804.1613.qdelta

# create a transformation table
show-coords -HTl 2804.1499.qdelta \
    | awk '$5>1000{if($4<$3){print $10"\t"$1"\t"$2"\t"$11"\t"$4"\t"$3"\t-\t"$9}else{print $10"\t"$1"\t"$2"\t"$11"\t"$3"\t"$4"\t+\t"$9}}' \
    > 1499.to.2804.trn
show-coords -HTl 2804.1613.qdelta \
    | awk '$5>1000{if($4<$3){print $10"\t"$1"\t"$2"\t"$11"\t"$4"\t"$3"\t-\t"$9}else{print $10"\t"$1"\t"$2"\t"$11"\t"$3"\t"$4"\t+\t"$9}}' \
    > 1613.to.2804.trn

# make directories
snakemake  -s findGenomeModifications.py

# manually copy the sequencing data and translation tables 

# find draft conversions
snakemake  -s findGenomeModifications.py  -j 16

# manually curate the draft SV calls in IGV

# produce final plots
snakemake  -s findGenomeModifications.py  -j 16
```





