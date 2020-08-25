#!/usr/bin/env snakemake

# Copyright (c) 2020 Michael Roach (Australian Wine Research Institute)


# Note: this was only designed as a once-off and will probably not work for 
    # your data without a lot of tweaking.


import gzip
import os
import numpy as np



# usage
# snakemake  -s findGenomeModifications.py  -j 16
# manually curate files in /OUT/MODS/ and rerun



# pipeline procedure
    # map reads for all samples to respective reference genomes
    # calculate windowed coverage
    # call snps and map to same windows
    # flag candidate gene conversion/duplication/deletion regions
        # (based on snp/coveraged changes compared to reference)
    # manually curated these candidate regions
    # transform to haploid reference 2804
    # create files for plotting
    # plot to visualize



# SAMPLE PREFIXES
s1499 = ['1499', '1499_50g_6A', '1499_100g_5A', '1499_100g_52C', '1499_11A', '1499_17A', '1499_52C']
s1613 = ['1613', '1613_50g_31B', '1613_100g_21B', '1613_100g_8B', '1613_10A', '1613_13A', '1613_14A']



# STATIC FILES
    # reference genomes
fasta1499 = 'REFS/1499.fasta'
fasta1613 = 'REFS/1613.fasta'
fasta2804 = 'REFS/2804.fasta'
    # translation tables
trn1499 = 'REFS/1499.to.2804.trn'
trn1613 = 'REFS/1613.to.2804.trn'
    # directory structure
dirs = ['BAMS', 'LOG', 'OUT', 'OUT/SNPCOV', 'MISC', 'MISC/MODS', 'MISC/TRN', 'OUT/MODS', 'REFS', 'SEQ']



# PROGRAMS + PARAMETERS, these should be in a config file
nglmrRun = 'ngmlr -x ont -t 12'                                                             # ... -r ref.fasta -q query.fastq.gz | samtools sort -@ 4 -m 1G -o query.bam -T tmp
samtoolsMpileupRun = 'samtools mpileup -A -x -B'                                            # ... -f ref.fasta  sample.bam
varscanRun = 'java -jar ~/VarScan.v2.3.9.jar mpileup2snp --p-value 1e-5 --output-vcf'       #  samtools mpileup | varscan |  grep -P '#|HET=1'  > output
bedtoolsMapRun = "bedtools map -c 4"                                                        # ... -a sample.bed -b sample.snp.bed
    # snp cov diff script
snpCovDiffRun = "perl scripts/snpCovDiff.pl"                                                        # ... -rd ref.seqCov  -qd query.seqCov  -rs ref.snp  -rc ref.cov  -qs evolv.snp  -qc evolv.cov
    # different cutoffs for 1499 and 1613, output still requires manual curation
genCnvAwk1499 = "awk '$5 >= -20 {{ print }}' | awk '$4 <= -30 {{ print }}'"                 # snpCovDiff | genCnvAwk | snpCovPostProc > output   # these output files need manual curation
genDupAwk1499 = "awk '$5 >= 30 {{ print }}'"
genDelAwk1499 = "awk '$5 <= -30 {{ print }}'"
genCnvAwk1613 = "awk '$5 >= -20 {{ print }}' | awk '$4 <= -50 {{ print }}'"                 # snpCovDiff | genCnvAwk | snpCovPostProc > output   # these output files need manual curation
genDupAwk1613 = "awk '$5 >= 50 {{ print }}'"
genDelAwk1613 = "awk '$5 <= -50 {{ print }}'"
    # post processing for easy copy-pasting into IGV
snpCovPostProc = "| sort -k1,1 -k2,2n | bedtools merge | sed 's/\t0\t/\t1\t/'"
#
bedSort = "awk '{{print $4\"\\t\"$5\"\\t\"$6}}' | sort -k1,1 -k2,2n"
transBed = "scripts/transform_BED.pl"                                     # ... transTable  bedfile
bed2geom = "scripts/bed2geom_ideo.pl"                                     # ... ref.fai  bed1  bed2 ...
plotIdeo = "scripts/GcGduGde.plot.Rscript"



# RULES
rule all:
    input:
        dirs,
        'plot.svg'

rule makeDirectories:
    output:
        dirs
    shell:
        "for i in {output}; do mkdir -p $i; done"

rule mapReads:
    input:
        q = "SEQ/{sample}.fastq.gz",
        r = "REFS/{ref}.fasta"
    output:
        "BAMS/{sample}.{ref}.bam"
    threads:
        16
    log:
        n = "LOG/nglmr.stderr",
        s = "LOG/samtools.stderr"
    shell:
        "{nglmrRun} -r {input.r} -q {input.q} 2> {log.n} | samtools sort -@ 4 -m 1G -o {output} 2> {log.s}"

rule bedtoolsMakewindows:
    input:
        "REFS/{ref}.fasta.fai"
    output:
        "REFS/{ref}.windows.bed"
    log:
        "LOG/bedtools.makeWind.{ref}.stderr"
    run:
        fai = open(input[0], 'r')
        genLen = 0
        for line in fai:
            l = re.split('\s+', line)
            genLen += int(l[1])
        fai.close()
        windowSize = np.ceil(genLen / 1000000) * 1000
        stepSize = int(windowSize / 2)
        shell("bedtools makewindows -g {input} -w {windowSize} -s {stepSize} 2> {log} > {output}")

rule samtoolsBedcov:
    input:
        b = "BAMS/{sample}.{ref}.bam",
        w = "REFS/{ref}.windows.bed",
        i = "BAMS/{sample}.{ref}.bam.bai"
    output:
        "OUT/SNPCOV/{sample}.{ref}.cov"
    log:
        "LOG/samtools.bedcov.{sample}.stderr"
    shell:
        "samtools bedcov {input.w} {input.b} 2> {log}\
        | awk '{{print $1\"\t\"$2\"\t\"$3\"\t\"($4 / ($3-$2))}}'\
        > {output}"

rule callSNPs:
    input:
        r = "REFS/{ref}.fasta",
        b = "BAMS/{sample}.{ref}.bam", 
        ib = "BAMS/{sample}.{ref}.bam.bai",
        ri = "REFS/{ref}.fasta.fai"
    output:
        "MISC/{sample}.{ref}.snp.bed"
    log:
        s = "LOG/samtools.mpileup.{sample}.stderr",
        v = "LOG/varscan.{sample}.stderr"
    threads:
        2
    shell:
        "{samtoolsMpileupRun} -f {input.r} {input.b} 2> {log.s} \
            | {varscanRun} 2> {log.v} \
            | grep -P '#|HET=1'  \
            | grep -v \# \
            | awk '{{print $1\"\t\"$2\"\t\"$2\"\t\"1}}' \
            > {output}"

rule bedtoolsMap:
    input:
        b = "REFS/{ref}.windows.bed",
        s = "MISC/{sample}.{ref}.snp.bed"
    output:
        "OUT/SNPCOV/{sample}.{ref}.snp"
    log:
        "LOG/bedtools.map.{sample}.stderr"
    threads:
        2
    shell:
        "{bedtoolsMapRun} -a {input.b} -b {input.s} 2> {log} \
            | awk '{{print $1\"\t\"$2\"\t\"$3\"\t-\"$4}}' \
            | sed 's/-\./0/' \
            > {output}"

rule sampleSeqCoverage:
    input:
        "SEQ/{sample}.fastq.gz"
    output:
        "MISC/{sample}.seqCov"
    threads:
        2
    run:
        totalLen = 0
        f = gzip.open(input[0])
        for l in f:
            s = f.readline()
            totalLen += len(s) - 1
            f.readline()
            f.readline()
        o = open(output[0], 'w')
        o.write(str(totalLen) + "\n")

rule snpCovDiff1499:
    input:
        rd = "MISC/1499.seqCov",
        qd = "MISC/{sample}.seqCov",
        rs = "OUT/SNPCOV/1499.1499.snp",
        rc = "OUT/SNPCOV/1499.1499.cov",
        qs = "OUT/SNPCOV/{sample}.1499.snp",
        qc = "OUT/SNPCOV/{sample}.1499.cov"
    output:
        cn = "OUT/MODS/{sample}.1499.GC.tsv.CurateMeAndRerun",
        dp = "OUT/MODS/{sample}.1499.GDup.tsv.CurateMeAndRerun",
        dl = "OUT/MODS/{sample}.1499.GDel.tsv.CurateMeAndRerun",
        r  = "MISC/MODS/{sample}.1499.snpCovDiff.tsv"
    log:
        "LOG/snpCovDiff.{sample}.stderr"
    threads:
        4
    shell:
        "{snpCovDiffRun} -rd {input.rd} -qd {input.qd} -rs {input.rs} -rc {input.rc} -qs {input.qs} -qc {input.qc} 2> {log} \
        | tee \
            >( {genCnvAwk1499} {snpCovPostProc} > {output.cn} ) \
            >( {genDupAwk1499} {snpCovPostProc} > {output.dp} ) \
            >( {genDelAwk1499} {snpCovPostProc} > {output.dl} ) \
            > {output.r}"

rule snpCovDiff1613:
    input:
        rd = "MISC/1613.seqCov",
        qd = "MISC/{sample}.seqCov",
        rs = "OUT/SNPCOV/1613.1613.snp",
        rc = "OUT/SNPCOV/1613.1613.cov",
        qs = "OUT/SNPCOV/{sample}.1613.snp",
        qc = "OUT/SNPCOV/{sample}.1613.cov"
    output:
        cn = "OUT/MODS/{sample}.1613.GC.tsv.CurateMeAndRerun",
        dp = "OUT/MODS/{sample}.1613.GDup.tsv.CurateMeAndRerun",
        dl = "OUT/MODS/{sample}.1613.GDel.tsv.CurateMeAndRerun",
        r  = "MISC/MODS/{sample}.1613.snpCovDiff.tsv"
    log:
        "LOG/snpCovDiff.{sample}.stderr"
    threads:
        4
    shell:
        "{snpCovDiffRun} -rd {input.rd} -qd {input.qd} -rs {input.rs} -rc {input.rc} -qs {input.qs} -qc {input.qc} 2> {log} \
        | tee \
            >( {genCnvAwk1613} {snpCovPostProc} > {output.cn} ) \
            >( {genDupAwk1613} {snpCovPostProc} > {output.dp} ) \
            >( {genDelAwk1613} {snpCovPostProc} > {output.dl} ) \
            > {output.r}"

rule sortTrn:
    input:
        '{trnFile}.trn'
    output:
        '{trnFile}.trn.isec'
    shell:
        "cat {input} | {bedSort} > {output}"

rule trans1499to2804:
    input:
        t = trn1499,
        i = trn1499 + '.isec',
        b = 'OUT/MODS/{sample}.1499.{blagh}.CurateMeAndRerun'
    output:
        'MISC/TRN/{sample}.1499.{blagh}.2804'
    log:
        l1 = 'LOG/bedtoolsIsec.{sample}.{blagh}.stderr',
        l2 = 'LOG/transBed.{sample}.{blagh}.stderr'
    shell:
        "bedtools intersect -a {input.b} -b {input.i} 2> {log.l1} > {input.b}.isec; \
        {transBed} {input.t} {input.b}.isec 2> {log.l2} > {output}"

rule trans1613to2804:
    input:
        t = trn1613,
        i = trn1613 + '.isec',
        b = 'OUT/MODS/{sample}.1613.{blagh}.CurateMeAndRerun'
    output:
        'MISC/TRN/{sample}.1613.{blagh}.2804'
    log:
        l1 = 'LOG/bedtoolsIsec.{sample}.{blagh}.stderr',
        l2 = 'LOG/transBed.{sample}.{blagh}.stderr'
    shell:
        "bedtools intersect -a {input.b} -b {input.i} 2> {log.l1} > {input.b}.isec; \
        {transBed} {input.t} {input.b}.isec 2> {log.l2} > {output}"

rule bedToGeomIdeo:
    input:
        f = fasta2804 + '.fai',
        gc1 = expand('MISC/TRN/{sample}.1499.GC.tsv.2804', sample=s1499),
        gc2 = expand('MISC/TRN/{sample}.1613.GC.tsv.2804', sample=s1613),
        du1 = expand('MISC/TRN/{sample}.1499.GDup.tsv.2804', sample=s1499),
        du2 = expand('MISC/TRN/{sample}.1613.GDup.tsv.2804', sample=s1613),
        de1 = expand('MISC/TRN/{sample}.1499.GDel.tsv.2804', sample=s1499),
        de2 = expand('MISC/TRN/{sample}.1613.GDel.tsv.2804', sample=s1613)
    output:
        fasta2804 + '.fai.geom_path',
        expand('MISC/TRN/{sample}.1499.{blagh}.2804.geom_area', sample=s1499, blagh=['GC.tsv','GDup.tsv','GDel.tsv']),
        expand('MISC/TRN/{sample}.1613.{blagh}.2804.geom_area', sample=s1613, blagh=['GC.tsv','GDup.tsv','GDel.tsv'])
    log:
        'LOG/bed2geom.1499.stderr'
    shell:
        "{bed2geom} {input.f} {input.gc1} {input.gc2} 2> {log}; \
         {bed2geom} {input.f} {input.du1} {input.du2} 2> {log}; \
         {bed2geom} {input.f} {input.de1} {input.de2} 2> {log}; "

rule mergeGeoms:
    input:
        gc = expand('MISC/TRN/{sample}.1499.GC.tsv.2804.geom_area', sample=s1499),
        du = expand('MISC/TRN/{sample}.1499.GDup.tsv.2804.geom_area', sample=s1499),
        de = expand('MISC/TRN/{sample}.1499.GDel.tsv.2804.geom_area', sample=s1499),
        gc2 = expand('MISC/TRN/{sample}.1613.GC.tsv.2804.geom_area', sample=s1613),
        du2 = expand('MISC/TRN/{sample}.1613.GDup.tsv.2804.geom_area', sample=s1613),
        de2 = expand('MISC/TRN/{sample}.1613.GDel.tsv.2804.geom_area', sample=s1613)
    output:
        gc = 'OUT/PLOT/GC.geom_poly',
        du = 'OUT/PLOT/GDup.geom_poly',
        de = 'OUT/PLOT/GDel.geom_poly'
    shell:
        "cat {input.gc} > {output.gc}; \
         cat {input.du} > {output.du}; \
         cat {input.de} > {output.de}; \
         cat {input.gc2} >> {output.gc}; \
         cat {input.du2} >> {output.du}; \
         cat {input.de2} >> {output.de}"

rule plotGeoms:
    input:
        fasta2804 + '.fai.geom_path',
        'OUT/PLOT/GC.geom_poly',
        'OUT/PLOT/GDup.geom_poly',
        'OUT/PLOT/GDel.geom_poly'
    output:
        "plot.svg"
    shell:
        '{plotIdeo} {output} {input}'


rule fastaIndex:
    input:
        "{fastaFile}"
    output:
        "{fastaFile}.fai"
    shell:
        "samtools faidx {input}"


rule bamIndex:
    input:
        "{bamFile}.bam"
    output:
        "{bamFile}.bam.bai"
    shell:
        "samtools index {input}"








