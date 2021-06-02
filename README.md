# Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms
Supporting scripts for the paper **"Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms"** by Guiglielmoni et al., available in [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.16.993428v1). 

In this paper, we compare assemblies of PacBio and Nanopore reads with seven assemblers, and evaluate their efficiency at producing properly collapsed haploid assemblies. We also try to reduce uncollapsed haplotypes by using a subset of filtered reads, and haplotig-purging tools. We also present the tool [HapPy](https://github.com/AntoineHo/Happy) to evaluate the correct haploid representation of a diploid genome.

## Table of contents

* [Requirements](https://github.com/nadegeguiglielmoni/overcoming_uncollapsed_haplotypes_in_long_read_assemblies/blob/main/README.md#requirements)
* [Datasets](https://github.com/nadegeguiglielmoni/overcoming_uncollapsed_haplotypes_in_long_read_assemblies/blob/main/README.md#datasets)
* [Assembly of full datasets](https://github.com/nadegeguiglielmoni/overcoming_uncollapsed_haplotypes_in_long_read_assemblies/blob/main/README.md#assembly-of-full-datasets)
* [Read filtering](https://github.com/nadegeguiglielmoni/overcoming_uncollapsed_haplotypes_in_long_read_assemblies/blob/main/README.md#read-filtering)
* [Impact of read-depth](https://github.com/nadegeguiglielmoni/overcoming_uncollapsed_haplotypes_in_long_read_assemblies/blob/main/README.md#impact-of-read-depth)
* [Assembly evaluation](https://github.com/nadegeguiglielmoni/overcoming_uncollapsed_haplotypes_in_long_read_assemblies/blob/main/README.md#assembly-evaluation)

## Requirements

### Long-read assemblers
* [Canu](https://github.com/marbl/canu)
* [Flye](https://github.com/fenderglass/Flye)
* [NextDeNovo](https://github.com/Nextomics/NextDenovo)
* [Ra](https://github.com/lbcb-sci/ra)
* [Raven](https://github.com/lbcb-sci/raven)
* [Shasta](https://github.com/chanzuckerberg/shasta)
* [wtdbg2](https://github.com/ruanjue/wtdbg2)

### Haplotigs-purging tools
* [HaploMerger2](https://github.com/mapleforest/HaploMerger2)
* [purge_dups](https://github.com/dfguan/purge_dups)
* [purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)

### Sampling depth subsets
* [BBmap](https://sourceforge.net/projects/bbmap/)

### Evaluation
* [assembly-stats](https://github.com/sanger-pathogens/assembly-stats)
* [BUSCO](https://busco.ezlab.org/)
* [KAT](https://github.com/TGAC/KAT)
* [tinycov](https://github.com/cmdoret/tinycov)
* [HapPy](https://github.com/AntoineHo/Happy)

### Other dependencies

Python packages:
* Biopython

## Datasets

We tested assemblies on two long-read datasets for the bdelloid rotifer *Adineta vaga*:
* PacBio dataset (23.5 Gb, N50 = 11.6 kb) designated as **pacbio_data**
* Nanopore dataset (17.5 Gb, N50 = 18.8 kb) designated as **ont_data**
* Illumina dataset (11.4 Gb, 2*250 bp) designated as **illumina_end1** and **illumina_end2**

## Assembly of full datasets

### PacBio reads 

We ran the assemblers on the PacBio data in the following manner. Assemblies were repeated 5 times to test the robustness of the output assemblies.

#### Canu
```bash
canu -d corrected -p corrected genomeSize=100m useGrid=false -pacbio-raw pacbio_data
for i in {1..5}
do
    canu -d out -p out genomeSize=100m useGrid=false \
        -pacbio-corrected corrected/corrected.correctedReads.fasta
    cp out/out.contigs.fasta canu_assembly_${i}.fasta
done
```

#### Flye
```bash
for i in {1..5}
do
    flye -o out -g 100m --pacbio-raw pacbio_data
    cp out/assembly.fasta flye_assembly_${i}.fasta
done
```

#### NextDeNovo

```bash
echo pacbio_data > input.fofn
seq_stat input.fofn -g 100Mb -d 150 > stats.txt
for i in {1..5}
do
    nextDenovo run.cfg
done
```

#### Ra
```bash
for i in {1..5}
do
    ra -x pb pacbio_data > ra_assembly_${i}.fasta
done
```

#### Raven
```bash
for i in {1..5}
do
    raven pacbio_data > raven_assembly_${i}.fasta
done
```

#### Shasta
```bash
for i in {1..5}
do
    shasta --input pacbio_data --Reads.minReadLength 0 --assemblyDirectory out \
        --Assembly.consensusCaller Modal --Kmers.k 12
    cp out/Assembly.fasta shasta_assembly_${i}.fasta
done
```

#### wtdbg2
```bash
for i in {1..5}
do
    wtdbg2 -x rs -g 100m -i pacbio_data -fo out
    wtpoa-cns -i out.ctg.lay.gz -o out.ctg.fa
    minimap2 -x map-pb -a out.ctg.fa pacbio_data | samtools sort > out.ctg.bam
    samtools view out.ctg.bam | wtpoa-cns -d out.ctg.fa -i - -fo wtdbg2_assembly_${i}.fasta
done
```

### Nanopore reads 

We ran the assemblers on the Nanopore data in the following manner. Assemblies were repeated 5 times to test the robustness of the output assemblies.

#### Canu
```bash
canu -d corrected -p corrected genomeSize=100m useGrid=false -nanopore-raw ont_data
for i in {1..5}
do
    canu -d out -p out genomeSize=100m useGrid=false \
        -nanopore-corrected corrected/corrected.correctedReads.fasta
    cp out/out.contigs.fasta canu_assembly_${i}.fasta
done
```

#### Flye
```bash
for i in {1..5}
do
    flye -o out -g 100m --nano-raw ont_data
    cp out/assembly.fasta flye_assembly_${i}.fasta
done
```

#### NextDeNovo

```bash
echo ont_data > input.fofn
seq_stat input.fofn -g 100Mb -d 150 > stats.txt
for i in {1..5}
do
    nextDenovo run.cfg
done
```

#### Ra
```bash
for i in {1..5}
do
    ra -x ont ont_data > ra_assembly_${i}.fasta
done
```

#### Raven
```bash
for i in {1..5}
do
    raven ont_data > raven_assembly_${i}.fasta
done
```

#### Shasta
```bash
for i in {1..5}
do
    shasta --input ont_data --Reads.minReadLength 0 --assemblyDirectory out 
    cp out/Assembly.fasta shasta_assembly_${i}.fasta
done
```

#### wtdbg2
```bash
for i in {1..5}
do
    wtdbg2 -x ont -g 100m -i ont_data -fo out
    wtpoa-cns -i out.ctg.lay.gz -o out.ctg.fa
    minimap2 -x map-ont -a out.ctg.fa ont_data | samtools sort > out.ctg.bam
    samtools view out.ctg.bam | wtpoa-cns -d out.ctg.fa -i - -fo wtdbg2_assembly_${i}.fasta
done
```

## Read filtering

The reads were filtered with a python script.

```bash
python3 filter_reads.py --min 15000 --fasta pacbio_data --out pacbio.min15000.fasta
python3 filter_reads.py --min 30000 --fasta ont_data --out ont.min30000.fasta  
```

The assemblies were run on these filtered datasets in the same manner as for assemblies of the full datasets.

## Haplotig purging 

Assemblies were post-processed with HaploMerger2, purge_dups, purge_haplotigs (one replicate per assembler). L, M, H are selected depending on the histogram produced by ```purge_haplotigs hist```.

### PacBio assemblies

```bash
for i in canu flye ra raven shasta wtdbg2
do
    minimap2 -ax map-pb ${i}_assembly_1.fasta pb_data --secondary=no | samtools sort -o ${i}_all.ali.sorted.bam -T tmp.ali
    samtools index ${i}_all.ali.sorted.bam
    samtools faidx ${i}_assembly_1.fasta
    purge_haplotigs hist -b ${i}_all.ali.sorted.bam -g ${i}_assembly_1.fasta
    purge_haplotigs cov -i ${i}_all.ali.sorted.bam -l L -m M -h H -o ${i}_all.cov_stats.csv
    purge_haplotigs purge -g ${i}_assembly_1.fasta -c ${i}_all.cov_stats.csv -o ${i}_assembly_1.purged.fasta
done
```

### Nanopore assemblies

```bash
for i in canu flye ra raven shasta wtdbg2
do
    minimap2 -ax map-ont ${i}_assembly_1.fasta ont_data --secondary=no | samtools sort -o ${i}_all.ali.sorted.bam -T tmp.ali
    samtools index ${i}_all.ali.sorted.bam
    samtools faidx ${i}_assembly_1.fasta
    purge_haplotigs hist -b ${i}_all.ali.sorted.bam -g ${i}_assembly_1.fasta
    purge_haplotigs cov -i ${i}_all.ali.sorted.bam -l L -m M -h H -o ${i}_all.cov_stats.csv
    purge_haplotigs purge -g ${i}_assembly_1.fasta -c ${i}_all.cov_stats.csv -o ${i}_assembly_1.purged.fasta
done
```

## Impact of read-depth

Subsets were randomly sampled from PacBio and Nanopore datasets to reach a certain read-depth. The script **reformat.sh** from BBmap was used to generate these subsets. The target number of bases if based on the estimated haploid genome size of 102.3 Mb. 

```bash
for i in {1..5}
do
    reformat.sh in=pacbio_data out=pacbio.subset.10X.0${i}.fasta samplebasestarget=1023000000
    reformat.sh in=pacbio_data out=pacbio.subset.20X.0${i}.fasta samplebasestarget=2046000000
    reformat.sh in=pacbio_data out=pacbio.subset.30X.0${i}.fasta samplebasestarget=3069000000
    reformat.sh in=pacbio_data out=pacbio.subset.40X.0${i}.fasta samplebasestarget=4092000000
    reformat.sh in=pacbio_data out=pacbio.subset.50X.0${i}.fasta samplebasestarget=5115000000
    reformat.sh in=pacbio_data out=pacbio.subset.60X.0${i}.fasta samplebasestarget=6138000000
    reformat.sh in=pacbio_data out=pacbio.subset.80X.0${i}.fasta samplebasestarget=8184000000
    reformat.sh in=pacbio_data out=pacbio.subset.100X.0${i}.fasta samplebasestarget=10230000000

    reformat.sh in=ont_data out=ont.subset.10X.0${i}.fasta samplebasestarget=1023000000
    reformat.sh in=ont_data out=ont.subset.20X.0${i}.fasta samplebasestarget=2046000000
    reformat.sh in=ont_data out=ont.subset.30X.0${i}.fasta samplebasestarget=3069000000
    reformat.sh in=ont_data out=ont.subset.40X.0${i}.fasta samplebasestarget=4092000000
    reformat.sh in=ont_data out=ont.subset.50X.0${i}.fasta samplebasestarget=5115000000
    reformat.sh in=ont_data out=ont.subset.60X.0${i}.fasta samplebasestarget=6138000000
    reformat.sh in=ont_data out=ont.subset.80X.0${i}.fasta samplebasestarget=8184000000
    reformat.sh in=ont_data out=ont.subset.100X.0${i}.fasta samplebasestarget=10230000000
done
```

The assemblies were run on these subsets in the same manner as for assemblies of the full datasets. For Canu assemblies, the parameter stopOnLowCoverage was set to 1 to allow runs on low-depth datasets.

## Assembly evaluation

### Basic statistics
Assembly size and N50 were computed with **assembly-stats**

```bash
assembly-stats assembly.fasta
```

### BUSCO completeness
The BUSCO completeness was computed using the dataset metazoa odb10. 

```bash
busco -i assembly.fasta -o busco_output -l metazoa_odb10 -m genome
```

### *k*-mer completeness
The *k*-mer completeness is calculated against a low-error rate Illumina dataset.

```bash
kat comp -o kat_output 'illumina_end1 illumina_end2' assembly.fasta
```

### tinycov
First we mapped the long reads to the assemblies, PacBio reads for PacBio assemblies, Nanopore reads for Nanopore assemblies.

```bash
minimap2 -ax map-pb assembly.fasta pacbio_data --secondary=no | samtools sort -o aligned.bam -T tmp.ali
samtools index mapping.bam
```
or

```bash
minimap2 -ax map-ont assembly.fasta ont_data --secondary=no | samtools sort -o aligned.bam -T tmp.ali
samtools index mapping.bam
```

Then the mapping data was used to study the coverage with **tinycov**.

```bash
tinycov covplot -r 20000 -t cov.txt aligned.bam
```

### HapPy

The mapping data is used again for HapPy, that also relies on coverage.

```bash
main.py coverage -d happy_output aligned.bam
main.py estimate -S 102M -O happy_out/happy.stats happy_output/aligned.bam.hist 
```
