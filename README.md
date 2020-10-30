# Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms
Supporting scripts for the paper **"Overcoming uncollapsed haplotypes in long-read assemblies of non-model organisms"** by Guiglielmoni et al., available in [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.16.993428v1). 

In this paper, we compare assemblies of PacBio and Nanopore reads with seven assemblers, and evaluate their efficiency at producing properly collapsed haploid assemblies. We also try to reduce uncollapsed haplotypes by using a subset of filtered reads, and haplotigs-purging tools. We also present the tool [HapPy](https://github.com/AntoineHo/Happy) to evaluate the correct haploid representation of a diploid genome.

## Requirements

### Long-read assemblers
* [Canu]
* [Flye]
* [NextDeNovo]
* [Ra]
* [Raven]
* [Shasta]
* [wtdbg2]

### Haplotigs-purging tools
* [HaploMerger2]
* [purge_dups]
* [purge_haplotigs]

### Sampling depth subsets
* [BBmap]

### Assembly evaluation
* [assembly-stats]
* [BUSCO]
* [KAT]
* [tinycov]
* [HapPy]

## Datasets

We tested assemblies on two long-read datasets for the bdelloid rotifer *Adineta vaga*:
* PacBio dataset (23.5 Gb, N50 = 11.6 kb) designated as **pacbio_data**
* Nanopore dataset (17.5 Gb, N50 = 18.8 kb) designated as **ont_data**

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