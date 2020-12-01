### PARAMETERS

# Datasets
pacbio_data="pacbio_reads.fasta"
ont_data="nanopore_reads.fasta"

# Assemblers
CANU="~/Tools/canu/"
FLYE="~/Tools/Flye/bin/flye"
RA="~/Tools/ra/build/bin/ra"
RAVEN="~/Tools/raven/build/bin/raven"
SHASTA="~/Tools/shasta-0.3.0"
WTDBG2="~/Tools/wtdbg2/wtdbg2"
WTPOA="~/Tools/wtdbg2/wtpoa-cns"

# Mapping tools

MINIMAP="~/Tools/minimap2-2.17_x64-linux/minimap2"

# Evaluation tools

ASSEMBLYSTATS="~/Tools/assembly-stats/build/assembly-stats"

# Directories
OUTDIR="full_assemblies"
PBDIR=${OUTDIR}/PacBio_assemblies
ONTDIR=${OUTDIR}/Nanopore_assemblies

# Other
THREADS=1
MAXMEMORY=10 # in GB

mkdir -p $OUTDIR
mkdir $PBDIR
mkdir $ONTDIR

### FULL ASSEMBLIES

for i in {1..5}
do

    # PacBio assemblies
    
    $CANU -d ${PBDIR}/canu_default_0${i} -p canu_default useGrid=false \
        maxMemory=${MAXMEMORY}
        genomeSize=100m -pacbio-raw $pacbio_data
    cp ${PBDIR}/canu_default_0${i}/canu_default.contigs.fasta ${PBDIR}/canu.all.0${i}.fasta

    $FLYE -o ${PBDIR}/flye_default_0${i} -t ${THREADS} -g 100m \
        --pacbio-raw $pacbio_data
    cp ${PBDIR}/flye_default_0${i}/assembly.fasta ${PBDIR}/flye.all.0${i}.fasta

    $RA -x pb -t ${THREADS} $pacbio_data > ${PBDIR}/ra.all.0${i}.fasta

    $RAVEN -t ${THREADS} $pacbio_data > ${PBDIR}/raven.all.0${i}.fasta

    $SHASTA --input $pacbio_data --assemblyDirectory ${PBDIR}/shasta_default_0${i} \
        --Reads.minReadLength 0 --Assembly.consensusCaller Modal --Kmers.k 12
    cp ${PBDIR}/shasta_default_0${i}/Assembly.fasta ${PBDIR}/shasta.all.0${i}.fasta

    $WTDBG2 -x rs -g 100m -i $pacbio_data -fo ${PBDIR}/out
    $WTPOA -i ${PBDIR}/out.ctg.lay.gz -o ${PBDIR}/out.ctg.fa
    $MINIMAP -x map-pb -a ${PBDIR}/out.ctg.fa $pacbio_data | samtools sort > ${PBDIR}/out.ctg.bam
    samtools view ${PBDIR}/out.ctg.bam | $WTPOA -d ${PBDIR}/out.ctg.fa -i - -fo ${PBDIR}/wtdbg2.all.0${i}.fasta

    # Nanopore assemblies

    $CANU -d ${ONTDIR}/canu_default_0${i} -p canu_default useGrid=false \
        maxMemory=${MAXMEMORY}
        genomeSize=100m -nanopore-raw $ont_data
    cp ${ONTDIR}/canu_default_0${i}/canu_default.contigs.fasta ${ONTDIR}/canu.all.0${i}.fasta

    $FLYE -o ${ONTDIR}/flye_default_0${i} -t ${THREADS} -g 100m \
        --nano-raw $ont_data
    cp ${ONTDIR}/flye_default_0${i}/assembly.fasta ${ONTDIR}/flye.all.0${i}.fasta

    $RA -x ont -t ${THREADS} $ont_data > ${ONTDIR}/ra.all.0${i}.fasta

    $RAVEN -t ${THREADS} $ont_data > ${ONTDIR}/raven.all.0${i}.fasta

    $SHASTA --input $ont_data --assemblyDirectory ${ONTDIR}/shasta_default_0${i} \
        --Reads.minReadLength 0 
    cp ${ONTDIR}/shasta_default_0${i}/Assembly.fasta ${ONTDIR}/shasta.all.0${i}.fasta

    $WTDBG2 -x ont -g 100m -i $ont_data -fo ${ONTDIR}/out
    $WTPOA -i ${ONTDIR}/out.ctg.lay.gz -o ${ONTDIR}/out.ctg.fa
    $MINIMAP -x map-ont -a ${ONTDIR}/out.ctg.fa $ont_data | samtools sort > ${ONTDIR}/out.ctg.bam
    samtools view ${ONTDIR}/out.ctg.bam | $WTPOA -d ${ONTDIR}/out.ctg.fa -i - -fo ${ONTDIR}/wtdbg2.all.0${i}.fasta

done

### ASSEMBLIES EVALUATION

$ASSEMBLYSTATS -t ${PBDIR}/*.all.0*fasta > ${PBDIR}/assembly-stats.tsv
$ASSEMBLYSTATS -t ${ONTDIR}/*.all.0*fasta > ${ONTDIR}/assembly-stats.tsv

for i in canu flye ra raven shasta wtdbg2
do

    for j in {1..5}
    do

        # BUSCO
        busco -i ${PBDIR}/${i}.all.0${j}.fasta -o busco_${i}_${j} -l metazoa_odb10 -m genome
        mv busco_${i}_${j} $PBDIR
        busco -i ${ONTDIR}/${i}.all.0${j}.fasta -o busco_${i}_${j} -l metazoa_odb10 -m genome
        mv busco_${i}_${j} $ONTDIR

    done

done