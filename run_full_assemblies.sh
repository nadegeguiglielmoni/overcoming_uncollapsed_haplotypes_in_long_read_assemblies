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

    $RA -x pb -t ${THREADS} pacbio_data > ${PBDIR}/ra.all.0${i}.fasta

    $RAVEN -t ${THREADS} pacbio_data > ${PBDIR}/raven.all.0${i}.fasta

    # Nanopore assemblies

    $CANU -d ${ONTDIR}/canu_default_0${i} -p canu_default useGrid=false \
        maxMemory=${MAXMEMORY}
        genomeSize=100m -nanopore-raw $ont_data
    cp ${ONTDIR}/canu_default_0${i}/canu_default.contigs.fasta ${ONTDIR}/canu.all.0${i}.fasta

    $FLYE -o ${ONTDIR}/flye_default_0${i} -t ${THREADS} -g 100m \
        --nano-raw $ont_data
    cp ${ONTDIR}/flye_default_0${i}/assembly.fasta ${ONTDIR}/flye.all.0${i}.fasta

    $RA -x ont -t ${THREADS} ont_data > ${ONTDIR}/ra.all.0${i}.fasta

    $RAVEN -t ${THREADS} ont_data > ${ONTDIR}/raven.all.0${i}.fasta

done