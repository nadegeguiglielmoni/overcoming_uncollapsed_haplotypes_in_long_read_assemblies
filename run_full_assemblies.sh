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
WTDBG2=""
WTPOA=""

# Other
THREADS=1
OUTDIR="full_assemblies"

mkdir -p $OUTDIR

for i in {1..5}
do

    # PacBio assemblies

    # Nanopore assemblies


done