# MitoSAlt-custom
A wrapper with enhancements around MitoSalt MT SV caller tool

# Quick start

Build Docker image

```bash
docker build -t mitosalt .
```

Build Singularity image if needed

```bash
singularity build mitosalt.sif docker-daemon://mitosalt:latest
```

Run the full MitoSAlt pipeline, please mind that pipeline expects you to mount to `/output`

```bash
# temp folders to explicitely provide to singularity 
mkdir -p tmp_home tmp

# output prefix
sample=$1
# read files basename (assuming they are located according to mount below)
R1=$2
R2=$3

singularity exec --no-home \
 --bind $(pwd)/tmp_home:/tmp_home \
 --bind $(pwd)/tmp:/tmp \
 --bind $(pwd)/data/reads/:/data \
 --bind $(pwd)/genome:/output/genome \ # genome ref and indices built according to MitoSAlt pipeline requirements
 --bind $(pwd)/output/${sample}:/output \ # this mount is needed
 --env HOME=/tmp_home,TMPDIR=/tmp \
 --pwd /output \
 mitosalt.sif  bash -c "cd /output && mkdir -p bam bin bw indel log plot tab \
   && perl /opt/MitoSAlt_1.1.1/MitoSAlt1.1.1.pl \
    /opt/MitoSAlt_1.1.1/config_human.txt \
    /data/${R1} \
    /data/${R2} \
    ${sample}"
```

## SAltShaker - Visualizer module for MitoSAlt results

A Python port of the original MitoSAlt R plotting and classification script, with enhanced event pattern classification and visualization. Read more [here](./saltshaker/README.md).
