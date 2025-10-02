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

```
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

or run the classifier and visualisation script separately on the outputs of MitoSAlt main pipeline (before the delplot.R script): 

```
# output prefix
sample=$1
# default settings from MitoSAlt config file
msize=16569
orihs=16081
orihe=407
orils=5730
orile=5763
sizelimit=10000
hplimit=0.01
MT_fasta=human_mt_rCRS.fasta
flank=15

singularity exec --no-home \
	--bind $(pwd)/data:/data \
	--bind $(pwd)/genome:/output/genome \
	--bind $(pwd)/output/${sample}:/output \
	--bind $(pwd)/output/${sample}:/tmp \
	mitosalt.sif bash -c "python /opt/MitoSAlt_1.1.1/mitosalt_visualizer.py \
			${msize} ${orihs} ${orihe} ${orils} ${orile} ${sizelimit} \
			/output/indel/${sample}.cluster /output/indel/${sample}.breakpoint \
			${sample} ${hplimit} /output/genome/${MT_fasta} ${flank} \
			--blacklist /output/genome/mt_blacklist.bed --output-dir /output"

```

## SAltShaker - MitoSAlt visualizer module

A Python port of the original MitoSAlt R plotting and classification script, with enhanced event pattern classification and visualization. Read more [here](./saltshaker/README.md).