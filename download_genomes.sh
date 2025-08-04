#!/bin/bash

# Human reference genome download and indexing script for MitoSAlt Docker container
# Run this script inside the container to download and index human reference genomes

set -e
set -u
set -o pipefail

# Variables from original setup script
TMPFILE=$(mktemp)
now=$(date)

# Human genome data URLs and filenames
HG19_URL=https://www.dropbox.com/s/e1xwzye9hieewxz/human_g1k_v37.fasta.gz
HG19_V=hg19_g1k.fasta
MTRCRS_V=human_mt_rCRS.fasta
HG19S_V=hg19_g1k.size

# Default thread count
THREADS=${THREADS:-4}

# Function to ask yes/no questions
ask_yes_no() {
    local question=$1
    local response
    while true; do
        echo -e "\x1b[97;41m $question (yes/no): \x1b[m"
        read response
        case $response in
            [Yy]* ) return 0;;
            [Nn]* ) return 1;;
            * ) echo "Please answer yes or no.";;
        esac
    done
}

echo -e "\x1b[30;43m MitoSAlt Human Reference Genome Setup \x1b[m"
echo -e "\x1b[30;43m Current working directory: $(pwd) \x1b[m"

# Ensure we're in the MitoSAlt directory - check for MitoSAlt Perl scripts
if [[ ! -f "MitoSAlt1.1.1.pl" ]] && [[ ! -f "MitoSAlt_SE1.1.1.pl" ]]; then
    echo "Error: MitoSAlt Perl scripts not found. Please run this script from the MitoSAlt directory."
    echo "Current directory contents:"
    ls -la
    exit 1
fi

# Ask user what to download
if ask_yes_no "Download human genome and extract Mt sequence?"; then
    HUMAN=yes
else
    HUMAN=no
fi

if ask_yes_no "Build human indexes?"; then
    HINDEX=yes
else
    HINDEX=no
fi

echo -e "\x1b[97;41m Number of threads to use for building indexes (default: $THREADS): \x1b[m"
read user_threads
if [[ -n "$user_threads" ]]; then
    THREADS=$user_threads
fi

echo -e "\n\n\x1b[30;43m $now: Starting human genome download and indexing with $THREADS threads \x1b[m"

# Download human genome and extract MT genome
if [ "$HUMAN" == "yes" ]; then
    echo -e "\x1b[97;41m $now: Downloading human genome... \x1b[m"
    wget -O $TMPFILE $HG19_URL
    gunzip -c $TMPFILE > genome/$HG19_V
    faSize -detailed genome/$HG19_V | egrep -v 'GL|MT' > genome/$HG19S_V
    rm $TMPFILE

    echo 'MT' > tmp.txt
    faSomeRecords genome/$HG19_V tmp.txt genome/$MTRCRS_V
    rm tmp.txt
    echo -e "\x1b[92m Human genome downloaded successfully \x1b[m"
fi

# Build human indexes
if [ "$HINDEX" == "yes" ]; then
    echo -e "\x1b[97;41m $now: Building human indexes... \x1b[m"
    hisat2-build -p $THREADS genome/$HG19_V genome/hg19_g1k
    lastdb -uNEAR genome/human_mt_rCRS genome/$MTRCRS_V
    samtools faidx genome/$HG19_V
    samtools faidx genome/$MTRCRS_V
    echo -e "\x1b[92m Human indexes built successfully \x1b[m"
fi

echo -e "\n\n\x1b[30;43m $now: Human reference genome setup finished! \x1b[m"
echo -e "\x1b[92m Available genomes: \x1b[m"
ls -la genome/
