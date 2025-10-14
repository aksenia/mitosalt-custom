FROM python:3.12-bullseye

ENV DEBIAN_FRONTEND=noninteractive \
    LANGUAGE=C.UTF-8 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8


# Install system dependencies and tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre wget unzip ca-certificates \
    python3 python3-pip python-is-python3 python3-venv \
    perl r-base r-cran-plotrix r-cran-rcolorbrewer \
    hisat2 last-align samtools make gcc g++ zlib1g-dev \
    libbz2-dev liblzma-dev expect \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install bedtools manually (latest stable release)
ENV BEDTOOLS_VERSION=2.31.1
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    tar -zxvf bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/bedtools /bin/bedtools && \
    cd / && \
    rm -rf bedtools2 bedtools-${BEDTOOLS_VERSION}.tar.gz

# Install sambamba (latest version)
ENV SAMBAMBA_VERSION=1.0.1
RUN wget -O /tmp/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz && \
    gunzip /tmp/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz && \
    mv /tmp/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static /usr/bin/sambamba && \
    chmod +x /usr/bin/sambamba

# Install BBMap (latest version)
ENV BBMAP_VERSION=38.06
RUN wget -qO /tmp/BBMap.tar.gz "https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download" && \
    mkdir -p /usr/bin/bbmap && \
    tar -xzf /tmp/BBMap.tar.gz -C /usr/bin/bbmap --strip-components=1 && \
    rm /tmp/BBMap.tar.gz
ENV PATH="/usr/bin/bbmap:${PATH}"

# Install UCSC utilities (v369) into /bin as expected by MitoSAlt
ENV UCSC_VERSION=v369
RUN wget -O /bin/bedGraphToBigWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.${UCSC_VERSION}/bedGraphToBigWig && \
    chmod +x /bin/bedGraphToBigWig && \
    wget -O /bin/faSomeRecords http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.${UCSC_VERSION}/faSomeRecords && \
    chmod +x /bin/faSomeRecords && \
    wget -O /bin/faSize http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.${UCSC_VERSION}/faSize && \
    chmod +x /bin/faSize

# Dependencies for R visualisations
RUN R -e "install.packages(c('plotrix', 'RColorBrewer'), repos='https://cloud.r-project.org/')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')"

# Download and unpack MitoSAlt
WORKDIR /opt

ENV MITOSALT_VERSION=1.1.1
RUN wget -O mitosalt.zip "https://sourceforge.net/projects/mitosalt/files/MitoSAlt_${MITOSALT_VERSION}.zip/download" && \
    unzip mitosalt.zip && \
    rm mitosalt.zip

WORKDIR /opt/MitoSAlt_${MITOSALT_VERSION}

# Correct the paths in the analysis and plotting R script to work with container
RUN sed -i \
    -e "s|plotfile<-paste(\"plot/\",filename,\".pdf\",sep=\"\")|plotfile<-paste(\"/output/plot/\",filename,\".pdf\",sep=\"\")|" \
    -e "s|textfile<-paste(\"indel/\",filename,\".tsv\",sep=\"\")|textfile<-paste(\"/output/indel/\",filename,\".tsv\",sep=\"\")|" \
    /opt/MitoSAlt_${MITOSALT_VERSION}/delplot.R

RUN sed -i "s#delplot\.R#/opt/MitoSAlt_${MITOSALT_VERSION}/delplot.R#g" MitoSAlt${MITOSALT_VERSION}.pl

# Copy saltshaker package source code
COPY saltshaker/ /opt/MitoSAlt_${MITOSALT_VERSION}/saltshaker/
COPY setup.py pyproject.toml README.md /opt/MitoSAlt_${MITOSALT_VERSION}/

# Copy other MitoSAlt scripts and config
COPY config_human.txt download_genomes.sh /opt/MitoSAlt_${MITOSALT_VERSION}/

# Install the saltshaker package
RUN pip install --no-cache-dir --break-system-packages .

# Make helper scripts executable
RUN chmod +x download_genomes.sh


RUN echo '#!/usr/bin/expect -f\n\
set timeout -1\n\
spawn ./setup.sh\n\
expect "Do you want basic installation*" { send "no\\r" }\n\
expect "Erase existing files*" { send "no\\r" }\n\
expect "Download and compile hisat2*" { send "no\\r" }\n\
expect "Download and compile bedtools2*" { send "no\\r" }\n\
expect "Download and compile bbmap*" { send "no\\r" }\n\
expect "Download and compile last*" { send "no\\r" }\n\
expect "Download and compile samtools*" { send "no\\r" }\n\
expect "Download and compile sambamba*" { send "no\\r" }\n\
expect "Get UCSC utilities*" { send "no\\r" }\n\
expect "Download human genome*" { send "no\\r" }\n\
expect "Build human indexes*" { send "no\\r" }\n\
expect "Download mouse genome*" { send "no\\r" }\n\
expect "Build mouse indexes*" { send "no\\r" }\n\
expect "Install R package Plotrix*" { send "no\\r" }\n\
expect "Install R package RColorBrewer*" { send "no\\r" }\n\
expect "Number of threads*" { send "4\\r" }\n\
expect eof' > mitosalt_setup.exp && chmod +x mitosalt_setup.exp

RUN ./mitosalt_setup.exp

# Add MitoSAlt and tools to PATH
ENV PATH="/opt/MitoSAlt_${MITOSALT_VERSION}:/usr/bin:/usr/local/bin:/bin:${PATH}"

# Set working directory so imports work
# WORKDIR /data

# Ensure Python path includes the directory
ENV PYTHONPATH="/opt/MitoSAlt_${MITOSALT_VERSION}:${PYTHONPATH}"
CMD ["/bin/bash"]
