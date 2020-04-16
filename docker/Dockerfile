FROM debian:stretch

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libopenblas-dev \
    zip \
    python-dev \
    python3-dev \
    python3-pip \
    default-jre \
    git \
    curl

# Install all software under /opt/software:
RUN mkdir -p /opt/software 


RUN pip3 install regex cutadapt

RUN cd /opt/software && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip

# Install bbmap:
RUN cd /opt/software && \
    wget https://sourceforge.net/projects/bbmap/files/BBMap_38.81.tar.gz && \
    tar -xzf BBMap_38.81.tar.gz 

# Install samtools:
RUN cd /opt/software && \
  wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
  tar -xjf samtools-1.9.tar.bz2 && \
  cd samtools-1.9 && \
  ./configure --prefix=/opt/software/samtools && \
  make && \
  make install
ENV PATH="/opt/software/samtools/bin:${PATH}"

# Get the Picard JAR:
RUN mkdir -p /opt/software/picard && \
  cd /opt/software/picard && \
  wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
ENV PICARD_JAR=/opt/software/picard/picard.jar

# BEDtools:
RUN cd /opt/software && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
    tar -zxvf bedtools-2.29.1.tar.gz && \
    cd bedtools2 && \
    make

# BWA
RUN mkdir -p /opt/software/ && \
  cd /opt/software/ && \
  wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz && \
  tar xzf v0.7.17.tar.gz && \
  cd bwa-0.7.17 && \
  make
ENV PATH="/opt/software/bwa-0.7.17:${PATH}"

# install gcloud/gsutil:
RUN curl -sSL https://sdk.cloud.google.com > /tmp/gcl \
    && bash /tmp/gcl --install-dir=/root/gcloud --disable-prompts
ARG gcloud_bin=/root/gcloud/google-cloud-sdk/bin
ENV PATH "$PATH:"${gcloud_bin}
ENV GCLOUD=${gcloud_bin}/gcloud
RUN "${gcloud_bin}"/gcloud components install beta --quiet

COPY extract_readnames.py /opt/software/
COPY split_seqs.py /opt/software/

ENTRYPOINT ["/bin/bash"]
