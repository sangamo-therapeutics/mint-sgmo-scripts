##########  BUILDING THE IMAGE:  ###########################################################################################
#

#########  BUILDING THE IMAGE:  ###########################################################################################
#
#  It is intended to create a base, cacheable point for faster rebuilding of the nonclin webapps image(s)
#    without modifing the underlyihg os or basic environment installations.
#
#

FROM python:3.11-slim-bookworm
USER root
RUN mkdir /configs /data /reference /src /logs /sequences /buildstuff /scripts

#RUN echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list

## add build tools and a some basic applications
RUN apt-get update \
  && apt-get -y install build-essential \
  && apt-get -y install zlib1g-dev liblzma-dev libbz2-dev libcurl4-gnutls-dev libssl-dev \
  && apt-get -y install git curl vim wget procps unzip \
  && apt-get clean

# python installations
RUN pip install -U pip==25.2
RUN pip install setuptools==57.0.0

# Install sequence 'deal with' apps
RUN apt-get update \
    && apt-get -y install samtools seqkit bedtools bowtie2 \
    && apt-get clean

#ENV BOWTIE2_VERSION=2.5.5
## 2.3.3.1
#WORKDIR /opt
#RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
#    &&  unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
#    && rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
#    && ln -s /opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2 /usr/local/bin \
#    && ln -s /opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2-build /usr/local/bin


RUN git clone --depth 1 https://github.com/sangamo-therapeutics/mint-sgmo-scripts.git /source_packages/mint-scripts  \
    --branch main
#    --branch ${BRANCH}

RUN mkdir -p /seqdb/bowtie2 /fastq_reads/aligned

RUN pip install pandas matplotlib logomaker
RUN pip install cutadapt python-levenshtein
RUN pip install pybedtools

# sets local timezone
ENV TZ="America/Los_Angeles"
RUN date

# override the python entrypoint
ENTRYPOINT bash


# This point in build should be sufficient for building updates using internal sgmo zf-web applications
# referencing this as the base 'from' image