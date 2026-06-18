##########  BUILDING THE IMAGE:  ###########################################################################################
#
#  This file will automatically
#  pull the codoe for the mint-sgmo-scripts and place them in a container
#  with the necessary applications and environment to run the scripts. You
#  will need the container to exist in an environment with access to an hg38 genome
#  in fasta format or, preferrably, a bowtie2 index of the same genome.
#
#  To build the image:
#  You must first configure the .env file for your local environment.  You will need to
#  determine where on your local machine you will store input data (fastq reads) and where you
#  want to write output data to.
#  If in a unix/linux environment, you will also want to determine your numerical UID
#  to ensure that permissions of files created are correct.
#
#  A sample file SAMPLE.env is included with the variables you will need to set
#  Set the proper local variables and save as '.env' in the docker_files directory.
#
#  With the .env file in the same directory as the Dockerfile and docker-compose.yml
#  issue the command:
#
#    $  docker compose build mint_environment
#
#  This may take a few minutes to construct the image but this only needs to be done once.
#  Once build, you can launch the image with the command.
#    $  docker compose run mint_environment
#
#  If successful, your container should launch and you will see a prompt similar to
#
#      root@3908bd0ca18a:/#
#
#  This indicates you are successfully operating in the container. The .env file should have
#  mounted the external data in and data out locations for access within the container.
#
#
#  All scripts are located within the directory /scripts
#  Please see the README for details on operating your scripts.

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

ARG BRANCH

RUN git clone --depth 1 https://github.com/sangamo-therapeutics/mint-sgmo-scripts.git /source_packages/mint-scripts  \
    --branch ${BRANCH} \
RUN cd /source_packages/mint-scripts && pip install ./

RUN mkdir -p /seqdb/bowtie2 /fastq_reads/aligned /figure_files/input /figure_files/output


ENV DEFAULT_FIGURE_OUT="/figure_files/output" \
    DEFAULT_FIGURE_IN="/figure_files/input" \
    PYTHONPATH='/source_pacakges/src/minter'


# # INSTALL THESE FROM pyproject.toml ??
#RUN pip install pandas matplotlib logomaker
#RUN pip install cutadapt python-levenshtein
#RUN pip install pybedtools
#RUN pip install logomaker scipy==1.13.1
#RUN pip install biopython==1.85

# sets local timezone
ENV TZ="America/Los_Angeles"
RUN date

# override the python entrypoint
ENTRYPOINT bash


# This point in build should be sufficient for building updates using internal sgmo zf-web applications
# referencing this as the base 'from' image