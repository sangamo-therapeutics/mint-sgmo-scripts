# MINT scripts

This repository contains python scripts and configuration files to regenerate plots and analyses from Fauser et al.
Nature Biotechnology 2026.

The raw data files can be found in the NCBI Sequence Read Archive under BioProject accession number PRJNA1450726

The motif plots use logomaker

Scripts used to create figures shown in the manuscript are located in the directories:

    mint-sgmo-scripts/figure_maker_scripts
and
    mint-sgmo-scripts/bacterial-selection-specificity

The scripts are provided with a Dockerfile that has instructions on creating a container with 
the necessary applications for running the scripts.

#### BUILDING THE IMAGE:  ###########################################################################################

  This file will automatically
  pull the code for the mint-sgmo-scripts and place them in a container
  with the necessary applications and environment to run the scripts. You
  will need the container to exist in an environment with access to an hg38 genome
  in fasta format or, preferrably, a bowtie2 index of the same genome.

  To build the image:
  You must first configure the .env file for your local environment.  You will need to
  determine where on your local machine you will store input data (fastq reads) and where you
  want to write output data to.
  If in a unix/linux environment, you will also want to determine your numerical UID
  to ensure that permissions of files created are correct.

  A sample file SAMPLE.env is included with the variables you will need to set
  Set the proper local variables and save as '.env' in the docker_files directory.

  With the .env file in the same directory as the Dockerfile and docker-compose.yml
  issue the command:

    $  docker compose build mint_environment

  This may take a few minutes to construct the image but this only needs to be done once.
  Once build, you can launch the image with the command.
    $  docker compose run mint_environment

  If successful, your container should launch and you will see a prompt similar to

      root@3908bd0ca18a:/#

  This indicates you are successfully operating in the container. The .env file should have
  mounted the external data in and data out locations for access within the container.


  All scripts are located within the directories

    /bacterial-selection-specificity
  and
    /genome-wide-specificity-analysis


  Please see the README within each for details on operating your scripts.