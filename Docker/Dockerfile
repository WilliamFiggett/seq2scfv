# building in ubuntu 20.04
#FROM ubuntu@sha256:b795f8e0caaaacad9859a9a38fe1c78154f8301fdaf0872eaf1520d66d9c0b98 AS builder

# Updated base image to fix compatability with R packages:
FROM rocker/tidyverse:latest

SHELL ["/bin/bash", "-c"]

ENV CONTAINER_NAME="NGS-AI Seq2scFv"
ENV CONTAINER_VERSION="0.1"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/usr/local/bin:$PATH"
ENV PATH="/usr/miniconda/bin:$PATH"

# turn off interactive program installation
ARG DEBIAN_FRONTEND=noninteractive
# version numbers
ARG VERSION_SEQKIT="2.3.1"
ARG VERSION_IGBLAST="1.20.0"
#ARG VERSION_PANDAS="2.2.2"
#ARG VERSION_BIOPYTHON="1.84"
ARG VERSION_CLUSTALO="1.2.4"
ARG VERSION_BIOAWK="1.0"

# install some necessary base tools for construction
RUN apt-get update -y && apt-get install -y unzip wget perl python3 make build-essential libbz2-dev zlib1g-dev liblzma-dev lzma liblzma-dev git libxml2 r-base

# install conda
RUN echo "Install Bioconda"
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
RUN bash miniconda3.sh -b -u -p /usr/miniconda
ENV PATH="/usr/miniconda/bin:$PATH"
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# accept channel tos
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# install python packages
RUN conda install -y pandas biopython

# conda installation of other tools
RUN conda install -y -c bioconda seqkit=${VERSION_SEQKIT} \
                                 clustalo=${VERSION_CLUSTALO} \
                                 bioawk=${VERSION_BIOAWK}

RUN pip3 install logomaker \
                 antpack

# install IgBLAST
ARG PATH_IGBLAST="/igblast"
ENV IGDATA=${PATH_IGBLAST}/ncbi-igblast-${VERSION_IGBLAST}
WORKDIR ${PATH_IGBLAST}
RUN wget --no-check-certificate https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${VERSION_IGBLAST}/ncbi-igblast-1.20.0-x64-linux.tar.gz &&\
    tar zxf ncbi-igblast-${VERSION_IGBLAST}-x64-linux.tar.gz &&\
    ln -t /usr/local/bin ncbi-igblast-${VERSION_IGBLAST}/bin/* &&\
    rm ncbi-igblast-${VERSION_IGBLAST}-x64-linux.tar.gz

COPY install2.r /usr/bin/install2.r
RUN chmod +x /usr/bin/install2.r

RUN mkdir /home/seq2scfv
# Use this forked repository (dev branch) to include script corrections
RUN git clone -b dev https://github.com/WilliamFiggett/seq2scfv.git /home/seq2scfv
RUN cp /home/seq2scfv/bin/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*

RUN mkdir /home/R
ENV R_LIBS_USER="/home/R"

# Install R packages
RUN install2.r --error \
    data.table \
    reshape2 \
    dplyr \
    tidyr \
    ggplot2 \
    gridExtra \
    stringr

CMD ["/bin/bash"]
