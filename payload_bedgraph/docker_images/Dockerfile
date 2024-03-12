FROM ubuntu:20.04

# Set paths.

ARG PATH=/root/miniconda3/bin:$PATH
ENV PATH=/root/miniconda3/bin:$PATH
ENV PATH=/root/miniconda3/envs/pl15s/bin:$PATH

ENV DEBIAN_FRONTEND=noninteractive

RUN \
    apt-get -y update \
    && apt-get -y install git wget vim libtiff5 libjpeg9 libarchive13 tzdata \
    && apt-get -y install build-essential
ENV TZ="America/New_York"

RUN \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && conda update -y -n base -c defaults conda \
    && conda install -y -n base -c conda-forge mamba \
    && mamba update -y -n base -c conda-forge mamba

RUN \
    mamba create -y -n pl15s -c conda-forge -c bioconda r-base=4.1.2 r-rmarkdown==2.21  \
    && mamba clean -y --all

RUN R -e 'install.packages("BiocManager", repos="https://cloud.r-project.org/")'
RUN R -e 'BiocManager::install(version = "3.14")'
RUN R -e 'BiocManager::install(c("GenomicRanges", "DT", "tidyverse", "gridextra"))'
RUN R -e 'install.packages(c("biomaRt", "XML"), repos="https://cran.r-project.org")'
RUN R -e 'BiocManager::install(c("Sushi"))'

# Install wrapper script
RUN \
    mkdir git_tmp \
    && git clone https://github.com/egenesis/seven_bridges_apps.git git_tmp/ \
    && cd git_tmp/payload_bedgraph/ \
    && git checkout d192f843a20c0702d86f5d289028a29b55c43aa6 \
    && chmod a+x coverage_bedgraph.Rmd \
    && mv coverage_bedgraph.Rmd /opt/ \
    && chmod a+x run_report.R \
    && mv run_report.R /opt/

CMD ["/bin/bash"]
