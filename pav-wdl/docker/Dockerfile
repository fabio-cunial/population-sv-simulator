FROM continuumio/miniconda3:4.11.0

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-pav/bin/:${PATH}

# test installation is correct
COPY test_install.py /tmp/
SHELL ["conda", "run", "-n", "lr-pav", "/bin/bash", "-c"]
RUN python /tmp/test_install.py

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.15.1
ARG BCFTOOLS_VERSION=1.15.1
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 make \
                 minimap2 \
                 tree \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends tree

ARG PAV_VERSION=1.2.1
RUN mkdir -p /opt/ && \
    git clone -b ${PAV_VERSION} --single-branch --recursive https://github.com/EichlerLab/pav.git /opt/pav
