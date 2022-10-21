FROM continuumio/miniconda3
MAINTAINER Fabio Cunial
ARG work_dir=/simulation
WORKDIR ${work_dir}


# --------------------------------- Versions -----------------------------------
ARG gcloud_version=405.0.1
ARG samtools_version=1.16.1
ARG bcftools_version=1.16
ARG minimap2_version=2.24
ARG pbsv_version=2.8.0
ARG sniffles1_version=1.0.12
ARG sniffles2_version=2.0.7
ARG hifiasm_version=0.16.1
ARG pav_version=2.1.0
# PBSIM1 has no version
# ------------------------------------------------------------------------------


# Basics
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        apt-transport-https \
        ca-certificates \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        default-jdk \
        default-jre \
        git \
        python3-pip

# GSUTIL
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}
        
# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# Bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help

# Minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v${minimap2_version}/minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && tar xjf minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && rm minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && cd minimap2-${minimap2_version}_x64-linux \
    && cp minimap2 /usr/local/bin/ \
    && cd ${work_dir} \
    && rm -rf minimap2-${minimap2_version}_x64-linux \
    && minimap2 --help

# PBSV
RUN conda update conda \
    && conda install -c bioconda pbsv=${pbsv_version} \
    && pbsv --help

# Sniffles 1
RUN wget https://github.com/fritzsedlazeck/Sniffles/archive/refs/tags/v${sniffles1_version}a.tar.gz \
    && tar -xzf v${sniffles1_version}a.tar.gz \
    && rm -f v${sniffles1_version}a.tar.gz \
    && cd Sniffles-${sniffles1_version}a \
    && mkdir -p build/ \
    && cd build/ \
    && cmake .. \
    && make \
    && cp ../bin/sniffles-core-${sniffles1_version}/sniffles /usr/local/bin/sniffles1 \
    && cd ${work_dir} \
    && rm -rf Sniffles-${sniffles1_version}a \
    && sniffles1 --help
    
# Sniffles 2
RUN conda update conda \
    && conda install -c bioconda -c conda-forge -c defaults sniffles=${sniffles2_version} \
    && sniffles --help

# Hifiasm
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${hifiasm_version}.tar.gz \
    && tar -xzf ${hifiasm_version}.tar.gz \
    && rm -f ${hifiasm_version}.tar.gz \
    && cd hifiasm-${hifiasm_version} \
    && make \
    && cp hifiasm /usr/local/bin/ \
    && cd ${work_dir} \
    && rm -rf hifiasm-${hifiasm_version} \
    && hifiasm -h

# PAV
COPY ./pav-environment.yml .
RUN conda env create -f ./pav-environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-pav/bin/:${PATH}
COPY pav_test_install.py .
SHELL ["conda", "run", "-n", "lr-pav", "/bin/bash", "-c"]
RUN echo "conda init bash" > ~/.bashrc
RUN python pav_test_install.py
RUN wget https://github.com/EichlerLab/pav/archive/refs/tags/v${pav_version}.tar.gz \
    && tar -xzf v${pav_version}.tar.gz \
    && rm -f v${pav_version}.tar.gz \
    && mv pav-${pav_version}/ pav/

# PBSIM 1
#
# Remark: we use PBSIM 1 in CLR mode, since that is the only way we can
# customize the error rate. PBSIM2 dropped support for CCS.
#
# Remark: we do not use BADREAD <https://github.com/rrwick/Badread> since it is
# too slow in practice.
RUN git clone https://github.com/pfaucon/PBSIM-PacBio-Simulator
# File <src/pbsim.cpp> was modified to output only a FASTA file. The original
# program outputs a FASTQ and a MAF file that contains read alignments to the
# reference: these take a lot of disk space and are useless for our purposes.
COPY pbsim.cpp PBSIM-PacBio-Simulator/src/
RUN cd PBSIM-PacBio-Simulator \
    && autoreconf -i \
    && ./configure \
    && make \
    && cp src/pbsim /usr/local/bin \
    && cp -r data/ ${work_dir}/pbsim_data \
    && cd ${work_dir} \
    && rm -rf PBSIM-PacBio-Simulator \
    && pbsim --help

# Simulation code
COPY haplotype2reads.sh .
COPY reads2svs.sh .
COPY pav.sh .
COPY PrintHaplotypes.java .
COPY SimulateHaplotypes.java .
COPY BuildModel.java .
RUN javac *.java
