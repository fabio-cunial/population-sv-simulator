FROM continuumio/miniconda3
MAINTAINER Fabio Cunial
ARG work_dir=/simulation
WORKDIR ${work_dir}


# --------------------------------- Versions -----------------------------------
ARG gcloud_version=405.0.1
ARG htslib_version=1.16
ARG samtools_version=1.16.1
ARG bcftools_version=1.16
ARG minimap2_version=2.24
ARG pbsv_version=2.8.0
ARG sniffles1_version=1.0.12
ARG sniffles2_version=2.0.7
ARG hifiasm_version=0.16.1
ARG quast_version=5.2.0
# - PAV has releases, but for simplicity we just clone the repo to resolve a 
# dependency on SVPOP. Should be done more rigorously.
# - PBSIM2 and PBSIM1 have no version
# ------------------------------------------------------------------------------


# OS
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
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
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

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
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

# BCFTOOLS
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

# MINIMAP2
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
    && conda clean --packages --tempfiles --logfiles \
    && pbsv --help

# SNIFFLES1
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
    
# SNIFFLES2
RUN conda update conda \
    && conda install -c bioconda -c conda-forge -c defaults sniffles=${sniffles2_version} \
    && conda clean --packages --tempfiles --logfiles \
    && sniffles --help

# HIFIASM
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${hifiasm_version}.tar.gz \
    && tar -xzf ${hifiasm_version}.tar.gz \
    && rm -f ${hifiasm_version}.tar.gz \
    && cd hifiasm-${hifiasm_version} \
    && make \
    && cp hifiasm /usr/local/bin/ \
    && cd ${work_dir} \
    && rm -rf hifiasm-${hifiasm_version} \
    && hifiasm -h

# QUAST
RUN conda update conda \
    && conda install -c bioconda -c conda-forge libgcc-ng joblib quast \
    && conda clean --packages --tempfiles --logfiles \
    && quast --version \
    && python -c 'import joblib'

# PAV
COPY ./pav_environment.yml .
RUN conda env create -f ./pav_environment.yml \
    && conda clean --packages --tempfiles --logfiles
COPY pav_test_install.py .
SHELL ["conda", "run", "-n", "lr-pav", "/bin/bash", "-c"]
RUN echo "conda init bash" > ~/.bashrc
RUN python pav_test_install.py
RUN git clone --recursive https://github.com/EichlerLab/pav.git
# We modified <pav/rules/call.snakefile> as follows:
# - We set MERGE_BATCH_COUNT=1, since this parameter is not accessible from the
#   configuration JSON.
# - We fixed an empty CSV bug in rule <call_merge_haplotypes>.
COPY pav_call.snakefile ./pav/rules/call.snakefile
# We fixed an empty BED bug in rule <vcf_write_vcf> of file
# <pav/rules/vcf.snakefile>.
COPY pav_vcf.snakefile ./pav/rules/vcf.snakefile

# PBSIM2
#
# Remark: PBSIM 2 is a CLR simulator (it dropped support for CCS), but the
# error rate can be customized. PBSIM 1 supported CCS, but it was probably an
# old model that is unrealistic now. PBSIM 3 seems to be under development
# <https://github.com/yukiteruono/pbsim3>: it seems to simulate CLR and to feed
# it to the CCS consensus algorithm. The paper does not seem to be accepted yet.
#
# Remark: we do not use BADREAD <https://github.com/rrwick/Badread> since it is
# too slow in practice.
RUN git clone https://github.com/yukiteruono/pbsim2.git
# We modified file <src/pbsim.cpp> to output only a FASTA file. The original
# program outputs a FASTQ and a MAF file that contains read alignments to the
# reference: these take a lot of disk space and are useless for our purposes.
COPY pbsim2.cpp pbsim2/src/pbsim.cpp
RUN cd pbsim2 \
    && ./configure \
    && make \
    && mv src/pbsim /usr/local/bin \
    && mv data/ ${work_dir}/pbsim_data \
    && make clean \
    && pbsim
# We also include a version of <src/pbsim.cpp> that outputs FASTA without
# errors. Clearly we could set errors to zero via command line arguments; this
# version prints the substring sampled from the reference, just to be sure.
COPY pbsim2noError.cpp pbsim2/src/pbsim.cpp
RUN cd pbsim2 \
    && make \
    && mv src/pbsim /usr/local/bin/pbsimNoError \
    && cd ${work_dir} \
    && rm -rf pbsim2 \
    && pbsimNoError

# PBSIM1
# RUN git clone https://github.com/pfaucon/PBSIM-PacBio-Simulator
# COPY pbsim1.cpp PBSIM-PacBio-Simulator/src/pbsim.cpp
# RUN cd PBSIM-PacBio-Simulator \
#     && autoreconf -i \
#     && ./configure \
#     && make \
#     && mv src/pbsim /usr/local/bin \
#     && mv data/ ${work_dir}/pbsim_data \
#     && cd ${work_dir} \
#     && rm -rf PBSIM-PacBio-Simulator \
#     && pbsim --help

# Simulation code
COPY haplotype2reads.sh .
COPY reads2svs.sh .
COPY *.java .
RUN javac *.java
