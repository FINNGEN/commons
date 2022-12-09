FROM ubuntu:22.04

LABEL name="BIOINFORMATICS"
LABEL maintainer="jkarjala@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

ARG HTSLIB_VER=1.16
ARG BCFTOOLS_VER=1.16
ARG R_VER=4.2.2

WORKDIR /tmp

RUN apt-get update && apt-get upgrade --yes

RUN apt-get install vim emacs jq less datamash wget gawk curl lsb-core make gcc g++ g++-9 cmake gfortran libreadline6-dev \
    libz-dev libbz2-dev liblzma-dev zlib1g-dev libpcre2-dev libpcre3-dev libssl-dev libcurl4-openssl-dev \
    libopenblas-dev default-jre unzip bzip2 libboost-all-dev \
    libpng-dev libcairo2-dev python3-pip libxml2-dev --yes

RUN curl -LO https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 && \
    tar -xvjf htslib-${HTSLIB_VER}.tar.bz2 && cd htslib-${HTSLIB_VER} && \
    ./configure && make && make install && cd .. && rm -rf htslib-${HTSLIB_VER}*

RUN wget https://code.enkre.net/qctool/zip/release/qctool.tgz && \
    unzip qctool.tgz && \
    cd qctool && python3 waf configure && python3 waf && \
    mv build/release/apps/qctool_v2* /usr/local/bin/qctool && cd .. && rm -rf qctool*

# bgen tools don't compile with g++ 11 so use 9
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && \
    tar xzf bgen.tgz && cd bgen.tgz && \
    CXX=/usr/bin/g++-9 python3 waf configure && python3 waf && \
    mv build/apps/bgenix build/apps/cat-bgen build/apps/edit-bgen /usr/local/bin/ && cd .. && rm -rf bgen.tgz

RUN curl -LO https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2 && \
    tar -xvjf bcftools-${BCFTOOLS_VER}.tar.bz2 && cd bcftools-${BCFTOOLS_VER} && \
    ./configure && make && make install && cd .. && rm -rf bcftools-*

RUN mkdir plink && cd plink && curl -O http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && \
    unzip plink_linux_x86_64_latest.zip && mv plink /usr/local/bin && rm * && \
    curl -O http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip && \
    unzip plink2_linux_avx2_latest.zip && mv plink2 /usr/local/bin && cd .. && rm -rf plink

RUN curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod a+x liftOver && mv liftOver /usr/local/bin && \
    curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz && \
    curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz && \
    mkdir /liftover && mv *.chain.gz /liftover/

ENV CLOUDSDK_INSTALL_DIR /usr/local/gcloud/
RUN curl -sSL https://sdk.cloud.google.com | bash
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin

RUN curl -O https://cloud.r-project.org/src/base/R-4/R-${R_VER}.tar.gz && \
    tar xvzf R-${R_VER}.tar.gz && \
    cd R-${R_VER} && \
    ./configure --with-x=no --enable-R-shlib=yes --with-blas="-lopenblas" && \
    make && mkdir -p /usr/local/lib/R/lib && make install && cd .. && rm -rf R-${R_VER}*

ADD install_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install_packages.R && \
    install_packages.R

ADD requirements.txt .
RUN pip3 install -r requirements.txt && rm requirements.txt

RUN apt-get remove default-jre --yes && apt autoremove --yes && apt-get clean

ADD print_versions.sh print_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/print_versions.sh /usr/local/bin/print_packages.R && \
    print_versions.sh

WORKDIR /
