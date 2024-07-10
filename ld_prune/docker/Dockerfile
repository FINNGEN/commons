FROM eu.gcr.io/phewas-development/tomahawk:beta-0.7.1-dirty-fg-v1-ubuntu22
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install gcsfuse
RUN export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s` && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update --fix-missing && \
    apt-get install -y wget gawk make zlib1g-dev g++ libopenblas-dev gfortran bedtools python-pip gcsfuse && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get -y  install apt-transport-https ca-certificates gnupg curl &&  \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" |\
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && apt-get install -y  google-cloud-cli

# HTSLIB (for bgzip)
ENV HTSLIB_VERSION 1.10.2
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --disable-bz2 --disable-lzma && \
    make && make install && \
    rm -r ~/htslib-${HTSLIB_VERSION}*


RUN mkdir /scripts 
COPY scripts /scripts
ENV PYTHONPATH=/scripts
ENV PATH="$PATH:/scripts"


ADD requirements.txt requirements.txt

RUN wget https://bootstrap.pypa.io/get-pip.py && python3 get-pip.py && pip install -r requirements.txt

