FROM ubuntu:16.04
ADD scripts/split_single_chrom.py split_single_chrom.py
RUN apt-get update && apt-get install mercurial --yes && \
apt-get install gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev \
bzip2 curl python make gcc g++ libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes && \
wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
tar -xvjf htslib-1.7.tar.bz2 && cd htslib-1.7 && ./configure && make && make install && \
hg clone -b beta https://gavinband@bitbucket.org/gavinband/qctool && \
cd qctool && ./waf-1.5.18 configure && ./waf-1.5.18 && \
hg clone https://bitbucket.org/gavinband/bgen && \
cd bgen && ./waf-1.8.13 configure && ./waf-1.8.13

FROM ubuntu:16.04
RUN apt-get update && apt-get install curl libcurl4-openssl-dev python liblzma-dev zlib1g-dev libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes 
COPY --from=0 /usr/local/bin/bgenix /usr/local/bin 
COPY --from=0 /usr/local/bin/cat-bgen /usr/local/bin
COPY --from=0 /usr/local/bin/qctool_v2.0-rc8 /usr/local/bin
COPY --from=0 /usr/local/bin/tabix /usr/local/bin 
COPY --from=0 /usr/local/bin/bgzip /usr/local/bin/
COPY --from=0 /usr/local/bin/htsfile /usr/local/bin
COPY --from=0 split_single_chrom.py /usr/local/bin
