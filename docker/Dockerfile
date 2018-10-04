FROM ubuntu:16.04
COPY scripts/lift.py scripts/split_single_chrom.py scripts/scrape_vcf_info.py scripts/liftOver data/hg19ToHg38.over.chain.gz ./
RUN apt-get update && apt-get install mercurial --yes && \
apt-get install gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev \
bzip2 curl python3 make gcc g++ libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes && \
wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
tar -xvjf htslib-1.7.tar.bz2 && cd htslib-1.7 && ./configure && make && make install && cd .. && \
hg clone -b beta https://gavinband@bitbucket.org/gavinband/qctool && \
cd qctool && ./waf-1.5.18 configure && ./waf-1.5.18 && cd ..  && \
hg clone https://bitbucket.org/gavinband/bgen && \
cd bgen && ./waf configure && ./waf && \
wget https://github.com/samtools/bcftools/releases/download/1.7/bcftools-1.7.tar.bz2 && tar -xvjf bcftools-1.7.tar.bz2 && cd bcftools-1.7 && \
./configure && make && make install

FROM ubuntu:16.04
RUN apt-get update && apt-get install wget gawk curl lsb-core libcurl4-openssl-dev python liblzma-dev zlib1g-dev libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes
#RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && apt-get update && apt-get install google-cloud-sdk --yes
# installing gcloud over doubles the image size. install later if needed

COPY --from=0 /usr/local/bin/tabix /usr/local/bin/
COPY --from=0 /usr/local/bin/bgzip /usr/local/bin/
COPY --from=0 /usr/local/bin/htsfile /usr/local/bin/
COPY --from=0 /usr/local/bin/bcftools /usr/local/bin/
COPY --from=0 split_single_chrom.py /usr/local/bin/
COPY --from=0 lift.py /usr/local/bin/
COPY --from=0 scrape_vcf_info.py /usr/local/bin/
COPY --from=0 /qctool/build/release/qctool_v2* /usr/local/bin/qctool
COPY --from=0 /bgen/build/apps/bgenix /usr/local/bin/
COPY --from=0 /bgen/build/apps/cat-bgen /usr/local/bin/
COPY --from=0 liftOver /usr/local/bin/
COPY --from=0 hg19ToHg38.over.chain.gz ./
