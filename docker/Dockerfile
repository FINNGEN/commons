FROM ubuntu:16.04
ARG SSH_PRIVATE_KEY
RUN apt-get update && apt-get install mercurial --yes && \
apt-get install git lsb-release gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev \
bzip2 curl python make gcc g++ libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes && \
wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
tar -xvjf htslib-1.7.tar.bz2 && cd htslib-1.7 && ./configure && make && make install && \
hg clone -b beta https://gavinband@bitbucket.org/gavinband/qctool && \
cd qctool && ./waf-1.5.18 configure && ./waf-1.5.18 && \
hg clone https://bitbucket.org/gavinband/bgen && \
cd bgen && ./waf-1.8.13 configure && ./waf-1.8.13 &&
echo "${SSH_PRIVATE_KEY}" >> /root/.ssh/id_rsa && chmod 600 /root/.ssh/id_rsa  && \
git clone git clone git@github.com:FINNGEN/CONVERT_AND_SPLIT_VCF.git


FROM ubuntu:16.04

RUN apt-get update && apt-get install curl lsb_release python liblzma-dev zlib1g-dev libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes && \\
&& COPY --from=0 /bgen/build/apps/* /usr/local/bin && COPY --from=0 /qctool/build/release/qctool_v2.0-rc8 /usr/local/bin && \
COPY --from=0 /qctool/build/release/tabix /usr/local/bin && COPY --from=0 /qctool/build/release/bgzip /usr/local/bin && \
COPY --from=0 /qctool/build/release/htsfile /usr/local/bin && \\
COPY --from=0 /CONVERT_AND_SPLIT_VCF/scripts/split_single_chrom.py /usr/local/bin
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \\
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add - && \
apt-get update && apt-get install google-cloud-sdk --yes && \
