FROM ubuntu:16.04

RUN apt-get update && apt-get install mercurial --yes && \
apt-get install curl python make gcc g++ libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes && \
hg clone -b beta https://gavinband@bitbucket.org/gavinband/qctool && \
cd qctool && ./waf-1.5.18 configure && ./waf-1.5.18 && \
ENV PATH $PWD/build/release/:$PATH && \
hg clone https://bitbucket.org/gavinband/bgen && \
cd bgen && ./waf-1.8.13 configure && ./waf-1.8.13 && \
ENV PATH /bgen/build/apps/:$PATH

FROM ubuntu:16.04
RUN apt-get update && apt-get install curl python libz-dev libbz2-dev liblzma-dev libreadline6 libreadline6-dev --yes && \\
&& COPY --from=0 /bgen/build/apps/* /usr/local/bin && COPY --from=0 /qctool/build/release/qctool_v2.0-rc8 /usr/local/bin
