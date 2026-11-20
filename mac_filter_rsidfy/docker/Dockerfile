FROM python:3.11-slim

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# System utilities required by the pipeline:
#   gawk  — annotation merge join + MAC/MAF filtering
#   gzip  — reading inputs
#   tabix — bgzip + tabix (htslib) for compressing and indexing outputs
RUN apt-get update && \
    apt-get install -y --no-install-recommends gawk gzip tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
