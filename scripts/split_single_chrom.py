#!/usr/bin/env python3
import argparse

from multiprocessing import Pool, TimeoutError
import multiprocessing


import subprocess as sb
from functools import partial
import time

QCTOOL="qctool_v2.0-rc8"
OFILETYPE=""

def run_chunk_vcf(inputfile, chr, start, end, output, bits, ofiletype, rounding_error):

    index_cmd = ['tabix', inputfile, "{}:{}-{}".format(chr, start, end)]

    qct_cmd = [QCTOOL, '-g -', '-vcf-genotype-field',' GP', '-filetype vcf', '-og', output, '-bgen-compression zlib',
               '-ofiletype', ofiletype, '-bgen-input-rounding-error', "{}".format(rounding_error)  ]

    #print("running pipe: {} | {}".format( " ".join(index_cmd), " ".join(qct_cmd) ))

    tab = sb.Popen(index_cmd,stdout=subprocess.PIPE)
    tp2 = subprocess.Popen(qct_cmd, stdin=tab.stdout)
    p2.communicate()
    return 0

def run_chunk_bgen(inputfile, chr, start, end, output):

    index_cmd = ['bgenix', '-g', '-incl-range', "{}:{}-{}".format(chr, start, end)]

    print("running command {}".format(index_cmd) )

    with open(output, 'w') as out:
        index_cmd = ['bgenix', '-g', '-incl-range', "{}:{}-{}".format(chr, start, stop)]
        sb.check_call(index_cmd, stdout=out)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="split a chromosome to predefined bgen chunks")

    parser.add_argument('inputfile', action='store', type=str,
                        help='vcf or bgen file to splot')
    parser.add_argument('chunks', action='store', type=str,
                        help='chunks as chr start stop')

    parser.add_argument('-nthreads', action='store', type=int,
                        help='how many threads to run.')

    parser.add_argument('-bgentype', default="bgen_v1.2", action='store', type=str,
                            help='how many threads to run.')

    parser.add_argument('-bgenbits', default=8, action='store', type=int,
                        help='how many threads to run.')

    parser.add_argument('-bgenrounderror', default=0.005, action='store', type=float,
                        help='how many threads to run.')

    args = parser.parse_args()

    ncpus = 1
    if ( args.nthreads is not None):
        ncpus = nthread
    else:
        ncpus = multiprocessing.cpu_count() - 1


    print("Running with {} cpus".format(ncpus) )
    threadpool = Pool(processes=ncpus )

    conv_func = None

    if (args.inputfile.endswith("bgen")):
        conv_func = partial(run_chunk_bgen, inputfile=args.inputfile)
    elif (args.inputfile.endswith("vcf.bgz") or args.inputfile.endswith("vcf.bg") or args.inputfile.endswith("vcf.gz")):
        conv_func = partial(run_chunk_vcf,inputfile=args.inputfile, bits=args.bgenbits, ofiletype=args.bgentype, rounding_error=args.bgenrounderror)
    else:
        raise Exception("Unsupported inputfile")

    n_processed = 0
    ## thread unsafety not critical
    def progress(val):
        global n_processed
        n_processed += 1
        if( n_processed % 10==0):
            print("Progress: {} chunks".format( n_processed))

    def error_callback(val):
        print("ERROR {}".format(val))

    with open( args.chunks, 'r') as chunks:
        for line in chunks:
            dat = line.rstrip("\n").split()
            outfile= "{}_{}_{}_{}.bgen".format(args.inputfile, dat[0], dat[1], dat[2] )
            threadpool.apply_async( partial(conv_func, chr=dat[0], start=dat[1], end=dat[2], output=outfile), callback=progress )

    threadpool.close()
    threadpool.join()