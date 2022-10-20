#!/usr/bin/env python3

import sys
import re
import gzip

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_config(fname):
    eprint("reading config from " + fname)
    with open(fname, 'rt') as f:
        header = f.readline().strip().split()
        if len(header) != 6:
            raise Exception
        h_idx = {h: i for i,h in enumerate(header)}
        data = [line.strip().split() for line in f]
        config = [{h: d[h_idx[h]] for h in h_idx} for d in data]
    return config

def get_significant_vars(fname, pval_thresh):
    eprint("reading variants p<" + str(pval_thresh) + " from " + fname)
    vars = {}
    with gzip.open(fname, 'rt') as f:
        header = f.readline().strip().split()
        h_idx = {h: i for i,h in enumerate(header)}
        for line in f:
            s = re.split('\t| ', line.strip())
            if float(s[h_idx["pval"]]) < pval_thresh:
                v = s[0].replace("chr", "").replace("X", "23") + ':' + s[1] + ':' + s[2] + ':' + s[3]
                vars[v] = True
    return vars

def read_sumstat(vars, fname):
    eprint("reading subset of sumstats from " + fname)
    with gzip.open(fname, 'rt') as f:
        header = f.readline().strip().split()
        h_idx = {h: i for i,h in enumerate(header)}
        data = {}
        for line in f:
            s = re.split('\t| ', line.strip())
            v = s[0].replace("chr", "").replace("X", "23") + ':' + s[1] + ':' + s[2] + ':' + s[3]
            if v in vars:
                datum = {h: s[h_idx[h]] for h in h_idx}
                data[v] = datum
    return data

def print_joined_data(uniq_vars, data, config, columns):
    header = "chr\tpos\tref\talt"
    for c in config:
        for col in columns:
            header += "\t" + c["tag"] + '_' + c[col]
    print(header)
    for v in uniq_vars:
        line = "\t".join(v.split(":"))
        for i,sumstat_subset in enumerate(data):
            if v in sumstat_subset:
                line += "\t" + "\t".join([sumstat_subset[v][config[i][col]] for col in columns])
            else:
                line += "\t" + "\t".join(["NA" for col in columns])
        print(line)

if __name__ == '__main__':
    config = read_config(sys.argv[1])
    filtered_vars = [get_significant_vars(c["file"], float(c["pval_thresh"])) for c in config]
    uniq_vars = sorted(list(set([v for vars in filtered_vars for v in vars])))
    data = [read_sumstat(set(uniq_vars), c["file"]) for c in config]
    print_joined_data(uniq_vars, data, config, ["pval_col", "beta_col", "af_col"])
    eprint("done")
