import sys
import gzip
import pandas as pd
from decimal import Decimal

filename=sys.argv[1]
p_threshold = float(sys.argv[2])

print('\t'.join(['#variant', 'pheno', 'pval', 'beta', 'sebeta', 'maf', 'maf_cases', 'maf_controls']))
with gzip.open(filename, mode='rt', encoding='utf8') as f:
    headers = f.readline().strip().split('\t')
    pval_h = list(filter(lambda x: x is not None, [i if h.startswith('pval') else None for i, h in enumerate(headers)]))
    for line in f:
        vals = line.strip().split('\t')
        for i in pval_h:
            if vals[i] != 'NA':
                pval = float(vals[i])
                if pval < p_threshold:
                    result = [headers[i][5:], '{0:.2E}'.format(Decimal(pval)), vals[i+1], vals[i+2], vals[i+3], vals[i+4], vals[i+5]]
                    print(vals[0].replace('X', '23') + ':' + vals[1] + ':' + vals[2] + ':' + vals[3] + '\t' + '\t'.join(result))

