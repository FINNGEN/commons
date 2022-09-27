## Join subsets of sumstats

This is a small python3 program that joins summary statistic text files on a subset of variants for use in e.g. the multiple Manhattan plot tool. 

Usage:
```
python3 mmss.py config.txt > joined_sumstats.tsv
```

Configuration is given in a text file like:

```
$ cat config.txt | column -t
file                       tag      pval_col  beta_col  af_col  pval_thresh
data/WETAMDSTRICT.gz       wet      pval      beta      af_alt  1e-6
data/DRYAMDSTRICTV2.gz     dry      pval      beta      af_alt  1e-6
data/finngen_R9_H7_AMD.gz  both_fg  pval      beta      af_alt  1e-6
```

The program will read all variants below the given significance threshold in each configured summary stat file.
It will then take the union of these variants and print out stats from each summary stat file for the union of the variants
(NA will be printed if a variant is missing in a summary stat):

```
chr  pos        ref  alt  wet_pval     wet_beta  wet_af_alt  dry_pval     dry_beta   dry_af_alt  both_fg_pval  both_fg_beta  both_fg_af_alt
10   121340977  G    A    0.00798932   0.448246  0.00490252  0.000282716  0.51359    0.00493807  6.89176e-07   0.498317      0.00532263
10   121383201  C    T    0.000181172  0.102532  0.399575    0.0081343    0.0621894  0.399493    3.7699e-07    0.0833442     0.400034
...
```

The program doesn't use any nonstandard python libraries.  

The input summary stats should have chrom, pos, ref, alt as the first four columns and the other columns can be set in the config file. The summary stats can be space or tab separated. Currently error checking is minimal

