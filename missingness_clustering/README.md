# Cluster phenotypes with low missingness together

Cluster phenotypes that have low missingness between them together. Required for e.g. regenie, where binary traits that are processed together should have low missingness, < 5%. In case the exact same missingness calculation has been performed earlier on that machine, it is loaded from cache.

## Methods

Missingness is calculated in one of two ways:
### XOR/AND (default)
Given two phenotypes, a and b, with some samples not assigned to cases or controls, encode a sample being a case/control as 1 and not being either as 0.
For example:
```
a = [1,1,1,1,0,0,1]
b = [1,1,1,1,0,0,0]
```
Then, missing samples are those that are in one or the other phenotype, but not in both. This can be efficiently calculated with XOR(a,b).
To get a missingness ratio, the amount of missing samples (SUM(XOR(a,b))) can be divided with either 
    a) the amount of samples common to both of the phenotypes, calculated efficiently with AND(a,b), or
    b) the amount of samples in either of the phenotypes, calculated with OR(a,b).  
These would then be
```
    and_missingness(a,b) = SUM(XOR(a,b))/SUM(AND(a,b)))
    or_missingness(a,b) = SUM(XOR(a,b))/SUM(OR(a,b)))
```
for the aforementioned phenotypes, these return a missingness ratio of 0.25 and 0.2, respectively.
### Clustering
After the script has calculated a pairwise missingness matrix between all of the available phenotypes, this is passed to scikit-learn's agglomerative clustering solver (see [here](https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering) and [here](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html#sklearn.cluster.AgglomerativeClustering)). the clustering is "complete", meaning that in every cluster, all of the members of that cluster have missingness < the specified missingness threshold.
## Installation
```
python3 -m pip install -r requirements.txt
```

## Usage
```
usage: cluster data into phenotype clusters with a specified limit on missingness
       [-h] [--max-missingness MAX_MISSINGNESS] [--max-size MAX_SIZE]
       --phenolistfile PHENOLISTFILE [--test] [--out OUT]
       [--missingness-method {and,or}] [--force]
       cov_file

positional arguments:
  cov_file              covariate file

optional arguments:
  -h, --help            show this help message and exit
  --max-missingness MAX_MISSINGNESS
                        maximum missingness between phenotypes in cluster
  --max-size MAX_SIZE   maximum cluster size. CLustering is done to get the
                        largest clusters possible, and these will be divided
                        if necessary.
  --phenolistfile PHENOLISTFILE
                        list of phenotypes to cluster
  --test                Test with only the first 10k samples
  --out OUT             output cluster file
  --missingness-method {and,or}
                        Missingness calculated with either XOR(A,B)/AND(A,B)
                        or XOR(A,B)/OR(A,B). and is stricter (larger
                        missingness values for same at least somewhat non-
                        overlapping phenotypes).
  --force               do not use cache even if missingness data is available

```

## Examples
Cluster phenotypes using XOR/AND missingness, max 4 phenos per cluster
```
python3 cluster_phenos.py COV_FILE.txt.gz --phenolistfile PHENOLIST --max-missingness 0.05 --max-size 4  --out OUTPUT 
```
Cluster phenotypes, but forgo cache
```
python3 cluster_phenos.py COV_FILE.txt.gz --phenolistfile PHENOLIST --max-missingness 0.05 --max-size 4  --out OUTPUT --force
```
Cluster phenotypes, use XOR/OR missingness calculation
```
python3 cluster_phenos.py COV_FILE.txt.gz --phenolistfile PHENOLIST --max-missingness 0.05 --max-size 4  --out OUTPUT --missingness-method or
```