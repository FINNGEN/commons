"""
Cluster phenotypes based on missingness
Missingness calculated as proportion of samples not common to both phenotypes/max(proportion of samples common to both phenotypes,1)

Calculating the missingness is cached if possible.
"""
import argparse
import hashlib
import os

from typing import List, Any, Optional
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn.cluster as cluster

HASH_BUFFER = 256*1024 #256kb buffer for hashing


def divide_chunks(l: List[Any], n: int) -> List[List[Any]]:
    """
    Divide list of values into m lists of n values
    """
    # looping till length l
    for i in range(0, len(l), n): 
        yield l[i:i + n]

def load_phenos(phenofile: str) -> List[str]:
    return list(np.loadtxt(phenofile, dtype=str))

def hash(pheno_file: str, cov_file: str, test: bool, method: str) -> str:
    """
    Return md5sum of file
    """
    hasher = hashlib.md5()
    hasher.update(bytes(test))
    hasher.update(bytes(method,encoding="utf-8"))
    with open(cov_file, 'rb') as f:
        while True:
            data = f.read(HASH_BUFFER)
            if not data:
                break
            hasher.update(data)
    with open(pheno_file, 'rb') as f:
        while True:
            data = f.read(HASH_BUFFER)
            if not data:
                break
            hasher.update(data)
    return hasher.hexdigest()

def load_covariate_matrix(fpath: str, phenotypes: List[str], nrows: Optional[int]) -> np.ndarray:
    """
    Load a covariate matrix from data
    Args:
        fpath: filepath to covariate file
        phenotypes: list of phenotypes to load
    Returns:
        (np.ndarray): MXN array, where M is the number of phenotypes and N is the number of samples
    """
    data = 1-pd.read_csv(fpath,sep="\t",usecols=phenotypes,nrows=nrows)[phenotypes].isna().astype(np.int8)
    return data.values.T

def calculate_missingness(data: np.ndarray, method: str) -> np.ndarray:
    """
    Calculate missingness matrix line-by-line since with whole releases,
    calculating it as a single matrix operation would require multiple terabytes of memory. 
    """
    if method == "and":
        divisor = np.logical_and
    elif method == "or":
        divisor = np.logical_or
    else: #this should be a programming error as argparse only allows correct values
        raise Exception(f"Error: invalid missingness method: {method}")

    out_matrix = np.zeros((data.shape[0],data.shape[0]))
    #iterate on columns
    for idx in range(0,data.shape[0]):
        print("\rProgress: {:>6.3g}%".format(100*idx/data.shape[0]),end="")

        div_vec = np.sum(divisor(data[idx,:],data),axis=1)
        #fix divide by 0 by filling zeros with 1
        div_vec = np.where(div_vec <=0,1,div_vec)
        xor_vec = np.sum(np.logical_xor(data[idx,:],data),axis=1)
        missingness_vec=xor_vec/div_vec
        out_matrix[idx,:]=missingness_vec
    print()
    return out_matrix

def cluster_phenotypes(missingness_matrix: np.ndarray, phenos: List[str], missingness_limit:float)->List[List[str]]:
    """
    Cluster phenotypes using agglomerative clustering with complete linkage
    Agglomerative clustering: merge clusters (initially clusters containing a single member) until a threshold (missingness limit between members of cluster) is reached
    complete linkage: all members of cluster have less missingness between them than the missingness limit
    """
    clusters = cluster.AgglomerativeClustering(n_clusters=None,affinity="precomputed",linkage="complete",distance_threshold=missingness_limit).fit(missingness_matrix)
    in_clusters = defaultdict(list)
    pheno_colmap = {idx:pheno for idx,pheno in enumerate(phenos) }
    for idx,group in enumerate(list(clusters.labels_)):
        in_clusters[group].append(pheno_colmap[idx])
    
    return list(in_clusters.values())

def size_clusters(protoclusters: List[List[str]],max_size: int)->List[List[str]]:
    """
    Take in a list of clusters, and divide too large clusters to new clusters that are at most size max_size
    """
    out = []
    for pc in protoclusters:
        if len(pc)> max_size:
            out.extend(divide_chunks(pc,max_size))
        else:
            out.append(pc)
    return out

def save_clusters(clusters:List[List[str]],fname:str):
    with open(fname,"wt") as f:
        for cluster in clusters:
            data = "\t".join(cluster)
            f.write(f"{data}\n")

class ResultCache():
    """
    Implement simple cache for caching calculated missingness matrices
    """

    def __init__(self, parent_dir:Path):
        self.path = parent_dir / ".cache"
        if self.path.exists():
            return
        else:
            os.makedirs(self.path)
            return

    def in_cache(self, hashstr:str)->bool:
        return os.path.exists(self.path / f"{hashstr}.npz")

    def load_covariate_matrix(self, hashstr:str)->np.ndarray:
        try:
            load_path = self.path / f"{hashstr}.npz"
            with np.load(load_path) as f:
                data=f["data"]
            return data
        except:
            print("Problem loading file from cache")
            raise

    def store_missingness(self, hashstr:str, in_data:np.ndarray):
        save_path = self.path / f"{hashstr}.npz"
        try:
            np.savez_compressed(save_path,data=in_data)
        except:
            print("error occurred in saving missingness data")
            raise

if __name__ == "__main__":
    parser=argparse.ArgumentParser("cluster data into phenotype clusters with a specified limit on missingness")
    parser.add_argument("cov_file",help="covariate file")
    parser.add_argument("--max-missingness",type=float,default=0.05,help="maximum missingness between phenotypes in cluster")
    parser.add_argument("--max-size",type=int,help="maximum cluster size. CLustering is done to get the largest clusters possible, and these will be divided if necessary.")
    parser.add_argument("--phenolistfile",type=str,required=True,help="list of phenotypes to cluster")
    parser.add_argument("--test",action="store_true",help="Test with only the first 10k samples")
    parser.add_argument("--out",help="output cluster file")
    parser.add_argument("--missingness-method",choices=["and","or"],default="and", help="Missingness calculated with either XOR(A,B)/AND(A,B) or XOR(A,B)/OR(A,B). and is stricter (larger missingness values for same at least somewhat non-overlapping phenotypes).")
    parser.add_argument("--force",action="store_true",help="do not use cache even if missingness data is available")

    args=parser.parse_args()

    pheno_file = args.phenolistfile
    cov_file = args.cov_file
    nrows = None
    if args.test:
        nrows = 10000
    missingness_limit = args.max_missingness
    output_fname = args.out
    cluster_size = args.max_size

    #get source file location for caching purposes
    source_path = Path(__file__).resolve()
    source_dir = source_path.parent

    cache = ResultCache(source_dir)

    
    phenos = load_phenos(pheno_file)

    #load data: load from cache or calculate missingess if not available
    print("hash inputs...")
    cachehash = hash(pheno_file,cov_file,args.test,args.missingness_method)
    print(f"DEBUG: hash is {cachehash}")

    if cache.in_cache(cachehash) and (not args.force) :
        print("data in cache, loading...")
        missingness_matrix = cache.load_covariate_matrix(cachehash)
    else:
        print("calculating missingness...")
        cov = load_covariate_matrix(cov_file,phenos,nrows)
        missingness_matrix = calculate_missingness(cov,args.missingness_method)

    #store missingness matrix in cache
    if (not cache.in_cache(cachehash)) or args.force:
        print("saving missingness matrix in cache...")
        cache.store_missingness(cachehash,missingness_matrix)
    #calculate clusters
    print("clustering phenotypes...")
    protoclusters = cluster_phenotypes(missingness_matrix,phenos,missingness_limit)
    #resize clusters (divide too large clusters)
    clusters = size_clusters(protoclusters, cluster_size)
    #output
    save_clusters(clusters,output_fname)
    