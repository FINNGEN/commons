"""Chunk a file into multiple files
"""

import argparse as ap
import itertools

def grouper_it(n, iterable):
    it = iter(iterable)
    while True:
        chunk_it = itertools.islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,), chunk_it)

def chunk_file(fname, n_lines):
    """Chunk a file into multiple lines
    Save those files as filename_start_end
    """
    with open(fname,"r") as f:
        for idx,chunk in enumerate(grouper_it(n_lines,f)):
            ch = list(chunk)
            start = n_lines*idx +1
            end = min(n_lines*(idx+1),(n_lines*idx + len(ch)))
            print(start,end)
            with open(f"{fname}_{start}_{end}","w") as writer:
                writer.writelines(ch)


if __name__=="__main__":
    parser = ap.ArgumentParser("Chunk a file into chunks of size n")
    parser.add_argument("file")
    parser.add_argument("--n-lines",type=int,help="lines per chunk")

    args = parser.parse_args()

    chunk_file(args.file, args.n_lines)
