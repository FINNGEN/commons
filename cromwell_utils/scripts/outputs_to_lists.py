import json
import argparse
import os
from typing import List

def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list):
            rt.extend(flatten(i))
        else:
            rt.append(i)
    return rt

def main(files:List[str],out_folder:str):
    data = {}
    for f_path in files:
        with open(f_path,"r") as f:
            json_data=json.load(f)
            outs = json_data["outputs"]
            data[f_path]=outs
    #data layout: data contains dictionaries with values being the outputs.
    #If we want to join them, we need to 1) get keys, 2) 
    out_types = []
    for v in data.values():
        out_types.extend(list(v.keys()))
    all_keys = sorted(set(out_types))
    outputs = {}
    for key in all_keys:
        outputs[key]=[]
    for key in all_keys:
        for batchdata in data.values():
            if key in batchdata:
                tmp=flatten(batchdata[key])
                outputs[key].extend(tmp)
    #create folder if necessary
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    #write outputs
    for key,out in outputs.items():
        #out might still need to be unrolled
        tmp = flatten(out)
        fname = os.path.join(out_folder,key)
        with open(fname, "w") as f:
            for line in tmp:
                f.write(f"{line}\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser("create output lists from outputs jsons")
    ap.add_argument("files",metavar="FILES",nargs="+")
    ap.add_argument("--out-folder",type=str, required=True)
    args = ap.parse_args()

    main(args.files, args.out_folder)