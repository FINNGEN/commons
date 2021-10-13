import argparse
import json
from collections import defaultdict

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("files",nargs = "+")
    args=parser.parse_args()
    calls = defaultdict(list)
    for fname in args.files:
        with open(fname) as f:
            data=json.load(f)
            for key,value in data["calls"].items():
                for call in value:
                    calls[key].append(call['subWorkflowId'])
    for key, values in calls.items():
        with open(key,"w") as outfile:
            for v in values:
                outfile.write(f"{v}\n")