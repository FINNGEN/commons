from typing import List,NamedTuple,Tuple,Any,Union
import json
import requests
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import time


class Output(NamedTuple):
    cpra:List[str]
    rsid:List[str]
    prediction:List[str]
    inverted:bool

class Error(NamedTuple):
    cpra:List[str]
    reason:str

OutputResult = Union[Output,Error]

CHR_MAP={
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "23": "NC_000023.11",
    "24": "NC_000024.10",
    "25": "NC_012920.1"
}

def chr_to_ref_chr(chr:str)->str:
    return CHR_MAP[chr]

def get_request(cpra)->OutputResult:
    chr_seq = chr_to_ref_chr(cpra[0])
    hgvs = chr_seq+":g."+cpra[1]+cpra[2]+">"+cpra[3]
    hgvs_2 = chr_seq+":g."+cpra[1]+cpra[3]+">"+cpra[2]
    url = f"https://rest.ensembl.org/variant_recoder/human/{hgvs}?fields=id,hgvsp&content-type=application/json"
    while True:
        try:
            html = requests.get(url)
        except:
            time.sleep(1)
            continue
        else:
            break
    inverted=False

    if html.status_code!= 200:
        a=1
        url = f"https://rest.ensembl.org/variant_recoder/human/{hgvs_2}?fields=id,hgvsp&content-type=application/json"
        while True:
            try:
                html = requests.get(url)
            except:
                time.sleep(1)
                continue
            else:
                break
        inverted=True
        if html.status_code != 200:
            return Error(cpra,"returned invalid http code")

    #parse html
    jsond = json.loads(html.content)
    try:
        jsondata = jsond[0][cpra[3]]
    except:
        jsondata = jsond[0][cpra[2]]
    preds = jsondata["hgvsp"] if "hgvsp" in jsondata else []
    predictions = [a for a in preds if a.startswith("NP_")]
    #get rsid
    ids = jsondata["id"] if "id" in jsondata else []
    rsids = [a for a in ids if a.startswith("rs")]
    return Output(cpra,rsids,predictions,inverted)

def main(filename: str, column: str, out: str):
    with open(filename) as f:
        header = f.readline().strip("\n").split("\t")
        header_ord = {a:i for i,a in enumerate(header)}
        data_col = header_ord[column]
        cpras = []
        for dataline in f:
            snp = dataline.split("\t")[data_col]
            if snp == "":
                continue
            cpra = snp.split("_")
            cpra[0] = cpra[0].replace("chr","").replace("X","23").replace("Y","24")
            cpras.append(cpra)
        
        data:List[Output] = []
        processes = []
        with ThreadPoolExecutor(max_workers=10) as executor:
            for cpra in cpras:
                processes.append(executor.submit(get_request, cpra))
        

        for task in as_completed(processes):
            res = task.result()
            if isinstance(res,Output):
                data.append(res)
            else:
                print(res)
        """
        for cpra in cpras:
            res = get_request(cpra)
            if isinstance(res,Output):
                data.append(res)
            else:
                print(res)
        """
    #output to --out
    with open(out,"w") as outf:
        outf.write("variant\trsids\tpredictions\tinverted_alleles_in_db\n")
        for item in data:
            cpra = "chr"+"_".join(item.cpra)
            rsids = ";".join(item.rsid)
            preds = ";".join(item.prediction)
            inverted = "true" if item.inverted else "false"
            outf.write("\t".join([cpra,rsids,preds,inverted])+"\n")


        

if __name__=="__main__":
    ap = argparse.ArgumentParser("Fetch HGVS protein annotation for variants")
    ap.add_argument("FILE",help="tsv file")
    ap.add_argument("column",help="variant column")
    ap.add_argument("--out",required=True,help="output tsv with variant, rsid, hgvs protein annotation and whether the variant needed ref and alt inverted to get data")

    args = ap.parse_args()

    main(args.FILE,args.column,args.out)