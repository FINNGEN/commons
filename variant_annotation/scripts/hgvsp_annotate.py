from typing import List,NamedTuple,Tuple,Any,Union
import json
import requests
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import hgvs

class Output(NamedTuple):
    cpra:List[str]
    rsid:List[str]
    prediction:List[str]

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


def vcf_to_vep_input(cpra:List[str])->str:
    """construct vep input, like 4:87845066-87845070:1/G from chr4_87845066_GGAAA_G
    """
    chrom = cpra[0].replace("chr","").replace("23","X")
    start = int(cpra[1])
    #if it's a deletion, then add to end
    if len(cpra[2]) > len(cpra[3]):
        end = start + max(0,len(cpra[2])-1)
    elif len(cpra[3]) > len(cpra[2]):#ins
        #TODO figure this out
        end = start
    else:
        end=start
    return f"{chrom}:{start}-{end}:1/{cpra[3]}"

def parse_vep_output(data:Any)->List[str]:
    """Parse VEP output and return vhgs coding sequences
    """
    out = []
    if "transcript_consequences" in data[0]:
        transcr_cons = data[0]["transcript_consequences"]
        for cons in transcr_cons:
            #take hgvsc out
            if "hgvsc" in cons:
                out.append(cons["hgvsc"])
    return out


def vv_request(cpra:List[str])->OutputResult:
    pvcf = ":".join(cpra)
    url=f"https://rest.variantvalidator.org/VariantFormatter/variantformatter/GRCh38/{pvcf}/all/None/False?content-type=application/json"
    while True:
        try:
            html = requests.get(url)
        except:
            time.sleep(1)
            continue
        if html.status_code in (200,404):
            break
        elif html.status_code in range(500,550):
            time.sleep(1)
            continue
        else:
            return Error(cpra, f"html returned the following stuff with code {html.status_code}: "+str(html.content,encoding="utf-8"))
    data = json.loads(html.content)
    inner_data = data[pvcf][pvcf]
    predictions = []
    if inner_data["genomic_variant_error"]== None:
        try:
            for _key,value in inner_data["hgvs_t_and_p"].items():
                if "p_hgvs_tlc" in value:
                    predictions.append(value["p_hgvs_tlc"])
        except:
            print(inner_data)
            print(data)
            raise
    predictions = [a for a in predictions if a != None]
    return Output(cpra,[],predictions)

    



def parse_recoder_data(cpra:List[str], data:Any)->OutputResult:
    key = list(data[0].keys())[0]
    jsondata=data[0][key]
    preds = jsondata["hgvsp"] if "hgvsp" in jsondata else []
    predictions = [a for a in preds if a.startswith("NP_")]
    #get rsid
    ids = jsondata["id"] if "id" in jsondata else []
    rsids = [a for a in ids if a.startswith("rs")]
    return Output(cpra,rsids,predictions)

def vep_request(cpra:List[str])->OutputResult:
    pass
    vep_input = vcf_to_vep_input(cpra)
    url_vep = f"https://rest.ensembl.org/vep/human/region/{vep_input}?vcf_string=true&hgvs=true&content-type=application/json"
    while True:
        try:
            html = requests.get(url_vep)
        except:
            time.sleep(1)
            continue
        else:
            break
    inputs = None
    if html.status_code == 200:
        inputs = parse_vep_output(json.loads(html.content))
    if not inputs:
        print(":".join(cpra),vep_input, str(html.content,encoding="utf-8"))
        return Output(cpra,[],[])
    data_to_parse = []
    for entry in inputs:
        url2 = f"https://rest.ensembl.org/variant_recoder/human/{entry}?fields=id,hgvsp&content-type=application/json"
        while True:
            try:
                html = requests.get(url2)
            except:
                time.sleep(1)
                continue
            else:
                break
        if html.status_code == 200:
            data_to_parse.append(json.loads(html.content))
        else:
            print(cpra,vep_input,entry,inputs, str(html.content,encoding="utf-8"))
            return Output(cpra,[],[])
    parsed = [parse_recoder_data(cpra,a) for a in data_to_parse]
    errors = [a for a in parsed if type(a)==Error]
    for e in errors:
        print(e)
    values = [a for a in parsed if type(a) == Output]
    #fold somehow.
    out_rsids = []
    out_hgvsp = []
    _ = [out_rsids.extend(a.rsid) for a in values ]
    _ = [out_hgvsp.extend(a.prediction) for a in values ]
    rsids = list(set(out_rsids))
    predictions = list(set(out_hgvsp))
    return Output(cpra,rsids,predictions)
    



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
        with ThreadPoolExecutor(max_workers=4) as executor:
            for cpra in cpras:
                processes.append(executor.submit(vv_request, cpra))
        

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
        outf.write("variant\trsids\tpredictions\n")
        for item in data:
            try:
                cpra = "chr"+"_".join(item.cpra)
                rsids = ";".join(item.rsid)
                preds = ";".join(item.prediction)
                outf.write("\t".join([cpra,rsids,preds])+"\n")
            except:
                print(item)

        

if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("FILE")
    ap.add_argument("column")
    ap.add_argument("--out",required=True)

    args = ap.parse_args()

    main(args.FILE,args.column,args.out)
