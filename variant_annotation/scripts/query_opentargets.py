#!/usr/bin/env python3
import requests
import json
import argparse
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import sys
import time


if sys.version_info[0] != 3 or sys.version_info[1]<7 or (sys.version_info[1]==7 and sys.version_info[2]<2):
    raise Exception(" Python Version >= 3.7.2 required")


query_string = """
            query PheWASQuery($variantId: String!) {
                pheWAS( variantId:$variantId ) {
                    totalGWASStudies
                    associations {
                      study {
                        traitReported
                        traitCategory
                        pubJournal
                        pubAuthor
                        pubTitle
                      }
                      beta
                      pval
                    }
                }
                indexVariantsAndStudiesForTagVariant(variantId:$variantId) {
                    associations {
                      indexVariant {
                        id
                        rsId
                      }
                      study {
                        traitReported
                        traitCategory
                        pmid
                        pubAuthor
                      }
                      pval
                      nTotal
                      overallR2
                      posteriorProbability
                      beta
                      pval
                    }
                }
                tagVariantsAndStudiesForIndexVariant(variantId:$variantId) {
                    associations {
                        tagVariant {
                            id
                            rsId
                        }
                        study {
                          traitReported
                          traitCategory
                          pmid
                          pubAuthor
                        }
                        beta
                        pval
                        overallR2
                    }
                }
            }
        """
base_url = "https://api.genetics.opentargets.org/graphql"

def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=tuple( x for x in requests.status_codes._codes if x != 400),
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def get_ot_anno(chrom:str, pos:str, ref:str, alt:str, retries=10, requests_session=None):

    var = f'{chrom.replace("chr","")}_{pos}_{ref}_{alt}'
    variables={"variantId":var}
    query= {"query": query_string, "variables": variables}

    print("####",file=sys.stderr)
    print(f'Querying {var}',file=sys.stderr)

    if not requests_session:
        requests_session = requests

    r = requests_session.post(base_url,json=query)
    if r.status_code!=200:
        raise Exception("Error requesting OpenTargets. Error " + r.text)

    response = json.loads(r.text)
    print(response["data"],file=sys.stderr)
    print("####",file=sys.stderr)
    return response["data"]

def run():
    parser = argparse.ArgumentParser(description='Add opentargets hit annotations and write to standard out.')
    parser.add_argument('result_file', action='store', type=str, help='The file to annotate')
    parser.add_argument('--chrom', default="chrom", type=str, help='')
    parser.add_argument('--pos', default="pos", type=str, help='')
    parser.add_argument('--ref', default="ref", type=str, help='')
    parser.add_argument('--alt', default="alt", type=str, help='')
    parser.add_argument('--sep', default="\t", type=str, help='')
    parser.add_argument('--p_threshold', default=5e-8, type=float, help='P value threshold to include results from OpenTargets')
    parser.add_argument('--r_threshold', default=0.6, type=float, help='r2 threshold to include tagging variants from opentargets')
    parser.add_argument('--ignore_studies', default="", type=str, help='comma separated list of studies to ignore. Eg. FINNGEN_R5')
    parser.add_argument('--ignore_categories', default="", type=str, help='comma separated list of categories to ignore. Eg. FINNGEN_R5')
    parser.add_argument('--parse_credsets_col', type=str, help='IF you have file with autoreporting style credsets (semicolon separated), providing the column name will trigger search with these variants')
    parser.add_argument('--max_cs_vars', default=3,type=int, help='Maximum number of TOP CS variants to use')
    parser.add_argument('--http_retries', default=10,type=int, help='')

    ##Uncategorised,measurement,UKB Neale v2

    args=parser.parse_args()

    ignorelist = [ ig.strip() for ig in args.ignore_studies.split(",") ]
    rsesh = requests_retry_session(retries=args.http_retries)

    with open(args.result_file) as infile:
        h_idx = { h:i for i,h in enumerate(infile.readline().rstrip("\n").split(args.sep)) }

        head = list(h_idx.keys())
        head.extend(["OT_phewas","OT_tagging"])
        print(args.sep.join(head))

        for l in infile:
            dat = l.rstrip("\n").split(args.sep)

            csvars = []
            if args.parse_credsets_col:
                varstosearch = [ csdat.split("|") for csdat in dat[h_idx[args.parse_credsets_col]].split(";") ]
                varstosearch = sorted(varstosearch, key=lambda x: x[1], reverse=True)
                varstosearch = [ v[0].split("_")  for v in varstosearch[0:min(args.max_cs_vars,len(varstosearch))] ]
            else:
                varstosearch = [(dat[h_idx[args.chrom]], dat[h_idx[args.pos]],
                    dat[h_idx[args.ref]], dat[h_idx[args.alt]] )]

            phewas = set()
            tagging = set()
            manual_retries = 3
            success = False
            for v in varstosearch:
                while not success and manual_retries>0:
                    try:
                        odat = get_ot_anno(v[0],v[1],v[2],v[3], rsesh)
                        success = True
                    except BaseException as err:
                        print(f'Error getting or retrying {err}', file=sys.stderr)
                        time.sleep(10)
                        manual_retries-=1
                        success = False

                if not success:
                    raise Exception("Connection failed multiple times.")

                for s in odat["pheWAS"]["associations"]:
                    if s["study"] is not None and not s["study"]["pubAuthor"] in ignorelist:
                        beta = "{:.2f}".format(s["beta"])
                        pval = "{:.2e}".format(s["pval"])
                        trait = s["study"]["traitReported"]
                        trait_cat = s["study"]["traitCategory"]
                        study = s["study"]["pubAuthor"]
                        if(s["pval"] is not None and s["pval"]<args.p_threshold):
                            print("adding" + ",".join([trait,trait_cat,beta,pval,study]), file=sys.stderr)
                            phewas.add(",".join([trait,trait_cat,beta,pval,study]))

                for t in odat["indexVariantsAndStudiesForTagVariant"]["associations"]:
                    if not t["study"]["pubAuthor"] in ignorelist:
                        beta = "{:.2f}".format(t["beta"]) if t["beta"] is not None else "NONE"
                        pval = "{:.2e}".format(t["pval"])
                        trait = t["study"]["traitReported"]
                        trait_cat = t["study"]["traitCategory"]
                        study = t["study"]["pubAuthor"]
                        r2 = "{:.2f}".format(t["overallR2"]) if t["overallR2"] is not None else "NONE"
                        var=t["indexVariant"]["id"]
                        study = t["study"]["pubAuthor"]
                        if(t["pval"] is not None and t["pval"]<args.p_threshold and (r2=="NONE" or float(r2)>args.r_threshold)):
                            tagging.add(",".join([trait,trait_cat,beta,pval,var,r2,study]))

                for t in odat["tagVariantsAndStudiesForIndexVariant"]["associations"]:
                    if not t["study"]["pubAuthor"] in ignorelist:
                        beta = "{:.2f}".format(t["beta"]) if t["beta"] is not None else "NONE"
                        pval = "{:.2e}".format(t["pval"])
                        trait = t["study"]["traitReported"]
                        trait_cat = t["study"]["traitCategory"]
                        study = t["study"]["pubAuthor"]
                        r2 = "{:.2f}".format(t["overallR2"]) if t["overallR2"] is not None else "NONE"
                        var=t["tagVariant"]["id"]
                        study = t["study"]["pubAuthor"]
                        if(t["pval"] is not None and t["pval"]<args.p_threshold and (r2=="NONE" or float(r2)>args.r_threshold)):
                            tagging.add(",".join([trait,trait_cat,beta,pval,var,r2,study]))

            print(f'Adding {phewas}', file=sys.stderr)
            dat.append(";".join(phewas))
            dat.append(";".join(tagging))
            print(args.sep.join(dat))

if __name__ == '__main__':
    run()
