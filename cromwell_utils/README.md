# Utilities for metadata

This folder contains utilities for working with cromwell-produced files (metadata, outputs). The following scripts are provided:
- `outputs_to_lists.py`, which extracts outputs to their own lists from an outputs.json file. Outputs of the same type are put in same files, and nested arrays are flattened. Those can be used in copying files, for example: `cat filelist|gsutil -m cp -I gs://output_bucket/`
- `extract_timing.py`, which extracts timing data from a list of workflows' metadata files and saves it as a TSV.
- `get_subworkflowids.py`, which extracts subworkflow IDs from a metadata json. This is useful for e.g. bulk downloading all subworkflow metadata

## Requirements

These scripts require python 3.6+

## Usage

Get a list of subworkflowids from a metadata file & download the associated metadata:

```bash
PROXY_PORT=5000
python3 get_subworkflowids.py metadata.json
#subworkflowids saved to subworkflowid named files, such as regenie.step2
#download metadata, 4 parallel processes
mkdir subworkflows
cat regenie.step2|xargs -I %  -P 4 bash -c 'curl --sS -X POST http://localhost:80/api/workflows/v1/{}/metadata -H "accept: application/json" -H "Content-Type: multipart/form-data" --socks5 localhost:PROXY_PORT -o subworkflows/{}.json'
```

extract task timing from the subworkflows:
```bash
python3 extract_timing.py subworkflows/*.json --tasknames "regenie_step2.step2" --inputfields "input1" "input2" "input3"  --out extracted_data_step2.tsv
```
this timing data can then be drilled down using e.g. R to see pre-emption rates, calculate approximate cost etc.

Download outputs from a workflow:
```bash
python3 outputs_to_lists.py outputs.json --out-folder output_lists
ls output_lists
>workflow.output1
>workflow.output2
>workflow.output3
```