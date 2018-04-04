

curl -X POST "http://0.0.0.0/api/workflows/v1/batch" -H  "accept: application/json" -H  "Content-Type: multipart/form-data" -F "workflowSource=@scripts/scrape_annot.wdl" -F "workflowInputs=@scripts/scrape_annot.run.wdl" --socks5 localhost:5000


curl -X GET "http://0.0.0.0/api/workflows/v1/94d0116b-7974-44e1-8724-1760d105c390/metadata?expandSubWorkflows=false" -H  "accept: application/json" --socks5 localhost:5000 | python -m json.tool
