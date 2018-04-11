

curl -X POST "http://0.0.0.0/api/workflows/v1/batch" -H  "accept: application/json" -H  "Content-Type: multipart/form-data" -F "workflowSource=@scripts/scrape_annot.wdl" -F "workflowInputs=@scripts/scrape_annot.wdl.json" --socks5 localhost:5000

curl -X GET "http://0.0.0.0/api/workflows/v1/07dd117a-07ea-4516-8bb7-f9b31315f9f1/metadata?expandSubWorkflows=false" -H  "accept: application/json" --socks5 localhost:5000 | python -m json.tool
