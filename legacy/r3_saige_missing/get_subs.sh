#!/bin/bash

for file in json_main/*.json; do
    while read line; do
	curl -X GET "http://0.0.0.0/api/workflows/v1/${line}/metadata?expandSubWorkflows=false" -H "accept: application/json" > json_sub/${line}.json
    done < <(python -m json.tool $file | grep subWorkflowId | cut -d':' -f2 | cut -d'"' -f2)
done
