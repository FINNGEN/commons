#!/bin/bash

while read line; do
    curl -X GET "http://0.0.0.0/api/workflows/v1/${line}/metadata?expandSubWorkflows=false" -H "accept: application/json" > json_main/${line}.json
done < workflows.txt

