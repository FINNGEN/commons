## Useful bash commands go here


### Get wdl workflow from a previously submitted workflow's metadata

Suppose you have a Cromwell server SOCKS proxy set up on port 5000:

```
curl --preproxy localhost:5000 http://0.0.0.0/api/workflows/v1/WORKFLOW_ID/metadata | jq -r '.submittedFiles.workflow' > WORKFLOW.wdl
```

Replace WORKFLOW_ID and WORKFLOW.

To get inputs:

```
curl --preproxy localhost:5000 http://0.0.0.0/api/workflows/v1/WORKFLOW_ID/metadata | jq -r '.submittedFiles.inputs' | jq -r '' > WORKFLOW.json
```

This is useful if you forget your laptop at home.
