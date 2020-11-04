#!/usr/bin/env node

const fs = require('fs')

const phenoHash = fs.readFileSync('R3_phenos_remaining_formanual_final.txt', 'utf8').trim().split(/\r?\n|\r/).reduce((acc, cur) => { acc[cur] = true; return acc}, {})
const workflows = fs.readFileSync('workflows.txt', 'utf8').trim().split(/\r?\n|\r/)
const jsons = workflows.map(wf => JSON.parse(fs.readFileSync(`json_main/${wf}.json`, 'utf8')))

const phenoNullHash = {}
jsons.forEach(json => {
    json.calls['saige.null'].forEach(nul => {
	const pheno = nul['inputs']['pheno']
	if (phenoHash[pheno] && !phenoNullHash[pheno] && nul['executionStatus'] == 'Done' && nul['returnCode'] == '0') {
	    phenoNullHash[pheno] = nul['outputs']['modelfile']
	}
    })
})

Object.keys(phenoNullHash).forEach(pheno => {
    const missing = {}
    for (let i = 1; i <= 23; i++) {
	missing[`gs://r3_data/bgen/fromdatateam/${i}R3.bgen`] = true
    }
    jsons.forEach(json => {
	json.calls['saige.test_combine'].forEach(test => {
	    if (test.inputs.pheno == pheno) {
		const sub = JSON.parse(fs.readFileSync(`json_sub/${test.subWorkflowId}.json`, 'utf8'))
		sub.calls['test_combine.test'].forEach(task => {
		    if (task['executionStatus'] == 'Done' && task['returnCode'] == '0') {
			missing[task.inputs['bgenfile']] = false
		    }
		})
	    }
	})
    })
	miss = Object.keys(missing).filter(m => missing[m])
	if (miss.length > 0) {
	    fs.writeFileSync(`path_null/${pheno}.null.txt`, phenoNullHash[pheno] + '\n')
	    fs.writeFileSync(`path_bgen/${pheno}.bgen.txt`, miss.join('\n') + '\n')
	}
})
