import json
import argparse
from typing import NamedTuple,List, Dict

class TimingData(NamedTuple):
    job: str
    status: str
    attempt: int
    start: str
    end: str
    n_phenos: int
    chrom: str
    cluster: str
    n_cpus: int

class Task(NamedTuple):
    status: str
    preemptible: bool
    attempt: int
    start: str
    end: str
    n_cpus: int
    mem: int
    data: Dict[str,str]

class WorkFlow(NamedTuple):
    id: str
    tasks: List[Task]


def extract_data(jsondata, tasklist:List[str], input_fields:List[str])->WorkFlow:
    """
    Extract tasks from a workflow.
    Args:
        jsondata (Any): the json data we parse
        tasklist (List[str]): List of tasks to extract
        input_fields (List[str]): List of input fields to extract from the tasks (e.g. chromosome, input file)
    """
    taskname = jsondata["id"]
    out_tasks = []
    for jobname,joblist in jsondata["calls"].items():
        if jobname in tasklist:
            for task in joblist:
                #filter out callcached results
                if "backendStatus" not in task.keys():
                    continue
                execstatus = task["executionStatus"]
                backendstatus = task["backendStatus"]
                if execstatus == "Done":
                    #load userAction
                    status = "Success"
                    actions = [a for a in task["executionEvents"] if "description" in a.keys()]
                    action = [a for a in actions if a["description"]=="UserAction"][0]
                elif execstatus == "RetryableFailure" and backendstatus == "Preempted":
                    status = "Preempted"
                    action = [a for a in task["executionEvents"] if a["description"]=="RunningJob"][0]
                elif execstatus == "Failed":
                    status = "Failed"
                    action = [a for a in task["executionEvents"] if a["description"]=="RunningJob"][0]
                else:
                    #couldn't identify the task type, might be callcaching, 
                    continue
                #get runtimeattributes
                start = action["startTime"]
                end = action["endTime"]
                preemptible = task["preemptible"]
                attempt = int(task["attempt"])
                #note: currently google adds at least 2 cpus for machines, change this if that becomes
                n_cpus = max(int(task["runtimeAttributes"]["cpu"]),2)
                mem = int(task["runtimeAttributes"]["memory"].split(" ")[0])
                #get input fields
                inputs = {}
                for field in input_fields:
                    try:
                        inputs[field] = task["inputs"][field]
                    except:
                        raise Exception(f"Field {field} not in inputs of task {jobname}")
                out_tasks.append(
                    Task(
                        status,
                        preemptible,
                        attempt,
                        start,
                        end,
                        n_cpus,
                        mem,
                        inputs
                    )
                )
    return WorkFlow(taskname,out_tasks)

def main(files,tasklist, input_fields, outfilename):
    # read files
    workflows = []
    for fname in files:
        with open(fname,"r") as f:
            d=json.load(f)
            workflows.append( extract_data(d,tasklist, input_fields) )
    #write output
    headercols = ["id","status","preemtible","attempt","start","end","n_cpus","mem"]+input_fields
    header = "\t".join( map(str,headercols) )+"\n"
    
    with open(outfilename,"w") as outfile:
        outfile.write(header)
        for workflow in workflows:
            for task in workflow.tasks:
                columns = [
                    workflow.id,
                    task.status,
                    task.preemptible,
                    task.attempt,
                    task.start,
                    task.end,
                    task.n_cpus,
                    task.mem,
                ]
                columns = columns + [task.data[i] for i in input_fields]
                line = "\t".join(map(str,columns))+"\n"
                outfile.write(line)

if __name__ == "__main__":
    parser=argparse.ArgumentParser("Extract timing data (including whether the task was pre-empted, and which attempt it was) from metadata file(s)")
    parser.add_argument("--out",required=True)
    parser.add_argument("files",nargs="+")
    parser.add_argument("--tasknames",nargs="+",help="list of task names to extract from workflow metadata")
    parser.add_argument("--inputfields",nargs="*",help="Extract input fields from ")
    args=parser.parse_args()
    main(args.files,args.tasknames, args.inputfields, args.out)