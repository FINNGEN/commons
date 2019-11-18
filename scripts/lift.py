#! /usr/bin/env python3

import argparse,gzip,sys,subprocess,os,shlex
from functools import partial
from tempfile import NamedTemporaryFile


def get_dat_var(line, index):
    d = line[index].split(":")
    if len(d)<4:
        print("WARNING: Not properly formatted variant id in line: " + line, file=sys.stderr )
        return None
    return d


def return_open_func(f):
    '''
    Detects file extension and return proper open_func
    '''
   
    file_path = os.path.dirname(f)
    basename = os.path.basename(f)
    file_root, file_extension = os.path.splitext(basename)
    
    if 'bgz' in file_extension:
        #print('gzip.open with rb mode')
        open_func = partial(gzip.open, mode = 'rb')
    
    elif 'gz' in file_extension:
        #print('gzip.open with rt mode')
        open_func = partial(gzip.open, mode = 'rt')

    else:
        #print('regular open')
        open_func = open      
    return open_func

def lift(args):
    '''
    Finds the proper header column names and calls the liftover.sh script
    '''
    open_func = return_open_func(args.file)
    with open_func(args.file) as res:
        #skip first line anyways
        header = res.readline().rstrip("\n").split()
        if args.var:
            print(args.var)
            if not args.numerical:
                args.var = header.index(args.var)
                
            print(args.var)
            joinsortargs =f"--var {args.var+1}"
            get_dat_func = partial(get_dat_var,index=args.var) 

        elif args.info:
            if not args.numerical:
                args.info = [header.index(elem) for elem in args.info]

            chrom,pos,ref,alt = args.info
            print(args.info)
            joinsortargs = f"--chr {chrom+1} --pos {pos+1} --ref {ref+1} --alt {alt+1}"
            get_dat_func = lambda line:  (line[chrom], line[pos], line[ref], line[alt])

        tmp_bed = NamedTemporaryFile(delete=True)
        with open(tmp_bed.name, 'w') as out:
            for line in res:
                vardat = get_dat_func(line.strip().split())
                string = "{}\t{}\t{}\t{}".format("chr"+vardat[0], str(int(vardat[1])-1), str(int(vardat[1]) + max(len(vardat[2]),len(vardat[3])) -1), ":".join([vardat[0],vardat[1],vardat[2],vardat[3]])) + "\n"
                out.write(string)
                
    cmd = f"liftOver {tmp_bed.name} {args.chainfile} variants_lifted errors"
    subprocess.run(shlex.split(cmd))
    
    joinsort = f"{os.path.join(args.scripts_path,'joinsort.sh')}"
    subprocess.run(shlex.split(f"chmod +x {joinsort}"))
    
    joincmd = f"{joinsort} {args.file} variants_lifted {joinsortargs}"
    subprocess.run(shlex.split(joincmd))

    if args.out:
        mv_cmd = f"mv -f ./{os.path.basename(args.file)}.lifted.gz ./{os.path.basename(args.file)}.lifted.gz.tbi variants_lifted errors {args.out}"
        subprocess.run(shlex.split(mv_cmd))
    
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Add 38 positions to summary file')
    parser.add_argument("file", help=" Whitespace separated file with either single column giving variant ID in chr:pos:ref:alt or those columns separately")
    parser.add_argument("--chainfile", help=" Chain file for liftover",required = True)
    parser.add_argument("--out", help=" Folder where to save the results",required = False)
        
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('--info',nargs =4, metavar = ('chr','pos','ref','alt'), help = 'Name of columns')
    group.add_argument("--var",help ="Column name if : separated")
    
    args = parser.parse_args()

    # checks if var/info are numerical or strings
    args.numerical = False
    if args.var and args.var.isdigit():
        args.var = int(args.var)
        args.numerical = True
    if args.info and  all(elem.isdigit() for elem in args.info):
        args.info = list(map(int,args.info))
        args.numerical = True

    if not args.out:
        args.out = '/'.join(args.file.split('/')[:-1])
    args.scripts_path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
        
    
    print(args)
    
    lift(args)

