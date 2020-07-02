# lift


Scripts to lift over a sumstat file to any other build through a chainfile.


## Usage


```
usage: lift.py [-h] --chainfile CHAINFILE [-o OUT] [--sep SEP]
               (--info chr pos ref alt | --var snpid snp_sep)
               file

Add liftover positions to summary file

positional arguments:
  file                  Whitespace separated file with either single column
                        giving variant ID in chr:pos:ref:alt or those columns
                        separately

optional arguments:
  -h, --help            show this help message and exit
  --chainfile CHAINFILE
                        Chain file for liftover
  -o OUT, --out OUT     Folder where to save the results
  --sep SEP             column separator in file to be lifted. Default tab
  --info chr pos ref alt
                        Name of columns
  --var snpid snp_sep   Column name of snpid and separator
```

Both `--info` and `--var` can be passed as column indexes (0 based) or column names.

