# Divide a list into many N-sized chunks

This script divides a file with many lines into many files with N or less lines.  The line order or line contents (including line endings) is not changed.

For example, to divide a phenotype list into files with 500 lines:

```
python3 chunk_list.py FILE --n 500
ls -1
> FILE_1_500
> FILE_501_1000
> FILE_1001_1501
#etc etc
```

## Usage
```
usage: Chunk a file into chunks of size n [-h] [--n-lines N_LINES] file

positional arguments:
  file

optional arguments:
  -h, --help         show this help message and exit
  --n-lines N_LINES  lines per chunk
```
