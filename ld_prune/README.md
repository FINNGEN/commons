# Script to LD cluster arbitrary text file rows containing chr, pos, ref, alt

## Installation

Python >=3.6 required. f-strings won't work on python < 3.6.
Install dependencies with pip install -r requirements.txt

## Running
Script takes tab-delimited, chr pos ordered `file` (order checked but not sorted, must pre-sort) with headers and considers adjacent lines for ld pruning if they are each `ld_w` basepairs apart. The `ld_w` is taken from the closest point of formed cluster (i.e. lines can get chained into single cluster this spanning larger area than `ld_w`).
All lines that have r2> `ld` will be clustered together, most important line output as is and selected columns (`prune_column_list`) will be output in additional column in the end.

LD source is FinnGen LD api  (i.e. FinnGen v3 imputation panel).
Hitting LD api is the slow part so consider running in batches if using a lot of data (or add threading to ld api access pull request:) )

`usage: ld_prune_lines.py [-h] [-pcol PCOL] [-chromcol CHROMCOL] [-poscol POSCOL] [-refcol REFCOL] [-altcol ALTCOL] [-ld LD] [-ld_w LD_W] [-prune_column_list PRUNE_COLUMN_LIST] file outfile`

Parameters:

  - file: input
  - outfile: output

  - -pcol: name for column containing priority for lines (must be convertable to float, smaller is more important e.g. pvalue)
  - -chromcol: chromosome column name
  - -refcol: ref column
  - -altcol: alt column
  - -ld: r2 above which rows are clustered
  - -ld_w: width of adjacent hits to consider for LD pruning. cluster width can be wider as this is for adjacent hits.
  - -prune_column_list: what columns to output from those lines that get merged to priority line
  - -ld_source: sisu3 or sisu4 currently supported. Default sisu4
  - fixed_ld_search_width: give this and the search region call to LD api is not made as small as possible but always uses this size.
  This is interim fix as currently LD server does not respect search boundaries exactly
