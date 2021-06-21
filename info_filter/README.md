# Filter variants based on a list

This WDL filters a list of summary statistics using a list of allowed variants. The variants should be written with c:p:r:a, i.e. one variant per row, with their chromosome, position, reference and alternate alleles, and with the fields separated with semicolon (:). For example:
```
1:1:A:T
1:5:C:G
2:5:CT:C
```
The variants do not have to be sorted. The summary statistic files need to be sorted by chromosome and position to allow tabix-indexing.

## Usage
Fill in `sumstats.loc` and `variant_file` to `filter.json` with the first being a list of summary staitistc locations and the second being the variant list.

The outputs can be accessed using the `outputs` endpoint in cromwell. The outputs are tabix-indexed.