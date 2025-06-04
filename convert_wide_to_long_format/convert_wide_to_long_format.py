import os
import argparse
import gzip
import re

COLUMN= 'COLUMN';
VALUE= 'VALUE';

def map_row_with_header_index(row, header_index, output_file):
    for header_row in header_index:
        new_row_dict = {
            key: row[idx] if isinstance(idx, int) and 0 <= idx < len(row) else idx
            for key, idx in header_row.items()
        }
        output_file.write("\t".join(str(v) for v in new_row_dict.values()) + "\n")

def run(input_file, output_file, custom_columns):
    # write the array of objects to a csv file
    if os.path.exists(output_file):
        os.remove(output_file)  # Overwrite the file with an empty DataFrame
        print(f"File '{output_file}' is already exist, so it is cleared first to remove the old contents")
    # read in the files
    try:
        with gzip.open(output_file,"wt",encoding="utf-8") as output:
            # Read the first line (header)
            print('start indexing')
            with gzip.open(input_file, 'rt') as input:
                match_dict = {}
                resulted_index_header = []
                header = input.readline().strip().split('\t')
                for column in header:
                    if column ==  'FINNGENID':
                        continue
                    else:
                        pattern = fr"^{column}_({custom_columns})?$"
                        default_dict = {}
                        default_dict[header[0]] = 0
                        default_dict[COLUMN] = column
                        default_dict[VALUE] = header.index(column)
                        match_dict = {name: idx for idx, name in enumerate(header) if re.match(pattern, name)}
                        if (len(match_dict) > 0):
                            combined_dict = {**default_dict, **match_dict}
                            print('.', end='', flush=True)
                            new_data = {key.replace(f"{column}_", ''): value for key, value in combined_dict.items()}
                            resulted_index_header.append(new_data)
                print('Finished indexing and reading row data!')
                first_keys = resulted_index_header[0].keys()
                is_header_index_matched = all(set(d.keys()) == set(first_keys) for d in resulted_index_header)
                if is_header_index_matched:
                    output_file_header = f"\t".join(first_keys)+"\n"
                    output.write(output_file_header)
                    # read rows in the input data
                    for line in input:
                        row = line.strip().split('\t')
                        print('.', end='', flush=True)
                        map_row_with_header_index(row, resulted_index_header, output)
        print('\nFinished successfully!')
    except Exception as e:  # Catch any exception
        print(f"An error occurred: {e}")
        raise

# check the argument 
parser = argparse.ArgumentParser(description="Convert the wide format into long format!")
parser.add_argument('input_file_path', type=str)
parser.add_argument('output_file_name', type=str)
parser.add_argument('filter_columns', type=str)

args = parser.parse_args()

if f"/{args.input_file_path}".endswith('.tsv.gz'):
    run(args.input_file_path, args.output_file_name, args.filter_columns)
else:
    print("Sorry, output file name should be tsv.gx. e-g string.tsv.gz")
    exit(1)
