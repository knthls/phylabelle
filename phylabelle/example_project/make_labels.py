"""
This script overwrites the obviously errorous labels with the manually
curated ones, resp. merges these columns to produce a reliable input file.

For an example, how to use the script, have a look at
example_project/sample_call
which provides an example, how to use this script ot reformat your input files,
in order to use it as a label source for your project.
"""

import argparse

from phylabelle.fileio import Table, Header

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('input_file', type=str, nargs=1,
                    help='The original input table')
parser.add_argument('-c', '--curate', metavar=['LABELS', 'CURATED'],
                    type=str, nargs=2,
                    help='overwrite labels in label-column with \
                    labels in curated-column, if possible')
parser.add_argument('-s', '--sep', type=str, default=['\t'], nargs=1,
                    help='The column seperator')
parser.add_argument('fields', type=str, nargs='+',
                    help='The columns which should be parsed')

args = parser.parse_args()

fields = args.fields

if args.sep[0] == '\\t':
    sep = '\t'
else:
    sep = args.sep[0]

with Table(args.input_file[0], 'r', seperator=sep) as tab:
    with Table('stdout', 'w', header=Header('\t', list_=fields)) as target:
        target.write_header()
        for line in tab:
            try:
                if args.curate:
                    cur = line[args.curate[1]]
                    if cur:
                        line[args.curate[0]] = cur
                
                target.write([line[x] for x in args.fields])
            except IndexError:
                pass
