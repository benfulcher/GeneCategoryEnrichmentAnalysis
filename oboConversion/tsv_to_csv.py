# Converts .txt files generated using https://github.com/sgrote/OboToTerm to csvs for easy importing into tables in an sqlite db.

import pandas as pd
import os

files = ['term.txt', 'term2term.txt']

for file in files:
    tsv = pd.read_csv(file, sep='\t', header=None)
    name = os.path.basename(file).split('.')[0]
    tsv.to_csv(f'{name}.csv', sep=',', index=False, header=False)
