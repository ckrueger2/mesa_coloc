#!/usr/bin/python

#import libraries
try:
    import pandas as pd
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install", "pandas"], check=True)
    import pandas as pd
    
import os
import argparse
import sys
import hail as hl
import subprocess

#arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="import hail table from phenotype accession number")
    parser.add_argument('--phecode', help='path to input file', required='True')
    parser.add_argument('--pop', help='path to input file', required='True')
    return parser.parse_args(args)

args = check_arg(sys.argv[1:])

phenotype_id = args.phecode
pop = args.pop

#check if phecode file exists
path = f"gs://fc-aou-datasets-controlled/AllxAll/v1/ht/ACAF/{pop}/phenotype_{phenotype_id}_ACAF_results.ht"

result = os.system(f"gsutil -u $GOOGLE_PROJECT ls {path} > /dev/null 2>&1")

if result == 0:
    print(f"Phenotype {phenotype_id} is in the All of Us database")
else:
    print(f"Phenotype {phenotype_id} is not in the All of Us database; enter valid phenotype ID")
    sys.exit(1)

#PULL TABLE FROM ALL OF US DATABASE   
#initialize hail
hl.init()

#define bucket to save to
bucket = os.getenv('WORKSPACE_BUCKET')
bucket # gs://fc-secure-bb61452f-d5e2-4d26-9227-6a9444241af8/

#find hail table and save to variable
ht = hl.read_table(f"gs://fc-aou-datasets-controlled/AllxAll/v1/ht/ACAF/{pop}/phenotype_{phenotype_id}_ACAF_results.ht")

#columns to keep
desired_columns = ['locus', 'alleles', 'BETA', 'SE', 'Het_Q', 'Pvalue', 'Pvalue_log10', 'CHR', 'POS', 'rank', 'Pvalue_expected', 'Pvalue_expected_log10']

#save the column headers
available_fields = set(ht.row)
key_fields = list(ht.key)
non_key_columns = [col for col in desired_columns if col in available_fields and col not in key_fields]
filtered_ht = ht.select(*non_key_columns)

#add Het_Q column if it doesn't exist
if 'Het_Q' not in available_fields:
    filtered_ht = filtered_ht.annotate(Het_Q=hl.null(hl.tfloat64))

#make sure columns are in the correct order
ordered_columns = [col for col in desired_columns if col in filtered_ht.row]
filtered_ht = filtered_ht.key_by()
ordered_ht = filtered_ht.select(*ordered_columns)

#save full table to bucket for S-PrediXcan input
ht_path = f'{bucket}/data/{pop}_full_{phenotype_id}.tsv'
ordered_ht.export(ht_path)

#CHECK IF FILES ARE SAVED TO BUCKET
try:
    check_filtered = subprocess.check_output(
        f"gsutil ls {bucket}/data/ | grep {ht_path}", 
        shell=True, 
        stderr=subprocess.DEVNULL
    )
    #if command succeeded 
    print("Full file successfully saved to bucket.\n")
except subprocess.CalledProcessError:
    #if command failed
    sys.exit(f"ERROR: File '{ht_path}' was not found in {bucket}/data/.\n")
