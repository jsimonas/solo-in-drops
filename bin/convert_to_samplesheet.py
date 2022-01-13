#!/usr/bin/env python3
import sys
import openpyxl
import argparse
import pandas as pd
import numpy as np

# script to covert extended sample sheet to standard illumina's sheet
# template of the extended_sample_sheet.xlsx can be found here:
# solo-in-drops/assets/extended_sample_sheet_template.xlsx

def convert_to_samplesheet(in_file, out_file):
    
    # read extended sample sheet
    ex_sheet = pd.read_excel(
        io = in_file,
        engine='openpyxl'
    )
    # validate if mandatory variables are not missing
    if not ex_sheet.loc[:, ['project_id', 'sample_id', 'index_seq']].isnull().any().any():
       
        # hardcode header of sample sheet
        header_dict={
            '0': ['[Header]', 'IEMFileVersion', 'Investigator Name', 'Experiment Name',
                  'Date', 'Workflow', 'Application', 'Assay', 'Description', 'Chemistry',
                  '', '[Reads]', '', '', '', '[Settings]', '', '[Data]', 'Sample_ID'],
            '1': np.repeat(['', 'Sample_Name'], [18,1]).tolist(),
            '2': np.repeat(['', 'Sample_Plate'],[18,1]).tolist(),
            '3': np.repeat(['', 'Sample_Well'], [18,1]).tolist(),
            '4': np.repeat(['', 'I7_Index_ID'], [18,1]).tolist(),
            '5': np.repeat(['', 'index'], [18,1]).tolist(),
            '6': np.repeat(['', 'Sample_Project'], [18,1]).tolist(),
            '7': np.repeat(['', 'Description'], [18,1]).tolist(),
            '8': np.repeat('', [19]).tolist(),
            '9': np.repeat('', [19]).tolist()
        }
        header_df = pd.DataFrame(header_dict)
        
        # add contents of extended sheet
        ex_sheet_dict = {
            '0': ex_sheet['sample_id'].tolist(),
            '5': ex_sheet['index_seq'].tolist(),
            '6': ex_sheet['project_id'].tolist()
        }
        sheet_df = pd.DataFrame(ex_sheet_dict)
        
        # combine and write sample sheet
        sheet = header_df.append(sheet_df).fillna('')
        sheet.to_csv(
            out_file,
            header = None,
            index = None
        )
    else:
        sys.exit('error: check if there are no empty entries in the extended_sample_sheet.xlsx')


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--file', type=str, metavar='FILENAME',
                        help='path to extended_sample_sheet.xlsx, template can be found in solo-in-drops/assets/')
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='path to output of the standard sample_sheet.csv')
    
    args = parser.parse_args()

    convert_to_samplesheet(
        in_file = args.file,
        out_file = args.out
    )

if __name__ == "__main__":
    main()
