from pathlib import Path
import openpyxl,json
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/sources'
for fn in ['ICPA_accepted_species_20261001.xlsx','ICPA_reference_sequences_20261001.xlsx']:
    w=openpyxl.load_workbook(P/fn,read_only=True,data_only=True)
    print('WORKBOOK',fn)
    for s in w:
        print('SHEET',s.title,s.max_row,s.max_column)
        for i,row in enumerate(s.iter_rows(values_only=True),1):
            if i<4 or any(str(v)=='Terrei' for v in row):
                print(i,[str(v)[:190] if v is not None else '' for v in row])
