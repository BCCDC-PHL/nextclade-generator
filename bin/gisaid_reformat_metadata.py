import os, sys 
import pandas as pd 

df = pd.read_csv(sys.argv[1], sep='\t')

df[['type','subtype']] = df['subtype'].str.split("_/_", expand=True)

df = df['name strain type subtype clade'.split()]

df.to_csv(sys.argv[2], index=False, sep='\t')
