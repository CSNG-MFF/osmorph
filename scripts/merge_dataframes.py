import pandas as pd

dfs = []
for fname in list(snakemake.input):
    df = pd.read_hdf(fname)
    df['fname']=fname
    dfs.append(df)
pd.concat(dfs).to_hdf(str(snakemake.output), key='first')
