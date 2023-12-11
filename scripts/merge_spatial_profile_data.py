import pandas as pd

data = pd.concat([pd.read_hdf(fname, key='first') for fname in list(snakemake.input)])

data.to_hdf(str(snakemake.output), key='first')
