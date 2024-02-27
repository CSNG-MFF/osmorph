from neurostim.opt_res_analysis import improved_find_xAPC50_over_lps
import pandas as pd

data = pd.read_hdf(str(snakemake.input))

result = improved_find_xAPC50_over_lps(data)

result.to_hdf(str(snakemake.output), key='first')
