from neurostim.opt_res_analysis import find_xAPCs_over_light_pwrs
import pandas as pd

data = pd.read_hdf(str(snakemake.input))

result = find_xAPCs_over_light_pwrs(data)

result.to_hdf(str(snakemake.output), key='first')
