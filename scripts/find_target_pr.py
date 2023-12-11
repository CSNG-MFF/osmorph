import numpy as np
import pandas as pd
from neurostim.analysis import find_target_pr, analyze_AP_count
from neurostim.utils import str_to_lvl2nestedOrderedDict

# snakemake cannot read a dictionary from a list of dictionaries, hence decode for snakemake:
cell_dict = str_to_lvl2nestedOrderedDict(snakemake.wildcards.cell_dict)
stimulator_dict = str_to_lvl2nestedOrderedDict(snakemake.wildcards.stimulator_dict)

print(cell_dict['cortical_depth']['L23'])

# parameters passed to simulate_spatial_profile in analysis.py
sim_spatial_profile_args = dict(
    cell_dict=cell_dict,
    stimulator_dict=stimulator_dict,
    radii_um = np.arange(*list(snakemake.params.radius_min_max_step)),
    angles_rad = np.arange(*list(snakemake.params.angle_min_max_step)),
    temp_protocol=snakemake.params.temp_protocol,
    seg_rec_vars=[
        ['time [ms]', 'V_soma(0.5)'],
        ['h._ref_t',  'h.soma(0.5)._ref_v']
        ],
    allseg_rec_var=None,
    sim_data_transform=None,
    scalar_result_names=['AP_count'],
    scalar_result_funcs=[analyze_AP_count],
    vector_result_func=None,
    interpol_dt_ms=snakemake.params.interpol_dt_ms,
    AP_threshold_mV=snakemake.params.AP_threshold_mV
)

li = find_target_pr(
    target_pr=float(snakemake.wildcards.target_pr), 
    tolerance=float(snakemake.params.tolerance), 
    stim_intensity_mWPERmm2_minmax=\
        snakemake.params.stim_intensity_mWPERmm2_minmax[
            cell_dict['cellname']][cell_dict['ChR_distribution']][cell_dict['ChR_soma_density']],
    simulate_spatial_profile_args=sim_spatial_profile_args
)

results = pd.DataFrame(
        data=[li],
        columns=['li@pr']
        )
results['neuron'] = cell_dict['cellname']
results['ChR_distribution'] = cell_dict['ChR_distribution']
results['ChR_soma_density'] = cell_dict['ChR_soma_density']
results['stim_diameter'] = stimulator_dict['diameter_um']
results['stim_NA'] = stimulator_dict['NA']
results['target_pr'] = snakemake.wildcards.target_pr

results.set_index(
        ['neuron', 'ChR_distribution', 'ChR_soma_density',
         'stim_diameter', 'stim_NA', 'target_pr']
)
results.to_hdf(str(snakemake.output), key='first')
