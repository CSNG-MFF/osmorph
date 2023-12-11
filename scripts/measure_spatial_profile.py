import numpy as np
import pandas as pd
from neurostim.analysis import simulate_spatial_profile
from neurostim.analysis import get_AP_count
from neurostim.analysis import multiply_seg_with_area


cortical_depth = {
    snakemake.params.neuron_data.loc[snakemake.wildcards.neuron].hoc:\
    snakemake.params.neuron_data.loc[snakemake.wildcards.neuron].cortical_depth,
}
# create variables and analysis
## record time and voltage at soma
seg_rec_vars = [
        ['time [ms]', 'V_soma(0.5)'],
        ['h._ref_t', 'h.soma(0.5)._ref_v']
]           
allseg_rec_var = None
sim_data_transform = None
def analyze_AP_count(sim_data, segs):
    return get_AP_count(
            sim_data, 
            interpol_dt_ms=float(snakemake.params.interpol_dt),
            t_on_ms=float(snakemake.params.light_delay),
            AP_threshold_mV=float(snakemake.params.AP_threshold)
    )
scalar_result_names = ['AP_count']
scalar_result_funcs = [analyze_AP_count]
vector_result_func = None

if bool(snakemake.params.record_conductances):
    ## evaluate ChR conductance in all segments
    allseg_rec_var = '._ref_gcat_chanrhod'
    # transform density conductance into abs conductance for all segments:
    sim_data_transform =  multiply_seg_with_area
    # define analysis on seg conductances:
    def melt(sim_data, segs):
        return pd.melt(
                sim_data, 
                id_vars=['time [ms]'], 
                value_vars = [str(seg) for seg in segs]
        ).set_index(['variable', 'time [ms]'])
    def mean_conductance_seg(sim_data, segs):
        mean_seg_cond_over_time = melt(sim_data, segs).groupby('variable').mean().reset_index()
        mean_seg_cond_over_time['value'] *= 1e12 # convert from S (Sievert) to pS
        mean_seg_cond_over_time = mean_seg_cond_over_time.rename(
                columns=dict(variable='segname', value='mean_conductance_pS')
        )
        return mean_seg_cond_over_time
    vector_result_func = mean_conductance_seg

# run simulation of spatial profile
results = simulate_spatial_profile(
    cell_dict=dict(
        cellname=snakemake.params.neuron_data.loc[
            snakemake.wildcards.neuron].hoc,
        cortical_depth=cortical_depth,
        ChR_soma_density=snakemake.params.neuron_data.loc[
            snakemake.wildcards.neuron].ChR_soma_density,
        ChR_distribution=snakemake.params.neuron_data.loc[
            snakemake.wildcards.neuron].ChR_distribution,
    ),
    stimulator_dict=dict(
        diameter_um=snakemake.params.stimulator_data.loc[
            snakemake.wildcards.stimulator].diameter,
        NA=snakemake.params.stimulator_data.loc[
            snakemake.wildcards.stimulator].NA,

    ),
    stim_intensity_mWPERmm2=float(snakemake.wildcards.intensity),
    radii_um=np.arange(*list(snakemake.params.radius_min_max_step)),
    angles_rad=np.arange(*list(snakemake.params.angle_min_max_step)),
    temp_protocol=dict(
        duration_ms=float(snakemake.params.light_duration),
        delay_ms=float(snakemake.params.light_delay),
        total_rec_time_ms=float(snakemake.params.tot_rec_time),
    ),
    seg_rec_vars=seg_rec_vars,
    allseg_rec_var=allseg_rec_var,
    sim_data_transform=sim_data_transform,
    scalar_result_names=scalar_result_names,
    scalar_result_funcs=scalar_result_funcs,
    vector_result_func=vector_result_func,
    interpol_dt_ms=float(snakemake.params.interpol_dt),
    AP_threshold_mV=float(snakemake.params.AP_threshold),
)
results.to_hdf(str(snakemake.output), key='first')
