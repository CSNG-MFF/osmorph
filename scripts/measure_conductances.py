import numpy as np
import pandas as pd
from neuron import h
from neurostim.cell import Cell
from neurostim.light_classes import LightSource, LightStimulation
from neurostim.utils import convert_polar_to_cartesian_xz, interpolate
from neurostim.polarmaps import get_AP_times

# NEURON setup
h.load_file("stdrun.hoc")
h.cvode_active(1)
# load cell and create stimulation object
cell = Cell(
    hoc_file="simneurostim/model/hoc/" + str(snakemake.wildcards.hoc_file) + ".hoc",
    cortical_depth=snakemake.params.cortical_depth,
    ChR_soma_density = float(snakemake.wildcards.chanrhod_expression),
    ChR_distribution=str(snakemake.wildcards.chanrhod_distribution),
    rm_mech_from_secs=None,
)

# additional recording variables (e.g. conductance)
segs = [seg for sec in h.allsec() for seg in sec][:-1] # exclude last seg as it is the light source
segcatpointernames = [str(seg) for seg in segs]
segcatpointers = [eval('h.'+str(seg)+'._ref_gcat_chanrhod') for seg in segs]

# simulate for all radius and angle combinations
results = []
for radius in np.arange(*list(snakemake.params.radius_min_max_step)):
    for angle in np.arange(*list(snakemake.params.angle_min_max_step)):
        light_x, light_y = convert_polar_to_cartesian_xz(radius, angle)
        light_z = 0  # cortical surface
        light_pos = (light_x, light_y, light_z)
        light_source = LightSource(
            model=str(snakemake.wildcards.light_model),
            position=(light_x, light_y, light_z),
            width=float(snakemake.wildcards.fiber_diameter),
            NA=float(snakemake.wildcards.fiber_NA)
        )
        light_stim = LightStimulation(
            cell=cell,
            light_source=light_source,
            delay=float(snakemake.params.light_delay),
            duration=float(snakemake.wildcards.stim_duration),
            light_power=float(snakemake.wildcards.light_power),
            record_all_segments=False,
        )
        # perform stimulation and record data
        measurement = pd.DataFrame(
            light_stim.simulate_and_measure(
                tot_rec_time=float(snakemake.params.tot_rec_time),
                extra_rec_var_names=segcatpointernames, extra_rec_var_pointers=segcatpointers
            )
        )
        # drop full row duplicates:
        measurement = measurement.drop_duplicates()
        # add 1e-12 ms to 2nd entry of 2 entries with the same time but different values


        measurement.loc[measurement["time [ms]"].diff() == 0, "time [ms]"] += 1e-12
        measurement = interpolate(
            df=measurement, interpolation_dt=float(snakemake.params.interpol_dt)
        )
        AP_times = get_AP_times(
            df=measurement,
            interpol_dt=float(snakemake.params.interpol_dt),
            t_on=float(snakemake.params.light_delay),
            AP_threshold=float(snakemake.params.AP_threshold)
        )
        # convert density conductance into real conductance for each segment
        for seg, seg_gcat in zip(segs, segcatpointernames):
            # eval area of segment and convert from um2 to cm2
            measurement[seg_gcat] *= eval('h.'+str(seg)+'.area()') * 1e-8
        t_sec_g = pd.melt(
                measurement, id_vars=['time [ms]'], value_vars=segcatpointernames
                ).set_index(['variable', 'time [ms]'])
        # take mean of conductance over time for each compartment
        mean_g_sec = t_sec_g.groupby('variable').mean()

        mean_g_sec["hoc_file"]= str(snakemake.wildcards.hoc_file)
        mean_g_sec["chanrhod_distribution"]= str(snakemake.wildcards.chanrhod_distribution)
        mean_g_sec["chanrhod_expression"]= int(np.round(float(snakemake.wildcards.chanrhod_expression),0))
        mean_g_sec["light_model"] = str(snakemake.wildcards.light_model)
        mean_g_sec["fiber_diameter"] = float(snakemake.wildcards.fiber_diameter)
        mean_g_sec["fiber_NA"] = float(snakemake.wildcards.fiber_NA)
        mean_g_sec["light_power"] = float(snakemake.wildcards.light_power)
        mean_g_sec["stim_duration [ms]"] = float(snakemake.wildcards.stim_duration)
        mean_g_sec["radius [um]"] = radius
        mean_g_sec["angle [rad]"] = angle
        mean_g_sec["AP_count"] = len(AP_times)
        mean_g_sec = mean_g_sec.reset_index()
        # convert from Sievert to pS (e-12)
        mean_g_sec['value'] *= 1e12
        mean_g_sec = mean_g_sec.rename(
            columns=dict(variable='segname', value='mean_conductance_pS'))
        results.append(mean_g_sec)
        
# save data
pd.concat(results).set_index(
    [
        "hoc_file",
        "light_model",
        "chanrhod_distribution",
        "chanrhod_expression",
        "fiber_diameter",
        "fiber_NA",
        "stim_duration [ms]",
        "light_power",
        "radius [um]",
        "angle [rad]",
        "segname",
    ]
).to_hdf(str(snakemake.output), key='first')
