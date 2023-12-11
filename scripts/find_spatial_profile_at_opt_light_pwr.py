import numpy as np
import pandas as pd
from neuron import h
from neurostim.cell import Cell
from neurostim.light_classes import LightSource, LightStimulation
from neurostim.utils import convert_polar_to_cartesian_xz, interpolate
from neurostim.polarmaps import get_AP_times
from neurostim.opt_res_analysis import *
from neurostim.old_spatial_res_analysis import find_APmax_50_10

# NEURON setup
h.load_file("stdrun.hoc")
h.cvode_active(1)
# load cell and create stimulation object
cell = Cell(
    hoc_file="simneurostim/model/hoc/" + str(snakemake.wildcards.hoc_file) + ".hoc",
    cortical_depth=snakemake.params.cortical_depth,
    ChR_soma_density=float(snakemake.wildcards.chanrhod_expression),
    ChR_distribution=str(snakemake.wildcards.chanrhod_distribution),
    rm_mech_from_secs=None,
)
def sim_all_radii_and_angles(light_power, radii):
    # simulate for all radius and angle combinations
    results = []
    for radius in radii:
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
                light_power=float(light_power),
                record_all_segments=False,
            )
            # perform stimulation and record data
            measurement = pd.DataFrame(
                light_stim.simulate_and_measure(
                    tot_rec_time=float(snakemake.params.tot_rec_time)
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
            result = {
                "hoc_file": str(snakemake.wildcards.hoc_file),
                "chanrhod_distribution": str(snakemake.wildcards.chanrhod_distribution),
                "chanrhod_expression": int(np.round(float(snakemake.wildcards.chanrhod_expression),0)),
                "light_model": str(snakemake.wildcards.light_model),
                "fiber_diameter": float(snakemake.wildcards.fiber_diameter),
                "fiber_NA": float(snakemake.wildcards.fiber_NA),
                "light_power": float(light_power),
                "stim_duration [ms]": float(snakemake.wildcards.stim_duration),
                "radius [um]": radius,
                "angle [rad]": angle,
                "AP_count": len(AP_times),
                "firint_rate [Hz]": len(AP_times) / float(snakemake.wildcards.stim_duration) * 1000, 
                # *1000 to convert stim_duration from ms to s
            }
            if len(AP_times) == 0:
                result["time-to-spike [ms]"] = np.nan
            else:
                result["time-to-spike [ms]"] = AP_times.iloc[0]

            results.append(result)
    # concatenate data
    results = pd.DataFrame(results).set_index(
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
        ]
    )
    # average over angles
    angle_averaged_df = pd.DataFrame(avrg_angles(results))
    profiles_df = find_APmax_50_10(angle_averaged_df, longform=False)
    assert len(profiles_df) == 1, "length of profiles df is not 1"
    APCmax = profiles_df['AP_max'].values[0]
    return APCmax, results

APC_desired = float(snakemake.wildcards.APC_desired)
lmin = float(snakemake.wildcards.l_min)
lmax = float(snakemake.wildcards.l_max)
cnt = 0
rad_min, rad_max, rad_step = list(snakemake.params.radius_min_max_step)
test_radii = np.arange(rad_min, 200+rad_step,rad_step)
other_radii = np.arange(200+rad_step, rad_max, rad_step)
while True:
    l_test = np.exp((np.log(lmin)+np.log(lmax))/2)
    APCmax, results = sim_all_radii_and_angles(light_power=l_test, radii=test_radii)
    print('----')
    print('tested l: ',l_test, 'APCmax: ',APCmax)
    cnt += 1
    if np.round(APCmax,0) == APC_desired:
        break
    elif np.round(APCmax,0) > APC_desired:
        lmax = l_test
    elif np.round(APCmax,0) < APC_desired:
        lmin = l_test
    if cnt > 30:
        results.iloc[:,:]='failed'
        break

APCmax, results_other_radii = sim_all_radii_and_angles(light_power=l_test, radii=other_radii)
results = results.append(results_other_radii)
# annotate APC_desired
results['APC_desired']=APC_desired
results.reset_index().set_index(
        [
            "hoc_file",
            "light_model",
            "chanrhod_distribution",
            "chanrhod_expression",
            "fiber_diameter",
            "fiber_NA",
            "stim_duration [ms]",
            "APC_desired",
            "radius [um]",
            "angle [rad]",
        ]
).to_hdf(str(snakemake.output), key='first')
