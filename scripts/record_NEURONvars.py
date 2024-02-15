from neuron import h
from neurostim.cell import Cell
from neurostim.models import *
from neurostim.stimulator import Stimulator
from neurostim.simulation import SimControl
from neurostim.utils import convert_polar_to_cartesian_xz
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# run simulation 
cell_dict=dict(
    cellmodel=snakemake.params.neuron_data.loc[
        snakemake.wildcards.neuron].model,
    ChR_soma_density=snakemake.params.neuron_data.loc[
        snakemake.wildcards.neuron].ChR_soma_density,
    ChR_distribution=snakemake.params.neuron_data.loc[
        snakemake.wildcards.neuron].ChR_distribution,
)
stimulator_dict=dict(
    diameter_um=snakemake.params.stimulator_data.loc[
        snakemake.wildcards.stimulator].diameter,
    NA=snakemake.params.stimulator_data.loc[
        snakemake.wildcards.stimulator].NA,

)
stim_intensity_mWPERmm2=float(snakemake.wildcards.intensity)
radii_um=np.arange(*list(snakemake.params.radius_min_max_step))
angles_rad=np.arange(*list(snakemake.params.angle_min_max_step))
temp_protocol=dict(
    duration_ms=float(snakemake.params.light_duration),
    delay_ms=float(snakemake.params.light_delay),
    total_rec_time_ms=float(snakemake.params.tot_rec_time),
)
seg_rec_vars=snakemake.params.rec_vars
interpol_dt_ms=float(snakemake.params.interpol_dt)

# NEURON setup
h.load_file("stdrun.hoc")
h.cvode_active(1)
# load cell
cell = Cell(
    model=eval(cell_dict['cellmodel']),
    ChR_soma_density = cell_dict['ChR_soma_density'],
    ChR_distribution=cell_dict['ChR_distribution'],
    rm_mech_from_secs=None,
)
# init stimulator
stimulator = Stimulator(
    diameter_um=stimulator_dict['diameter_um'],
    NA=stimulator_dict['NA'],
)
results = []
for radius in radii_um:
    for angle in angles_rad:
        # stim position in cartesian coords
        stim_x_um, stim_y_um = convert_polar_to_cartesian_xz(radius, angle)
        stim_z_um = 0  # cortical surface
        # init simulation
        simcontrol = SimControl(
            cell=cell,
            stimulator=stimulator
        )
        # run simulation
        sim_data = simcontrol.run(
            temp_protocol=temp_protocol,
            stim_location=(stim_x_um, stim_y_um, stim_z_um),
            stim_intensity_mWPERmm2=stim_intensity_mWPERmm2,
            rec_vars=seg_rec_vars,
            interpol_dt_ms=interpol_dt_ms
        )
        sim_data['radius'] = radius
        sim_data['angle'] = angle
        results.append(sim_data)
pd.concat(results).to_hdf(str(snakemake.output), key='first')
