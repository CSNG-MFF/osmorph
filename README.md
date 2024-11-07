# morph_impact_on_opt_stim
Impact of neuronal morphology on spatial precision of optogenetic stimulation with an optical fiber placed on the cortical surface.

# software requirements
* jupyterlab
* matplotlib
* neuron
* numpy
* pandas
* scipy
* seaborn
* snakemake
* tqdm
* gcc_linux-64
* gxx_linux-64
* gfortran_linux-64
# download simulation code and install the contained neurostim python package

* inside repository run "git clone git@github.com:dberling/simneurostim.git"
* pip install /simneurostim/base-neurostim

# compile mod files

* run "nrnivmodl simneurostim/model/mod/"

# check out simulation_playground.ipynb for example simulation

# optional: test snakemake execution to run workflows

* run "snakemake --snakefile workflows/Snakefile_simdata
