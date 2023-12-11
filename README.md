# morph_impact_on_opt_stim
Impact of neuronal morphology on spatial precision of optogenetic stimulation with an optical fiber placed on the cortical surface.

# install miniconda

# create conda env

* run "conda create -n mioos -c bioconda -c conda-forge jupyterlab conda-build matplotlib neuron numpy pandas scipy seaborn snakemake tqdm gcc_linux-64 gxx_linux-64 gfortran_linux-64"

# activate conda env

* run "conda activate mioos"

# download simulation code

* inside repository run "git clone git@github.com:dberling/simneurostim.git"

# compile mod files

* run "nrnivmodl simneurostim/model/mod/"

# add local python package via conda-develop

* run "conda-develop simneurostim/base-neurostim/"

# check out simulation_playground.ipynb for example simulation

# optional: test snakemake execution to run workflows

* run "snakemake --snakefile workflows/Snakefile_simdata

# Errors due to separation of simulation & model code from this particular analysis

Errors which might appear and how to fix:
* old module name "optostim" may appear in code, replace by "neurostim"
* paths of neuron / cell models may be given for the old location ("model/hoc/...") but are now located at "simneurostim/model/hoc/..."
