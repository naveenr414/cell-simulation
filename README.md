# Cell Simulation
This work is based on the Cell_Evolution_stickymoves repository (https://github.com/escolizzi/Cell_Evolution_stickymoves) which contains the experiments of the paper "Evolution of multicellularity by collective integration of spatial information". We enhance usability of the original code and provide wrapper functions that enable execution of experiments in a Python environment. Furthermore, we implement surrogate models to investigate several questions regarding the properties of the original simulation.

# Running the Simulation

Note that all scripts must be run from the scripts folder. To execute a single simulation run use the command 
```bash 
sh bash_scripts/run_simulation.sh -<config> -<identifier>
```
Parameters:

* ```<config>``` is the name of the configuration file as located in /data/parameters containing the simulation parameters.
* ```<identifier>``` is an integer used to identify the simulation's outputs. After the run is finished they will be stored in /data/output/data_cellcount_\<identifier\>.txt

# Python Wrapper Functions
Experiments can be set up and run from Python using the /scripts/dataprocessing/start_experiment.py module. It enables the creation of a config file based on updating parameters from a baseline parameter set read from baseline_parameters.par. Furthermore, provides an interface for running the simulation given a config file. 

## Post-processing of Simulation Results
The /scripts/dataprocessing/read_data.py data module allows to read in the output of an experimental run and to extract many of properties of the results. These include the number of multi-cellular structures ("blobs") at each time step and their sizes. 

# Building surrogate models
Methods realting to surrogate modelling are implemented in the /scripts/dataprocessing/gaussian.py module. These include methods for running latin hypercube sampling on the simulation, the fitting of a gaussian process (GP) wrapped by an emukit model, the visualization of a GP, and calculation of the sobol indices for sensitivity analysis. 

In the data folder we provide the parameters and results of almost 300 simulation runs. The purpose of this is to empower users to initialize and build their own surrogate models.

Tutorials are provided in the form of jupyter notebooks in the /scripts/notebooks folder. Furthermore, it contains the source code of several surrogate models as well as their analysis that enables a deeper understanding of several properties of the original simulation.

# Details on the Simulation 
* Cells move at 0.0055 px/time step. 
* The simulation is run based on 81 parameters. The subsection below, gives a brief overview over the most important ones.

## Parameters

* **mcs** : "monte carlo steps" defines the maximum time steps for the simulation; takes about 1hour per 100000 timesteps  
* **season_duration** : how many time steps a season lasts
* **gridx** and **gridy** : * define the grid size (default: 500 * 500)
* **mut_rate** : mutation rate of a single gene in the key or lock of a cell; typical key and lock size is 10 + 10, therefore, a mut_rate=0.025 leads to 20 * 0.025 = 0.5 an expected one mutation every two time steps 
* **gradnoise** : controls the noise of the (food source) gradient detection by individual cells (default: 0.9); deterministic if 0.0
