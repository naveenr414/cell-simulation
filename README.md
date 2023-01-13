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
Experiments can be set up and run from Python using the start_experiment module. It enables the creation of a config file based on updating parameters from a baseline parameter set read from baseline_parameters.par. Furthermore, provides an interface for running the simulation given a config file. 

# Post-processing of Simulation Results
The read_data data module allows to read in the output of an experimental run and to extract many of properties of the results. These include the number of multi-cellular structures ("blobs") at each time step and their sizes. 



# Building surrogate models
notebook sensetivity a b c 

In the data folder we provide the parameters and results of almost 300 simulation runs. The purpose of this is to empower users to initialize and build their own surrogate models.