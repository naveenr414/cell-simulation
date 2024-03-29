from emukit.core import ParameterSpace, ContinuousParameter, DiscreteParameter
from emukit.core.initial_designs import RandomDesign
from emukit.core.initial_designs.latin_design import LatinDesign
from GPy.models import GPRegression
from emukit.model_wrappers import GPyModelWrapper
import GPy
from emukit.sensitivity.monte_carlo import MonteCarloSensitivity
from data_processing.start_experiment import *
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from data_processing.read_data import *
from os.path import exists

def run_simulation_latin(parameter_space, num_samples, starting_identifier, constant_parameters):
    """Run num_samples simulations based on parameters from a LatinDesign
    
    Arguments:
        parameter_space: Object from the Parameter Space class
        num_samples: Number of samples to take from LatinDesign
        starting_identifier: Which number to start on when running simulations
        constant_parameters: Dictionary containing whichever parameters are constant 

    Returns: List of Parameter combinations run
    
    Side Effects: Runs num_samples simulations
    """
    
    parameter_names = [i.name for i in parameter_space.parameters]
    parameter_types = []
    for parameter in parameter_space.parameters:
        if 'Discrete' in str(type(parameter)):
            parameter_types.append(int)
        else:
            parameter_types.append(float)
    
    design = LatinDesign(parameter_space) 
    X = design.get_samples(num_samples)
    
    for i in range(X.shape[0]):
        new_parameters = deepcopy(constant_parameters)
        
        for j,name in enumerate(parameter_names):
            new_parameters[name] = parameter_types[j](X[i][j])
            
        identifier = starting_identifier+i
        config_file = "latin_{}.par".format(identifier)
        
        create_config(config_file,new_parameters)
        execute_experiment(config_file,identifier)

    return X
        
def get_emukit_model(X,Y,lengthscale=1,variance=20,noise_var=1e-10):
    """Get an Emukit Model by fitting X,Y from Gaussian Process Regression
    
    Arguments:
        X: numpy array that represents the parameters
        Y: some output response variable
        lengthscale: Optional argument which is lengthscale for RBF
        variance: Optional argument, variance for RBF
        noise_var: Optional argument for GPRegression
        
    Returns:
        Emukit model based on training a GP model
    """
        
    gpy_model = GPy.models.GPRegression(X, Y, GPy.kern.RBF(X.shape[1], lengthscale=lengthscale, variance=variance), noise_var=noise_var)
    emukit_model = GPyModelWrapper(gpy_model)
    
    return emukit_model

def plot_gaussian_process(index,parameter_space, emukit_model):
    """Plot a Gaussian Process using matplotlib by randomly sampling parameters, and seeing predictions
    Currently uses a LatinDesign to sample parameters
    
    Arguments:
        x_plot: np.linspace that contains which x should be in the plot
        index: Number from 0-len(emukit_model.X); which dimension of the input should we plot
        emukit_model: Which model contains the Gaussian Process

    Returns: Nothing
    
    Side Effects: Creates matplotlib plot
    """
    
    x_plot = np.linspace(emukit_model.X[:,index].min(),emukit_model.X[:,index].max(),100)
    
    design = LatinDesign(parameter_space) 
    X = design.get_samples(len(x_plot))
    X[:,index] = x_plot
    
    x_plot = x_plot.reshape((len(x_plot),1))
    mu_plot, var_plot = emukit_model.predict(X)
    
    mu_plot = mu_plot.flatten().reshape((len(mu_plot),1))
    var_plot = var_plot.flatten().reshape((len(var_plot),1))
    plt.plot(emukit_model.X[:,index], emukit_model.Y, "ro", markersize=10, label="Observations")
    plt.plot(x_plot, mu_plot, "C0", label="Model")
    plt.fill_between(x_plot[:, 0],
                     mu_plot[:, 0] + np.sqrt(var_plot)[:, 0],
                     mu_plot[:, 0] - np.sqrt(var_plot)[:, 0], color="C0", alpha=0.6)
    plt.fill_between(x_plot[:, 0],
                     mu_plot[:, 0] + 2 * np.sqrt(var_plot)[:, 0],
                     mu_plot[:, 0] - 2 * np.sqrt(var_plot)[:, 0], color="C0", alpha=0.4)
    plt.fill_between(x_plot[:, 0],
                     mu_plot[:, 0] + 3 * np.sqrt(var_plot)[:, 0],
                     mu_plot[:, 0] - 3 * np.sqrt(var_plot)[:, 0], color="C0", alpha=0.2)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$f(x)$")
    plt.xlim(min(x_plot),max(x_plot))
    plt.legend()
    plt.grid(True)

def get_sobol_indices(emukit_model,parameter_space):
    """Perform sensitivity analysis by retrieving the sobol indices from an emukit model and parameter space
    
    Arguments:
        emukit_model: Gaussian Process model trained on simulation data
        parameter_space: ParameterSpace object with information on each parameter
        
    Returns: 
        Two things: 
            main_effects: Dictionary with the main effects indices for each parameter
            total_effects: Dictionary with the total effects indices for each parameter
    """
    
    senstivity = MonteCarloSensitivity(model = emukit_model, input_domain = parameter_space)
    main_effects, total_effects, _ = senstivity.compute_effects(num_monte_carlo_points = 10000)
    
    return main_effects, total_effects


def get_timestep_rewards(simulation_f, reward_function):
    """Get rewards from a simulation, according to some reward_function
    
    Arguments:
        simulation: a file name representing the output of a single simulation
        reward_function: Function that takes in a list of cells and returns a list of rewards
        
    Returns:
        Numpy array of rewards
    """
    if not exists("../data/output/{}".format(simulation_f)):
        return np.array([])
    
    all_cells = read_data(simulation_f)

    rewards = reward_function(all_cells)
    return np.array(rewards)
    
def get_rewards(simulation_list,reward_function):
    """Get rewards from a list of simulations, according to some reward_function
    
    Arguments:
        simulation_list: String list of file names that represent different simulations
        reward_function: Function that takes in a list of cells and returns some float
        
    Returns:
        Numpy array of rewards
    """
    
    rewards = np.empty(shape=[0, 1])
    
    for i in range(len(simulation_list)):
        if not exists("../data/output/{}".format(simulation_list[i])):
            continue
            
        all_cells = read_data(simulation_list[i])
        
        rewards = np.vstack([rewards, reward_function(all_cells)])
        
    return rewards


if(__name__ == "__main__"):
    """
    run latin hypercube sampling with configuration below
    """

    parameters = [
        DiscreteParameter("season_duration", range(100, 100001)),
        ContinuousParameter("mut_rate", 0, 0.5),

        DiscreteParameter("T", range(1, 101)),
        DiscreteParameter("target_area", range(0, 101)),
        ContinuousParameter("gradnoise", 0, 1),

    ]
        



        

