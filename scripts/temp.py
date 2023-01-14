from data_processing.start_experiment import *
from data_processing.read_data import *
from data_processing.gaussian import *

import matplotlib.pyplot as plt
import numpy as np
from emukit.core import ParameterSpace, ContinuousParameter, DiscreteParameter
from emukit.core.initial_designs import RandomDesign
from emukit.core.initial_designs.latin_design import LatinDesign
from GPy.models import GPRegression
from emukit.model_wrappers import GPyModelWrapper
import GPy
from emukit.sensitivity.monte_carlo import MonteCarloSensitivity
import time

num_experiments = 10
parameter_list = ['season_duration']
constant_parameters = {
    "mcs" : 50000,
    "cell_size" : 50,
    "gamma" : 5,
    "T" : 16,
    "storage_stride":10000
}
parameter_space = ParameterSpace([DiscreteParameter('season_duration',list(range(1000, 10001, 25)))])
run_simulation_latin(parameter_space, num_experiments, 12, constant_parameters)
