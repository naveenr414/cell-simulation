import os

def execute_experiment(config_file,identifier):
    """Run an experiment based on a config file

    Arguments:
    config_file: .par file located in the data/parameters folder
    identifier: Unique number defining the experiment

    Returns: Nothing

    Side Effects: Runs the simulation and creates the corresponding data
    Copies the config and output file to the output folder, with the same name
    """
    
    output_file = "data_cellcount_{}.txt".format(identifier)
    backup_folder = "backup_{}".format(identifier)
    images_folder = "data_film{}".format(identifier)
    
    os.system("bash bash_scripts/run_simulation.sh {} {} &".format(config_file,identifier))

def create_config(file_name,data_dictionary):
    """Create a new config file based on updating parameters from a baseline parameter set

    Arguments:
    file_name: .par file which will be stored in the data/parameters folder
    data_dictionary: Dictionary, with keys being parameter names, and values being values

    Returns: Nothing

    Side Effects: Writes new file_name at data/parameters
    """

    baseline_parameters = {}
    baseline_file = "../data/parameters/{}".format("baseline_parameters.par")
    f = open(baseline_file).read().strip().split("\n")

    for line in f:
        line = line.replace("\r\n","\n")
        split_line = line.split(" = ")
        baseline_parameters[split_line[0].strip()] = split_line[1].strip()

    parameters = baseline_parameters

    for i in data_dictionary:
        if i == 'gamma':
            parameters['keylock_list_filename'] = '../data/keylock/KL_same_{}.dat'.format(data_dictionary[i])
        else:
            parameters[i] = data_dictionary[i]

    w = open("../data/parameters/{}".format(file_name),"w")
    for i in parameters:
        w.write("{} = {}\n".format(i,parameters[i]))
    w.close()
    
if __name__ == "__main__":
    file_name = "basic_chemotax2.par"
    output_file = "data_cellcount_testing.txt"
    
    mcs = 10
    season_experiment = 1
    season_duration = 10000

    create_config(file_name, {'mcs': mcs,
                              'season_experiment': season_experiment,
                              'season_duration': season_duration})
    execute_experiment(file_name,output_file)
