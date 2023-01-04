import os

#if run in other directory: os.chdir("<your_path>\\Cell_Evolution_stickymoves\\src")
#assuming compiled ./cell_evolution

src_dir = os.getcwd()

cpp_shell = "wsl" #whatever shell you use on your os to execute cpp 

def execute_experiment(config_file):
    """Run an experiment based on a config file

    Arguments:
    config_file: .par file located in the data/parameters folder

    Returns: Nothing

    Side Effects: Runs the simulation and creates the corresponding data
    """
    
    os.system("bash bash_scripts/run_simulation.sh {}".format(config_file)) 

#creates experiment config
def create_config(file_name, mcs, season_experiment, season_duration):    
    text = """
mcs = {mcs}
season_experiment = {season_experiment}
season_duration = {season_duration}
sizex = 500
sizey = 500
scaling_cell_to_ca_time = 10
T = 16
target_area = 50
half_div_area = 500
half_div_area_2 = -1
target_length = 0
lambda = 5
lambda2 = 0
startmu = 3.
persduration = 50
init_chemmu = 1.0
zero_persistence_past_theline = 0
keylock_list_filename = ../data/KL_l24_14_16_g6.dat
Jmed_rule_input = 8o4_3_2_1_1_1
key_lock_length = 24
howmany_makeit_for_nextgen = 50
popsize = 100
n_init_cells = 100
size_init_cells = 25
the_line = 50
divisions = 0
conn_diss = 0
vecadherinknockout = false
chemotaxis = 0
extensiononly = false
border_energy = 100
min_area_for_life = 4
neighbours = 2
periodic_boundaries = false
init_cell_config = 0
mut_rate = 0.0
evolsim = 1
evolreg = 0
n_chem = 0
food_influx_location = specified_experiment
gradnoise = 0.9
gradscale = 1.0
is_there_food = 0
initial_food_amount = 0
foodinflux = 0.05
eatprob = 1.0
growth = 0.
ardecay = 0.0
min_contact_duration_for_preying = 25
frac_contlen_eaten = 1.
metabolic_conversion = 0.5
init_maintenance_fraction = 1.
init_k_mf_0 = 1.
init_k_mf_A = 0.
init_k_mf_P = 0.
init_k_mf_C = 0.
init_k_ext_0 = 1.
init_k_ext_A = 0.
init_k_ext_P = 0.
init_k_ext_C = 0.
init_k_ext_0t = 0.
init_k_ext_Pt = 0.
init_k_chem_0 = 1.
init_k_chem_A = 0.
init_k_chem_P = 0.
init_k_chem_C = 0.
init_weight_for_chemotaxis = 1.0
rseed = -1
chancemediumcopied = 0.0001
subfield = 1.0
relaxation = 0
storage_stride = 1000
graphics = false
store = true
datadir = data_film2
save_text_file_period = 1000
datafile = data_cellcount.txt
save_backup_period = 100000
backupdir = backup
readcolortable = false
colortable_filename = ../data/circular.ctb    
    """.format(mcs=mcs, season_experiment=season_experiment, season_duration=season_duration)
    file_path = "../data/parameters/{}".format(file_name)
    with open(file_path, "w", newline = "\n") as f:
        f.write(text)
    
if __name__ == "__main__":
    file_name = "basic_chemotax2.par"
    mcs = 1000
    season_experiment = 1
    season_duration = 10000

    create_config(file_name, mcs, season_experiment, season_duration)
    
    execute_experiment("basic_chemotax.par")
