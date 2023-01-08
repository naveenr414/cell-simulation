import matplotlib.pyplot as plt
import numpy as np

class Cell:
    def __init__(self):
        self.time = 0
        self.identifier = 0
        self.cell_type = 0
        self.pos = [0,0]
        self.persistent_migration = [0,0]
        self.chemotaxis = [0,0]
        self.time_since_birth = 0
        self.key = ""
        self.lock = ""
        self.mu = 0
        self.num_particles = 0

        self.maintenance_k = {'fraction': 0, 'k_0': 0, 'k_A': 0,
                              'k_P': 0, 'k_C': 0}
        self.ext_k = {'fraction': 0, 'k_0': 0, 'k_A': 0, 'k_P': 0,
                      'k_C': 0, 'k_0t': 0, 'k_Pt': 0}
        self.chem_k = {'fraction': 0, 'k_0': 0, 'k_A': 0, 'k_P': 0,
                       'k_C': 0}

        self.neighbors = []
        
def read_parameters(parameter_file):
    """Function to return parameters in a dictionary from a given file, located in data/parameters
    
    Arguments:
        parameter_file: String for the parameter file name
        
    Returns:
        Dictionary representing the values of each of the parameters
    """
    
    baseline_parameters = {}
    baseline_file = "../data/parameters/{}".format(parameter_file)
    f = open(baseline_file).read().strip().split("\n")

    for line in f:
        line = line.replace("\r\n","\n")
        split_line = line.split(" = ")
        baseline_parameters[split_line[0].strip()] = split_line[1].strip()

    return baseline_parameters

def get_parameter_array(file_names,parameter_names):
    """Get each parameter from each parameter_file name, returned as a numpy array
    
    Arguments:
        file_names: List of file name strings, with each file in data/parameters
        parameter_names: List of parameter name strings
        
    Returns:
        Numpy data array of size (len(file_names),len(parameter_names))
    """
    
    parameter_array = np.zeros((len(file_names),len(parameter_names)))
    
    for i,file in enumerate(file_names):
        parameters = read_parameters(file)
        
        for j,p in enumerate(parameter_names):
            if p == 'gamma':
                if 'KL_same' not in parameters['keylock_list_filename']:
                    raise Exception("To use parameter gamma, keylock_list_filename must be in format KL_same_numbe")
                
                value = parameters['keylock_list_filename'].replace(".dat","")
                value = float(value.split("_")[-1])
            elif p in parameter_names:
                value = float(parameters[p])
            else:
                raise Exception("Parameter {} not a valid parameter".format(p))

            parameter_array[i][j] = value
            
    return parameter_array
            

def read_data(cell_file):
    """Read data from the cell file, located in data/output/cell_file

    Arguments:
    cell_file: String representing some data_cellcount.txt file

    Returns: List of Cell Objects, with processed information from each line
    """

    f = open("../data/output/{}".format(cell_file)).read().strip().split("\n")
    all_cells = []

    for i in f:
        line = i.strip().split(" ")
        cell = Cell()
        cell.time = int(line[0])
        cell.identifier = int(line[1])
        cell.cell_type = int(line[2])
        cell.pos = [float(line[3]),float(line[4])]
        cell.persistent_migration = [float(line[5]),float(line[6])]
        cell.chemotaxis = [float(line[7]),float(line[8])]
        cell.time_since_birth = int(line[9])
        cell.key = line[10]
        cell.lock = line[11]
        cell.mu = float(line[12])
        cell.num_particles = int(line[13])

        cell.maintenance_k = {'fraction': float(line[14]),
                              'k_0': float(line[15]),
                              'k_A': float(line[16]),
                              'k_P': float(line[17]),
                              'k_C': float(line[18])}
        cell.ext_k = {'fraction': float(line[19]),
                      'k_0': float(line[20]),
                      'k_A': float(line[21]),
                      'k_P': float(line[22]),
                      'k_C': float(line[23]),
                      'k_0t': float(line[24]),
                      'k_Pt': float(line[25])}
        cell.chem_k = {'fraction': float(line[26]),
                       'k_0': float(line[27]),
                       'k_A': float(line[28]),
                       'k_P': float(line[29]),
                       'k_C': float(line[30])}

        for i in range(31,len(line),2):
            cell_type = int(line[i])
            energy_difference = float(line[i+1])
            cell.neighbors.append({'cell_type': cell_type,
                                   'energy_difference': energy_difference})

        
        all_cells.append(cell)

    return all_cells

def cell_dist(cell_a,cell_b):
    """Find the distance between two cells
    
    Arguments:
    cell_a: Object from the cell class
    cell_b: Object from the cell class
    
    Returns:
    Float, representing the distance between cell_a and cell_b
    """
    
    pos_a = np.array(cell_a.pos)
    pos_b = np.array(cell_b.pos)
    
    return np.linalg.norm(pos_a-pos_b)

def create_graph(cell_data, time_period):
    """Create an adjacency matrix based on a list of objects from Cell
    
    Arguments:
    cell_data: List of objects from cell 
    time_period: Which time period to take cells from 
    
    Returns: 
    numpy matrix, 0-1 adjacency matrix
    """
    
    cells_at_time = [i for i in cell_data if i.time == time_period] 
    num_cells = len(cells_at_time)
    
    adjacency_matrix = np.zeros((num_cells,num_cells))
    
    for i,cell in enumerate(cells_at_time):
        neighbors = cell.neighbors
        baseline_energy_difference = [i for i in neighbors if i['cell_type'] == 0]
        
        if len(baseline_energy_difference) == 0:
            continue
        
        baseline_energy_difference = baseline_energy_difference[0]['energy_difference']
        
        adhering_neighbors = [i for i in neighbors if i['cell_type'] == 1 and                               baseline_energy_difference-i['energy_difference']/2 > 0]
        
        num_neighbors = len(adhering_neighbors)
        previous_neighbors = sum(adjacency_matrix[i])
        num_new_neighbors = num_neighbors-previous_neighbors
        num_new_neighbors = max(num_new_neighbors,0)
        
        all_distances = [(j,cell_dist(cell,cell_b)) for j,cell_b in enumerate(cells_at_time) if j != i]
        
        all_distances = sorted(all_distances, key=lambda k: k[1],reverse=True)
        
        while num_new_neighbors > 0:
            cell_num,closest_cell = all_distances.pop()
            
            if adjacency_matrix[i][cell_num] == 0:
                adjacency_matrix[i][cell_num] = 1
                adjacency_matrix[cell_num][i] = 1
                num_new_neighbors -= 1

    return adjacency_matrix 

def merge_blobs(blob_list, c1, c2):
    """merges the blobs at c1 and c2 in the blob list (updates the blob representation for all its members)
    
    Arguments:
    blob_list : List[Set[Int]] contains sets that represent the blob of the cell and any given index
    c1 : index of cell1
    c2 : index of cell2
    """
    if(c2 in blob_list[c1]):
        return

    result = blob_list[c1].union(blob_list[c2])
    for i in blob_list[c1]:
        blob_list[i] = result
    for i in blob_list[c2]:
        blob_list[i] = result


def find_blobs(adjacency_matrix):
    """extracts blobs from adjacency matrix

    Arguments
    adjacency_matrix : np.array (2D)
    
    Returns
    List of blobs represented as lists of ints: List[List[Int]]
    """
    blob_list = [set([i]) for i in range(len(adjacency_matrix))]

    for i in range(len(adjacency_matrix)):
        for j in range(i+1, len(adjacency_matrix)):
            if(adjacency_matrix[i,j]):
                merge_blobs(blob_list, i, j)
    
    blob_list = set([" ".join([str(i) for i in sorted(list(blob))]) for blob in blob_list])
    blob_list = [[int(i) for i in blob.split(" ")] for blob in blob_list]
    return blob_list

def weighted_average_blobsize(blob_list):
    n_cells = 0
    result = 0
    for blob in blob_list:
        n_cells += len(blob)
        result += len(blob)**2
    return result/n_cells

def weighted_blobsize_cells(all_cells,t):
    """Given a list of cells, find their weighted average blob size at time t

    Arguments:
        all_cells: A list of cell objects
        t: An integer representing the time

    Returns:
        The weighted average blob size at time t
    """
    
    return weighted_average_blobsize(find_blobs(create_graph(all_cells,t)))
    
def num_blobs_cells(all_cells,t):
    """Given a list of cells, find the number of blobs at time t

    Arguments:
        all_cells: A list of cell objects
        t: An integer representing the time

    Returns:
        The number of blobs at time t
    """

    return len(find_blobs(create_graph(all_cells,t)))

if __name__ == "__main__":
    all_cells = read_data("data_cellcount_testing.txt")
    print(len(all_cells))
