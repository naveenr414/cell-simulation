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

if __name__ == "__main__":
    all_cells = read_data("data_cellcount_testing.txt")
    print(len(all_cells))
