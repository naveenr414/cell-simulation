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
        
        

cell_count_location = "../../src/data_cellcount.txt"

f = open(cell_count_location).read().strip().split("\n")
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
