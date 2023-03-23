from VOI_Functions import *

File = 'Braila_v3.inp'  # File used for the model
network = Network(File) # Creation of the object 
network.GraphNet() # To graph the network 
network.CalculateParameters(Threshold=0.5,Leakmin=4,Leakmax=20,Num_Leakages=20) 
#Main function to calculate pressure variation and other paramater, for more information see VOI_Functions.py

TP =network.TimeDetection 
# This function create a list of arrays
# The number of arrays is equal to the number of nodes where pressure is read (Sensor_Nodes)
# Each Array is composed for dataframes, the number of dataframes is equal to the number of Leaks evaluated
# The rows of the dataframe indicates the hour evalauted (by default 24 rows, 1 day evaluated)
# The columns of each dataframe indicates the nodes where leaks are placed (Leak_Nodes).
# A leak is detected if the change of pressure is bigger than the threshold

fp = network.FlowRequired 
# This function create and matrix where columns are the nodes were pressure are read (Sensor_Node) 
# and rows are nodes where leaks are placed (Leak_Node), each value in the matrix is the minimum leak place in 
# a leak-node to cause a pressure drop bigger than the treshold at the sensor-node

PD = network.Perc_Detection # List with the percentage of cases detected 
# (number of leaks x number of leak-nodes) per each node where pressures are read

Solutions =[]   
PriorB = [] # List of priority nodes, Nodes where VOI is more relevant, by default all nodes have the same VOI

# =============================================================================
# Optimisation
# =============================================================================

A_S = 4 # Number of sensors to be placed
F_Nodes = [] # Fixed nodes, nodes where the user wants to install sensors, by default there is not fixed nodes

if len(F_Nodes) != 0:
    for i in range(len(F_Nodes)):
        network.SensorNodes.pop(network.SensorNodes.index(F_Nodes[i]))      

def Optimal(x):
    y = np.append(x[0],F_Nodes).tolist()    
    O_1 = network.VOI(N=y,L=network.LeakNodes,PriorityNodes=PriorB,VOIMultiplier=50)            
    return [O_1]

problem = Problem(1,1) 
problem.types[:] = Subset(network.SensorNodes,A_S) # define randon nodes in the whole system
problem.function = Optimal # Function created in previous steps.
problem.directions[0] = Problem.MAXIMIZE # VOI has to be maximized, by default Objetives are minimized 
algorithm = NSGAII(problem,population_size=5000) 
algorithm.run(1000) 
Best_S = unique(nondominated(algorithm.result))    
Solutions.append(Best_S)
     
# =============================================================================
# Graphs of Solutions
# =============================================================================
i=1
Soluciones = []
for i in range(len(Solutions)):
    Points = [s.variables[0] for s in Solutions[i]]
    Soluciones.append(Points)
 
for i in Points:   
    print(network.VOI(N=i,L=network.LeakNodes),len(network.NodesDetected),\
          network.Pro_Detection)
    print(network.VOI(N=F_Nodes,L=network.LeakNodes),len(network.NodesDetected),\
      network.Pro_Detection)
        
    
for i in Points:           
    plt.figure()
    ax = plt.gca()
    NodesDe = network.VOI(N=i,L=network.LeakNodes)   
    nn = network.NodesDetected     
    wntr.graphics.plot_network(network.OpenNet(),node_size=10,node_cmap='Pastel2_r',node_attribute=network.SensorNodes,\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax) 
    wntr.graphics.plot_network(network.OpenNet(),node_size=20,node_cmap='hsv',node_attribute=nn,\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax,node_colorbar_label='MLR')
    wntr.graphics.plot_network(network.OpenNet(),node_size=80,node_cmap='hsv',node_attribute=i,\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax)   
    wntr.graphics.plot_network(network.OpenNet(),node_size=35,node_cmap='Wistia',node_attribute=i,\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax)
 



    
    

    