from Functions import *

File = 'Braila_v3.inp'  # File used for the model
network = Network(File) # Creation of the object 
network.GraphNet() # To graph the network 


network.CalculateParameters(Threshold=1,Leakmin=8.5,Leakmax=20) #Main fucntion
# to calculate pressure variation and other paramaters once leaks are added

# =============================================================================
# Optimisation
# =============================================================================

def f(X):
    X =X.tolist()
    S = []
    for i in X:
        i = int(i)
        S.append(network.SensorNodes[i])      
    O_1 = network.VOI(N=S,L=network.LeakNodes)    
    return -np.sum(O_1)


varbound=np.array([[0,len(network.SensorNodes)-1]]*4)

algorithm_param = {'max_num_iteration': 1000,\
                   'population_size':100,\
                   'mutation_probability':0.1,\
                   'elit_ratio': 0.01,\
                   'crossover_probability': 0.5,\
                   'parents_portion': 0.3,\
                   'crossover_type':'uniform',\
                   'max_iteration_without_improv':15}


model=ga(function=f,dimension=4,variable_type='int',variable_boundaries=varbound,\
         algorithm_parameters=algorithm_param)

model.run()

# =============================================================================
# Graph of Solution
# =============================================================================

a = [172,12,36,76] # Insert here results of model.run(). They are index values of the sensor nodes
Nodos = []
for i in a:
    Nodos.append(network.SensorNodes[i])

for i in Nodos:           
    plt.figure()
    ax = plt.gca()
    NodesDe = network.VOI(N=[i],L=network.LeakNodes)   
    nn = network.NodesDetected     
    wntr.graphics.plot_network(network.OpenNet(),node_size=10,node_cmap='Pastel2_r',node_attribute=network.SensorNodes,\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax)
    wntr.graphics.plot_network(network.OpenNet(),node_size=80,node_cmap='hsv',node_attribute=nn,\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax)   
    wntr.graphics.plot_network(network.OpenNet(),node_size=35,node_cmap='Wistia',node_attribute=[i],\
                           link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax)
        
  
print(network.VOI(N=Nodos,L=network.LeakNodes),len(network.NodesDetected),\
      network.Pro_Detection)



    
    

    