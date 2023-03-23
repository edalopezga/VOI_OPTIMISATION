# =============================================================================
# Packages Required
# =============================================================================
'''
Most of the libraries needed to run this utility are installed with python itself.

Some of the packages not included and to be installed are:
    
Platypus: Platypus is a framework for evolutionary computing in Python with a 
focus on multiobjective evolutionary algorithms (MOEAs).
more information: https://platypus.readthedocs.io/en/latest/index.html 

WNTR: The Water Network Tool for Resilience (WNTR, pronounced winter) is a Python 
package designed to simulate and analyze resilience of water distribution 
networks. More information in:
https://wntr.readthedocs.io/en/latest/overview.html

Plotly: Plotly's Python graphing library makes interactive, publication-quality 
graphs. More information in: https://plotly.com/python/
'''  
  
# =============================================================================
# Libraries Required.
# =============================================================================
import wntr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from platypus import NSGAII, Problem, Subset, unique, nondominated
import time
from collections import Counter

import plotly.express as px
from plotly.offline import plot
import random


# =============================================================================
# Fuctions Required:
# =============================================================================
'''
To save space in the code and reuse the instructions in multiple scenarios, 
an specific class called Network and different functions were developed.

In addition, some of the functions incorporated in the different packages have 
a difficult structure to memorize and read, for this reason they were rewritten
'''

def RunNet (x):
    '''
    This function receives as variable an instance of a network created with the 
    OpenNet() function.
    This returns an object with all the results of the model, the queries have to 
    be related to the instance created with this function. 
    How to do the querries, visit: 
    https://wntr.readthedocs.io/en/latest/waternetworkmodel.html
    
    Example: Results = RunNet(MyInstance)'''        
    sim = wntr.sim.EpanetSimulator(x)
    sim = sim.run_sim(version=2.2)
    return sim

class Network:
    def __init__(self, file):
        self.file = file
        
    def OpenNet(self):
        '''
        This function receives the name of the file in string format with the extention
        .inp, it creates an instance of the network, it is important to open the folder
        where the file is located. 
        
        Example:  MyInstance =  OpenNet('Mynetwork.inp')
        
        With this instance we can modify physical charateristics of the network an
        options of modelation before Running the model.        
        
        Note: By now is important that the network created with EPANET should be 
        saved using version 2.2, with SI units, 24 hours of simulation, 1 hour interval,
        demand model = PDD, these parameters should be modified in EPANET for EASIER 
        use of the tool.
        '''
        return wntr.network.WaterNetworkModel(self.file)     

    def GraphNet(self):
        '''
        This function returns a graph of the network.      
        Example:  class.GraphNet()'''             
        wntr.graphics.plot_network(self.OpenNet(), title= 'Network')
   
    def BaseDemand (self, y = "No_defined"):
        '''
        This function rturns the base demand of the nodes y (in list format) to 
        be consulted. If any node is specified then all the base-demands are showed.       
        '''  
        x = self.OpenNet()
        BD = []
        if y == "No_defined":
            for i in x.junction_name_list:        
                Node = x.get_node(i)
                if Node.base_demand>0:            
                    BD.append(Node.base_demand)  
        else:  
            for i in y:
                Node = x.get_node(i)
                if Node.base_demand>0: 
                        BD.append(Node.base_demand) 
        return BD
    
# =============================================================================
# MAIN FUCTION
# =============================================================================

    def CalculateParameters(self,Leak_Nodes='All',Sensor_Nodes='All',Time=1, Method ='PDD',\
          Leakmin='None',Leakmax='None',Num_Leakages='None',Threshold=0.5): 

        start = time.time()  #start of the Timer, just for evaluate code performance.    
        Net = self.OpenNet()  #Instance of the network    
        if Leak_Nodes == 'All':
            Leak_Nodes = Net.junction_name_list #Nodos where leaks are going to be placed
    
        if Sensor_Nodes == 'All':
            Sensor_Nodes = Net.junction_name_list     
  
        Net.options.hydraulic.demand_model=Method   #Using PDD as demand method
        Net.options.time.duration = Time*24*3600    #Time of simulation
        Net_Sim = RunNet(Net)  #Run network, base case.
        St = int(Net.options.time.duration/Net.options.time.hydraulic_timestep)    
        BPM = Net_Sim.node['pressure'].loc[1:St*3600,Net.junction_name_list] #Base pressures
        BD = sum(self.BaseDemand())/len(self.BaseDemand()) #Average Demand of the system in m3/s

        if Leakmin == 'None':
            Leakmin = 0.1*BD*1000  #Minimum Leak added to nodes L/s
        if Leakmax == 'None':
            Leakmax = BD * 1000 #Maximum Leak added to nodes L/s
            if Leakmax <= Leakmin:
                Leakmax = 2*Leakmin
        if Num_Leakages == 'None': 
            Num_Leakages = 9
            Leakinterval = (Leakmax - Leakmin)/Num_Leakages
        if Num_Leakages != 'None':
            Leakinterval = (Leakmax - Leakmin )/(Num_Leakages)          
        LeakFlows = np.arange(Leakmin,Leakmax+0.001,Leakinterval)
    
        LPM = []        #Leak pressure matrix        
        LM = []         #Leakage Matrix
        DM = []         #Divergence matrix
    
        for i in Leak_Nodes:       
            for k in range(len(LeakFlows)):  
                    
              LeakFlow = [LeakFlows[k]/1000]*(24*Time+1) #array of the leak flow (m3/s)         
              Net.add_pattern(name ='New', pattern = LeakFlow) #Add New Patter To the model
              Net.get_node(i).add_demand(base= 1, pattern_name='New') #Add leakflow
              Net.options.time.duration = 24*Time*3600 #Time of simulation
              Net_New = RunNet(Net) # Run new model                        
              Net2 = Net_New.node['pressure'].loc[1:St*3600,Sensor_Nodes].\
                 rename_axis('Node_'+i+', '+str(round(LeakFlows[k],2))+'LPS',axis=1) # Give name to the dataframe
              LPM.append(Net2) # Save pressure results
              Difference = BPM[Sensor_Nodes].sub(Net2, fill_value=0) # Create divergence matrix 
              DM.append(Difference.abs().rename_axis('Node_'+ i+', '\
                    +str(round(LeakFlows[k],2))+'LPS', axis=1)) # Save Divergence M.                             
              LM.append(pd.DataFrame([k*1000 for k in LeakFlow[1:]], columns = ['LeakFlow']\
                            , index =list(range(3600,St*3600+3600,3600))).\
                            rename_axis('Node: ' + i, axis=1)) #Save leakflows used
              
              # Restore initial condictions:              
              Net = self.OpenNet()
              Net.options.hydraulic.demand_model = Method # Change type of simulation to PDD
              Net.options.time.duration = 24*Time*3600
              Net_Sim = RunNet(Net)      
        
        TM_ = []   # time when the leak was identified
        WLM = []   # Water loss Matrix (L/s), how much water is wasted
    
        for i in range(len(DM)):    
            TMtemp = [] 
            WLMtemp = []
            for j in Sensor_Nodes:
                WLMtemp2 = [] 
                for k in range(len(DM[0])):
                    if DM[i][j][(k+1)*3600] <= Threshold:                
                        WLMtemp2.append(LM[i].LeakFlow[(k+1)*3600]*3600)              
                    else:                
                        WLMtemp2.append(LM[i].LeakFlow[(k+1)*3600]*3600)
                        break
                TMtemp.append(k+1)
                WLMtemp.append(sum(WLMtemp2))        
            TM_.append(TMtemp)
            WLM.append(WLMtemp) 
        '''
    These matrices indicate the time in which each sensor could detect a leakage
    here we created a dataframe indicating where the perturbation
    was placed and the time of detection of each sensor.
        '''         
        TM = []
        for i in np.arange(0,len(DM),len(LeakFlows)):
            TM.append(pd.DataFrame(TM_[i:i+len(LeakFlows)],columns=Sensor_Nodes,\
                     index=LeakFlows). rename_axis('Node'+ Leak_Nodes[int(i/len(LeakFlows))]))     

        FP = []    #Flow required for causing a detectable perturbation.   
        
        for k in range(len(TM)):
            FPtemp=[]
            for j in Sensor_Nodes:       
                for i in range(len(LeakFlows)):
                    if TM[k][j][LeakFlows[i]] == Time*24:
                        pass
                    else:
                        FPtemp.append(LeakFlows[i])
                        break
                if TM[k][j][LeakFlows[i]] == Time*24:
                    FPtemp.append(1000000)
            FP.append(FPtemp)        
        FP = pd.DataFrame(FP, columns=Sensor_Nodes, index=Leak_Nodes)
        
        '''
    In this part the probability of detection of each sensor is calculate, then
    we counted how many times in the multiple scenarios a sensor could detect the
    leakage, this gave us a matrix of probability of detection and not detection
        '''
        Mes_Leak = []
        
        for i in Sensor_Nodes:
            Not_Detect = []
            for j in range(len(Leak_Nodes)):
                a = TM[j][i].value_counts().to_dict()
                if (Time*24 in a) == True:
                    Not_Detect.append(a[Time*24])
                else:
                    Not_Detect.append(0)
            Pro_D = 1-sum(Not_Detect)/(len(Leak_Nodes)*len(LeakFlows))
            Mes_Leak.append([Pro_D,(1-Pro_D)])
        
        Mes_Leak = pd.DataFrame(Mes_Leak,columns=['Detection','No Detection'],index=Sensor_Nodes)    
       
        self.LeakNodes = Leak_Nodes
        self.SensorNodes = Sensor_Nodes
        self.Leakmin = Leakmin
        self.Leakmax = Leakmax
        self.Leakinterval = Leakinterval
        self.LeakFlow = LeakFlows
        self.TimeDetection = TM
        self.FlowRequired = FP
        self.Perc_Detection = Mes_Leak        
        self.LeakPressureMatrix = LPM
        self.BasePressureMatrix = BPM
        self.DivergenceMatrix = DM
        self.NumLeakages = Num_Leakages
        
        
        print ("it took", time.time() - start, "seconds.") #Time taken by the algorithm
# =============================================================================
# Objective Functions 
# =============================================================================

# Value Of Information (VOI)

    def VOI(self,FalseAlarm=0.15,Cas=[[-5,-15],[-15,0]],\
            N='List of Sensor Nodes',L= 'List of Leak Nodes', \
            PriorityNodes ='EmptyList',VOIMultiplier=3,PriorBelief=0.5):
      
      NodesDetected = []
      QNodesDetected = [] 
      FlowRequiredD  = []  
      for i in L:        
          if (self.FlowRequired[N].loc[i] <= self.LeakFlow[-1]).any():
              NodesDetected.append(i)
              Lf = self.FlowRequired[N].loc[i].min()                        
              a = self.LeakFlow.tolist().index(Lf)             
              QNodesDetected.append(len(self.LeakFlow)-a) 
              FlowRequiredD.append(a+1)
          else:
              FlowRequiredD.append(self.FlowRequired[N].loc[i].min())                   
      Pro_Detection = sum(QNodesDetected)/(len(self.LeakNodes)*len(self.LeakFlow))     
      Mes_Leak = [Pro_Detection, 1 - Pro_Detection]
      Mes_Noleak = [FalseAlarm,1-FalseAlarm] 
      
      VOI = []       
      for i in range(len(L)):     
          S = [PriorBelief,1-PriorBelief]            
          A = [S[0] * k for k in Mes_Leak]
          B = [S[1] * k for k in Mes_Noleak]
          C = (A[0]+B[0])
          D = (A[1]+B[1])
          E = A[0]/C
          F = A[1]/D
          PSM =[[E,F],[1-E,1-F]]                 
          G = FlowRequiredD[i] # New leak flow used       
          H = [Cas[0][0]-G,Cas[0][1]-G*3] # Recalculation of CAS matrix
          I = [S[0]*l for l in H]  #expected utility decision Go to check
          J = [S[1]*m for m in Cas[1]] #expected utility decision Not to Go t
          U0 = max(I[0]+J[0],I[1]+J[1]) #Best option without extra information
        #Matrix of new utility with information:
          Cos = [[H[0],Cas[1][0]],[H[1],Cas[1][1]]]
          U1 = np.dot(Cos,PSM)
          U_1 = [max([n[0] for n in U1]),max([n[1] for n in U1])]
         #Value of each message:
          V =[n - U0 for n in U_1] 
          if A[0]==0: #we evaluate if the sensor is avaible to detect a least 
          # one leak event otherwise the value of information is zero
           VOI.append(0)
          else:
            if L[i] in PriorityNodes:
              VOI.append((C*V[0]+D*V[1])*VOIMultiplier)
            else:   
              VOI.append(C*V[0]+D*V[1])
      self.NodesDetected = NodesDetected
      self.Pro_Detection = Pro_Detection       
      return(sum(VOI))   

    def DelNodes(self,List):
        Netw = self.OpenNet()
        Nodes = Netw.junction_name_list
        return([x for x in Nodes if x not in List]) 

        
    def GraphDetection(self,LeakNode,LeakNumber ='All'):
        if LeakNumber =='All':
            LeakN = list(np.arange(len(self.LeakFlow)))
        else:
            LeakN = LeakNumber       
        LeakNodeIndex = self.LeakNodes.index(LeakNode)
        for i in LeakN:
            Leak = self.LeakFlow[i]
            TimeDetection = self.TimeDetection[LeakNodeIndex].loc[Leak]
            AffectedNodes = TimeDetection[TimeDetection<24]
            plt.figure()
            ax = plt.gca()
            wntr.graphics.plot_network(self.OpenNet(),node_size=10,node_cmap='Pastel2_r',node_attribute=self.SensorNodes,\
                                   link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax,add_colorbar=False,title='Leak '+str(Leak)+' L/s')
            wntr.graphics.plot_network(self.OpenNet(),node_size=20,node_cmap='gray',node_attribute=AffectedNodes,\
                                   link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax,add_colorbar=False)
            wntr.graphics.plot_network(self.OpenNet(),node_size=80,node_cmap='hsv',node_attribute=[LeakNode],\
                                   link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax,add_colorbar=False)
    


    def GraphPump(self, Leak, Node): 
        Time = 1
        Pumps = ['Apollo', 'GA-Braila', 'RaduNegru2', 'RaduNegruMare']    
        Net = self.OpenNet()
        Nresults = RunNet(Net)    
        Net.options.hydraulic.demand_model='PDD'   #Using PDD as demand method    
        LeakFlow = [Leak/1000]*(24*Time+1) #array of the leak flow (m3/s)         
        Net.add_pattern(name ='New', pattern = LeakFlow) #Add New Patter To the model
        Net.get_node(Node).add_demand(base= 1, pattern_name='New') #Add leakflow
        Net.options.time.duration = 24*Time*3600 #Time of simulation
        Net_New = RunNet(Net) # Run new model
        
        PumpDemandBase = (Nresults.node['demand'].loc[:,Pumps]*-1000).rename_axis\
               ('Node: '+Node+', '+str(Leak)+' L/s',axis=1)
                         
        PumpDemands = (Net_New.node['demand'].loc[:,Pumps]*-1000).rename_axis\
               ('Node: '+Node+', '+str(Leak)+' L/s',axis=1) # Give name to the dataframe
        time = np.arange(0,25,Time)
        
        fig = plt.figure(figsize =(6,8))  
        
        ax1 = plt.subplot2grid((7,2),(4,0))     
        plt.plot(time,PumpDemands[Pumps[0]],color='b',linewidth=1.2)
        plt.plot(time,PumpDemandBase[Pumps[0]], color = 'k', linestyle = '--',linewidth=0.5)    
        Error = 100*(PumpDemands[Pumps[0]].sum()-PumpDemandBase[Pumps[0]].sum())/PumpDemandBase[Pumps[0]].sum()
        ax1.title.set_text(Pumps[0] + ' -Dif: ' + "{:.2f}".format(Error)+'%')
        ax1.set_ylabel('Flow L/s' ,fontsize=10 )
        ax1.set_xlabel('Time - Hours' ,fontsize=10 )
        plt.grid(linestyle ='--', linewidth=0.5, color='gray',alpha=0.4)
    
        ax2 = plt.subplot2grid((7,2),(4,1)) 
        plt.plot(time,PumpDemands[Pumps[1]],color='r',linewidth=1.2)
        plt.plot(time,PumpDemandBase[Pumps[1]],color = 'k', linestyle = '--',linewidth=0.5)
        Error = 100*(PumpDemands[Pumps[1]].sum()-PumpDemandBase[Pumps[1]].sum())/PumpDemandBase[Pumps[1]].sum()
        ax2.title.set_text(Pumps[1] + ' -Dif: ' + "{:.2f}".format(Error)+'%') 
        ax2.set_xlabel('Time - Hours' ,fontsize=10 )
        plt.grid(linestyle ='--', linewidth=0.5, color='gray',alpha=0.4)
        
        ax3 = plt.subplot2grid((7,2),(6,0))
        plt.plot(time,PumpDemands[Pumps[2]],color='g',linewidth=1.2)
        plt.plot(time,PumpDemandBase[Pumps[2]],color = 'k', linestyle = '--',linewidth=0.5) 
        Error = 100*(PumpDemands[Pumps[2]].sum()-PumpDemandBase[Pumps[2]].sum())/PumpDemandBase[Pumps[2]].sum()
        ax3.title.set_text(Pumps[2] + ' -Dif: ' + "{:.2f}".format(Error)+'%')
        ax3.set_ylabel('Flow L/s' ,fontsize=10 )
        ax3.set_xlabel('Time - Hours' ,fontsize=10 )
        plt.grid(linestyle ='--', linewidth=0.5, color='gray',alpha=0.4)
        
        ax4 = plt.subplot2grid((7,2),(6,1))
        plt.plot(time,PumpDemands[Pumps[3]],color='y',linewidth=1.2)
        plt.plot(time,PumpDemandBase[Pumps[3]],color = 'k', linestyle = '--',linewidth=0.5) 
        Error = 100*(PumpDemands[Pumps[3]].sum()-PumpDemandBase[Pumps[1]].sum())/PumpDemandBase[Pumps[3]].sum()
        ax4.title.set_text(Pumps[3] + ' -Dif: ' + "{:.2f}".format(Error)+'%')
        ax4.set_xlabel('Time - Hours' ,fontsize=10 )
        plt.grid(linestyle ='--', linewidth=0.5, color='gray',alpha=0.4)
        
        ax5 = plt.subplot2grid((7,2),(0,0),colspan=2, rowspan=4)  
        
        wntr.graphics.plot_network(Net,node_size=80,node_cmap='bwr',node_attribute=[Pumps[0]],\
                               link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax5,title='Leak '+str(Leak)+' L/s at '+Node )
        wntr.graphics.plot_network(Net,node_size=80,node_cmap='bwr_r',node_attribute=[Pumps[1]],\
                               link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax5,)
        wntr.graphics.plot_network(Net,node_size=80,node_cmap='Set3_r',node_attribute=[Pumps[2]],\
                               link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax5)
        wntr.graphics.plot_network(Net,node_size=80,node_cmap='Accent',node_attribute=[Pumps[3]],\
                               link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax5)        
        wntr.graphics.plot_network(Net,node_size=80,node_cmap='cool',node_attribute=[Node],\
                               link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax5)        
        wntr.graphics.plot_network(Net,node_size=40,node_cmap='cool_r',node_attribute=[Node],\
                               link_cmap='Paired',link_width=0.5,link_alpha=0.1,ax=ax5)



def Rank(A):
    B = np.reshape(A,(1,np.shape(A)[0]*np.shape(A)[1]))
    C = Counter(B[0]).items()
    D = pd.DataFrame.from_dict(C).sort_values(by=[1],ascending=(False)\
                                              ,ignore_index=(True))
    return (D.rename(columns={0: 'Node', 1:'Selected'}))  

def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        words = fileObj.read().splitlines() #puts the file into an array
        fileObj.close()
        return words








    





