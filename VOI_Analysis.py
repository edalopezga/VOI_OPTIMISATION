# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 18:23:42 2021

@author: Andres L
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pylab

VOI = []
PriorBelief = np.arange(0.01,1,0.01)
Pro_detection = np.arange(0.01,1,0.01)
FirstD = []
DetectUP = []
NoDetectUP = []


for i in range(len(PriorBelief)):
    VOI_Temp = []
    FirstD_Temp = []
    DetectUP_Temp = []
    NoDetectUP_Temp = []
    for j in range(len(Pro_detection)):
        Cas=[[-5,-15],[-15,0]]
        S = [PriorBelief[i], 1-PriorBelief[i]]
        FalseAlarm = 0.15
        Pro_Detection = Pro_detection[j] 
        Mes_Leak = [Pro_Detection, 1 - Pro_Detection]
        Mes_Noleak = [FalseAlarm,1-FalseAlarm] 
        FlowRequiredD=[1,2,3,4,5,6,7,8,9,10]         
        A = [S[0] * k for k in Mes_Leak]
        B = [S[1] * k for k in Mes_Noleak]
        C = (A[0]+B[0])
        D = (A[1]+B[1])
        E = A[0]/C
        F = A[1]/D
        PSM =[[E,F],[1-E,1-F]]  
               
        G = FlowRequiredD[0] # New leak flow used 
        
        H = [Cas[0][0]-G,Cas[0][1]-G*3] # Recalculation of CAS matrix
        I = [S[0]*l for l in H]  #expected utility decision Go to check
        J = [S[1]*m for m in Cas[1]] #expected utility decision Not to Go t
        U0 = max(I[0]+J[0],I[1]+J[1]) #Best option without extra information
        if I[1]+J[1] > I[0]+J[0]:
            FirstD_Temp.append(0)
        else:
            FirstD_Temp.append(1)            
        #Matrix of new utility with information:
        Cos = [[H[0],Cas[1][0]],[H[1],Cas[1][1]]]
        U1 = np.dot(Cos,PSM)
        U_1 = [max([n[0] for n in U1]),max([n[1] for n in U1])]       
      
         #Value of each message:
        V =[n - U0 for n in U_1] 
        if A[0]==0: #we evaluate if the sensor is avaible to detect a least 
          # one leak event otherwise the value of information is zero
         voi=0
        else:
         voi=(C*V[0]+D*V[1])      
        VOI_Temp.append(voi)
        
        if (U1[0][0] > U1[1][0] and voi>0.1) or I[1]+J[1] < I[0]+J[0] :
            DetectUP_Temp.append(1)
        else:
            DetectUP_Temp.append(0)
        if (U1[0][1] > U1[1][1] and voi>0.1) or I[1]+J[1] > I[0]+J[0] :
            NoDetectUP_Temp.append(1)
        else:
            NoDetectUP_Temp.append(0)         
    FirstD.append(FirstD_Temp)
    VOI.append(VOI_Temp)
    DetectUP.append(DetectUP_Temp)
    NoDetectUP.append(NoDetectUP_Temp)


# Map of VOI gray scale

fig, ax1 = plt.subplots()
ax2 = ax1.twiny()
c = ax1.pcolor(Pro_detection,PriorBelief, VOI, cmap = "binary")
ax2.pcolor(Pro_detection[::-1],PriorBelief, VOI, cmap = "binary")
ax2.invert_xaxis()
ax1.set_xlabel('Probability'+' $q_{m1,s1}$',fontsize=13)
ax2.set_xlabel('Probability'+' $q_{m2,s1}$',fontsize=13)
ax1.set_ylabel('Prior Belief '+' $\Pi_{s}$' +' for $s_{1}$',fontsize=13 )
ax1.set_xticks(np.arange(0, 1.1, step=0.1))
ax2.set_xticks(np.arange(0, 1.1, step=0.1))
ax1.set_yticks(np.arange(0, 1.1, step=0.1))
clb = fig.colorbar(c,ax=ax1)
clb.ax.set_title('VOI')
ax1.tick_params(labelsize=11)
ax2.tick_params(labelsize=11)
plt.axvline(x=0.85,ymin=0,ymax=1,linestyle='--', linewidth=1.2,color='k')
plt.axhline(y=0.52,xmin=0,xmax=1,linestyle='--', linewidth=1.2,color='k')
pylab.text(0.55,0.8, 'I',fontsize=13)
pylab.text(0.65,0.75, 'Intervention')
pylab.text(0.935,0.8, 'II',fontsize=13)
pylab.text(0.95,0.2, 'III', fontsize=13)
pylab.text(0.55,0.2, 'IV', fontsize=13)
pylab.text(0.65,0.15, 'No Intervention')
pylab.text(0.99,0.47, 'Line (x)')
pylab.text(0.83,0.04, 'Line (y)',rotation=90)



# Map of Go / not Go without information
fig,ax = plt.subplots()
ax.pcolor(PriorBelief,Pro_detection, FirstD,cmap = "Pastel1")
ax.set_ylabel('Prior Belief '+' $\Pi_{s1}$' + ' - leak state' )
ax.set_yticks(np.arange(0, 1.1, step=0.1))
ax.axes.xaxis.set_visible(False)
ax2=ax.twinx()
ax2.set_yticks(np.arange(0, 1.1, step=0.1))
ax2.invert_yaxis()
ax2.set_ylabel('Prior Belief '+' $\Pi_{s2}$' + ' - Not leak state' )


   
# Map of Go / not Go if detected
plt.pcolor(PriorBelief,Pro_detection, DetectUP,cmap = "binary")
plt.xlabel('Conditional Probability'+' $q_{m,s1}$')
plt.ylabel('Prior Belief '+' $\Pi_{s1}$' + ' - leak state' )
plt.xticks(np.arange(0, 1.1, step=0.1))
plt.yticks(np.arange(0, 1.1, step=0.1))
plt.show()

   
    
# Map of Go / not Go if Not detected
plt.pcolor(PriorBelief,Pro_detection, NoDetectUP,cmap = "binary_r")
plt.ylabel('Detection Probability')
plt.xlabel('Prior Belief '+' $\Pi_s$' )
plt.xticks(np.arange(0, 1.1, step=0.1))
plt.yticks(np.arange(0, 1.1, step=0.1))
plt.show()  

    
    
    
    
