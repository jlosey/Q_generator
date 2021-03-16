#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Load modules
import glob
from zipfile import ZipFile
import io as iio
import numpy as np
import scipy as sp
#Load local files needed for calulations
from custom_constants import *
from helperFunctions import *
from GeneratorFunctions import *


# After compiling the trajectories in transition matricies located in the Data folder, the Q matrix is generated in the cell below.

# In[13]:


zip_name = "Compiled/bootstrap_data.zip"
method_list = ["DA","WA","QOG","CWO","QP","PA"]
function_list = [diagonalAdjustment,weightedAdjustment,QOG,CWO,qp_rate,polySolver]

# opening the zip file in READ mode 
with ZipFile(zip_name, 'r') as z: 
    for file in z.namelist():        
        fsplit1 = file.split('/')[-1]
        fsplit2 = fsplit1.split('_')
        div = fsplit2[0].split('tau')[0]
        tau = fsplit2[1].split('.')[0]
        seed = fsplit2[2]
        print(div,tau,seed)
        #loop over estimators
        for fi,mi in zip(function_list,method_list):
            #open each file in zip
            with z.open(file) as myfile:
                mat_string = myfile.read().decode('UTF-8')
                arr = np.fromstring(mat_string,dtype=np.float32,sep=' ')
                mat = arr.reshape((47,47))
                #run each estimator on the matrix            
                Q = fi(mat,int(tau))
                #write matrix to boostrap folder
                writeMatrix(Q, "Data/bootstrap/"+div+"tau"+f"_{tau}"+".mat"+f"_{mi}"+f"_{seed}")
                #input()
                


# In[5]:


method_list = ["DA","WA","QOG","CWO","QP","PA"]
function_list = [diagonalAdjustment,weightedAdjustment,QOG,CWO,qp_rate,polySolver]

for i in glob.glob("Compiled/*.mat"):
    file_name = i.split("/")[-1]
    tau = int(file_name.split("tau_")[-1].split(".")[0])
    mat=readMatrix(i)
    for fi,mi in zip(function_list,method_list):
        Q = fi(mat,tau)
        writeMatrix(Q, "Data2/"+file_name+f"_{mi}")


# In[7]:


SWITCH=0
if SWITCH==1:
    for i in glob.glob("Compiled/*.mat"):
        a=readMatrix(i)
        b=diagonalAdjustment(a)
        writeMatrix(b, i+"_DA")
        a=readMatrix(i)
        b=weightedAdjustment(a)
        writeMatrix(b, i+"_WA")
        a=readMatrix(i)
        b=QOG(a)
        writeMatrix(b, i+"_QOG")
        a=readMatrix(i)
        b=CWO(a)
        writeMatrix(b, i+"_CWO")
    #EM and MLE are done on supercomputers/BSLGPUs


# EM and MLE were generated using HPC resources as these algorithms are computationally expensive.
# Running EM and MLE was done using the 

# In[ ]:





# In[ ]:


#Calc MFPT


# In[ ]:


#Plot Potential

