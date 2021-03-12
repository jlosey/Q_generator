#!/usr/bin/python3

import sys
import random
import math
import numpy as np
from scipy import integrate
from scipy.linalg import expm, logm
from scipy import interpolate
import os
from helperFunctions import *
from datetime import datetime


#derivative functions
def potDerBS(x,k=10): #derivative of BS potential
    return k * (x**2 - 1.0) * x


def potDerABS(x,k=10): #derivative of ABS potential
    return 4 * k/20 * x**3 - 2 * k/7.69 * x + k/20


def potDer2HEMS(x,k=10): #derivative of 2HEMS potential
        return 2 * k/6.99 * x - 4 * k/31.5 * x**3 + 6 * k/526 * x**5


def potDerA2HEMS(x,k=10): #derivative of A2HEMS potential
        return 2 * k/6.99 * x - 4 * k/31.5 * x**3 + 6 * k/526 * x**5 + 0.2




def produceTrajectory(numTrajectories, width, potentialFunction, filename, start=0, dt=timestep, T=1, k=10, R=0.001, xc=0, K=0, k_=0, damping=1, beta=beta): #potential function - write to filename
#vectorized production of trajectories as per langevin dynamics.
#Euler Maruyama method
#Arguments:#
#numTrajectories - number of trajectories desired
#width - this is a holder, helps separate between the potentials with width 2 and width 4 for generating random starting areas
#potentialFunction - this is the potDerBS or potDerABS above.  Should be a derivative function
#filename - output filename for data generated
#start - start time
#dt - timestep
#T - end time
#k - force constant for potential
#R - 
#xc - 
#K - 
#k_ -
#damping - 1 is overdamped langevin
#beta - thermodynamic beta
    if filename != "test":
        f = open(filename, "a")
    #vectorized
    g = math.sqrt(2/(beta*dt*damping))
    t=0
    print(dt)
    _R = np.random.RandomState(numTrajectories)
    #randomSampling = _R.random_sample(numTrajectories)#random numbers from 0-1
    #randomSampling = np.linspace(0.0, 1.0, num=numTrajectories) #Not so random, testing a theory
    randomSampling = float(start)
    if width == 2:
        X = (randomSampling - 0.5) * 4
    elif width == 4:
        X = (randomSampling - 0.5 ) * 8
    else:
        print("error")

    v = 0
    traj = np.zeros((int(T/dt)+2, numTrajectories))
    #deviateAvailable = False
    step=0
    while (t<T):
        #randomd, deviateAvailable = randn(0,1,deviateAvailable)
        randomd = np.random.normal(loc=0, scale=1, size=numTrajectories)
        v = g * randomd - potentialFunction(X) / damping
        #v =  g* randomd - (k * (x**2 - 1) * x)
        X += dt * v

        if (R>0 and step<len(traj)):# and (t/dt)%(R/dt)==0):
            #bound X  from 0 to numberBins-1.
            value =  np.round(((X+width) * (numberBins-1) / 2 / width ))
            if value > numberBins-1 or value < 0:
                traj[step]=-1
            #####
            else:
                traj[step] = value % numberBins
            step = step + 1
        t += dt
        if(t%(dt*10000)<dt):
            print("%s ::: Currently: %f"%(datetime.now().strftime("%H:%M:%S"),int(t/dt)))
    if filename != "test":
        for i in traj:
            f.write(str(i))
            f.write("\n")
        f.close()
    return np.transpose(traj.astype(int))




def binTraj(trajectories, filePrefix='test', width=2):
#Bin trajectory(ies) into NxN transition matrix
#Arguments:#
#trajectories - Long Np Arrays for counting
#filePrefix - fileprefix for output
#width - same as above, used to differentiate math between width 2 and width 4 potentials
#returns:#
#outputs files for each of the bins counted.
#Todo: Vectorize this step (definitely slowest step right now)
#Todo: make the matrices something that can be put in as input variable.
#Div for paper [1,10,50,100,500,1000]
#Tau for paper[1,5,10,50,100,500,1000]
    matrix1=np.zeros((numberBins,numberBins))
    matrix5=np.zeros((numberBins,numberBins))
    matrix10=np.zeros((numberBins,numberBins))
    matrix50=np.zeros((numberBins,numberBins))
    matrix100=np.zeros((numberBins,numberBins))
    matrix500=np.zeros((numberBins,numberBins))
    matrix1000=np.zeros((numberBins,numberBins))

    counter = 0
    print("Trajectories Shape:")
    print(trajectories.shape)
    for trajectory in trajectories.astype(int):
        j=0
        trajectory = np.trim_zeros(trajectory, 'b')
        counter=0
        for i in range(len(trajectory)):
            if trajectory[i]==-1:
                continue
            if counter > 0:
                if trajectory[i-1]!=-1:
                    matrix1[trajectory[i-1],trajectory[i]]+=1
            if counter > 5:
                if trajectory[i-5]!=-1:
                    matrix5[trajectory[i-5],trajectory[i]]+=1
            if counter > 10:
                if trajectory[i-10]!=-1:
                    matrix10[trajectory[i-10],trajectory[i]]+=1
            if counter > 50:
                if trajectory[i-50]!=-1:
                    matrix50[trajectory[i-50],trajectory[i]]+=1
            if counter > 100:
                if trajectory[i-100]!=-1:
                    matrix100[trajectory[i-100],trajectory[i]]+=1
            if counter > 500:
                if trajectory[i-500]!=-1:
                    matrix500[trajectory[i-500],trajectory[i]]+=1
            if counter > 1000:
                if trajectory[i-1000]!=-1:
                    matrix1000[trajectory[i-1000],trajectory[i]]+=1

#Div for paper [1,10,50,100,500,1000]
#Tau for paper[1,5,10,50,100,500,1000]
            counter+=1
            if j==0 and counter==int(len(trajectory)/1000):
                print("Matrices written @ %s"%counter)
                writeMatrix(matrix1, "%s_%s_1.mat"%(filePrefix,1000))
                writeMatrix(matrix5, "%s_%s_5.mat"%(filePrefix,1000))
                writeMatrix(matrix10, "%s_%s_10.mat"%(filePrefix,1000))
                writeMatrix(matrix50, "%s_%s_50.mat"%(filePrefix,1000))
                writeMatrix(matrix100, "%s_%s_100.mat"%(filePrefix,1000))
                writeMatrix(matrix500, "%s_%s_500.mat"%(filePrefix,1000))
                writeMatrix(matrix1000, "%s_%s_1000.mat"%(filePrefix,1000))
                j += 1
            if j==1 and counter==int(len(trajectory)/500):
                print("Matrices written @ %s"%counter)
                writeMatrix(matrix1, "%s_%s_1.mat"%(filePrefix,500))
                writeMatrix(matrix5, "%s_%s_5.mat"%(filePrefix,500))
                writeMatrix(matrix10, "%s_%s_10.mat"%(filePrefix,500))
                writeMatrix(matrix50, "%s_%s_50.mat"%(filePrefix,500))
                writeMatrix(matrix100, "%s_%s_100.mat"%(filePrefix,500))
                writeMatrix(matrix500, "%s_%s_500.mat"%(filePrefix,500))
                writeMatrix(matrix1000, "%s_%s_1000.mat"%(filePrefix,500))
                j += 1
            if j==2  and counter==int(len(trajectory)/100):
                print("Matrices written @ %s"%counter)
                writeMatrix(matrix1, "%s_%s_1.mat"%(filePrefix,100))
                writeMatrix(matrix5, "%s_%s_5.mat"%(filePrefix,100))
                writeMatrix(matrix10, "%s_%s_10.mat"%(filePrefix,100))
                writeMatrix(matrix50, "%s_%s_50.mat"%(filePrefix,100))
                writeMatrix(matrix100, "%s_%s_100.mat"%(filePrefix,100))
                writeMatrix(matrix500, "%s_%s_500.mat"%(filePrefix,100))
                writeMatrix(matrix1000, "%s_%s_1000.mat"%(filePrefix,100))
                j += 1
            if j==3  and counter==int(len(trajectory)/50):
                print("Matrices written @ %s"%counter)
                writeMatrix(matrix1, "%s_%s_1.mat"%(filePrefix,50))
                writeMatrix(matrix5, "%s_%s_5.mat"%(filePrefix,50))
                writeMatrix(matrix10, "%s_%s_10.mat"%(filePrefix,50))
                writeMatrix(matrix50, "%s_%s_50.mat"%(filePrefix,50))
                writeMatrix(matrix100, "%s_%s_100.mat"%(filePrefix,50))
                writeMatrix(matrix500, "%s_%s_500.mat"%(filePrefix,50))
                writeMatrix(matrix1000, "%s_%s_1000.mat"%(filePrefix,50))
                j += 1
            if j==4 and counter==int(len(trajectory)/10):
                print("Matrices written @ %s"%counter)
                writeMatrix(matrix1, "%s_%s_1.mat"%(filePrefix,10))
                writeMatrix(matrix5, "%s_%s_5.mat"%(filePrefix,10))
                writeMatrix(matrix10, "%s_%s_10.mat"%(filePrefix,10))
                writeMatrix(matrix50, "%s_%s_50.mat"%(filePrefix,10))
                writeMatrix(matrix100, "%s_%s_100.mat"%(filePrefix,10))
                writeMatrix(matrix500, "%s_%s_500.mat"%(filePrefix,10))
                writeMatrix(matrix1000, "%s_%s_1000.mat"%(filePrefix,10))
#I->J was used in embedability paper we published.
    print("Final Matrices written @ %s"%counter)
    print(matrix1.sum())
    writeMatrix(matrix1, filePrefix+"_1"+".mat")
    writeMatrix(matrix5, filePrefix+"_5"+".mat")
    writeMatrix(matrix10, filePrefix+"_10"+".mat")
    writeMatrix(matrix50, filePrefix+"_50"+".mat")
    writeMatrix(matrix100, filePrefix+"_100"+".mat")
    writeMatrix(matrix500, filePrefix+"_500"+".mat")
    writeMatrix(matrix1000, filePrefix+"_1000"+".mat")



def printHist(trajectories, width):
#print a histogram for a trajectory
#used for simple validation
#Arguments:#
#trajectories - np array of trajectory data
#width - used for binning as above.
#returns:#
#None
#todo, return the histogram, print from a normal print(), rename to just getHist

    #traj2 = np.zeros(trajectories.shape)
    #for t in range(traj2.shape[0]):
    #    for i in range(traj2.shape[1]):
    #               traj2[t,i] = int(int((trajectories[t,i]+2) * 12 ) % 48)
    (unique, counts) = np.unique(trajectories, return_counts=True)
    for i in counts:
        print(i)





