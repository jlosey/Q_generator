#!/usr/bin/python3


import numpy as np
import os
from helperFunctions import *
import glob

#Div for paper [1,10,50,100,500,1000]
#Tau for paper[1,5,10,50,100,500,1000]

#file convention is mats/prefix_divisor_tau.mat

taus = ['1', '5', '10', '50', '100', '500',  '1000']
divs = ['10', '50', '100', '500', '1000']


#TODO: generalize
bins = 48
fileGlob = "mats/*BS_"

for div in divs:
    for tau in taus:
        totalMat = np.zeros((bins,bins))
        print(totalMat.shape)
        sumd=0
        f=0
        for i in glob.glob(fileGlob+div+"_"+tau+".mat"):
            addMat = readMatrix(i)
            totalMat += addMat
            sumd+=addMat.sum()
            f=f+1
        writeMatrix(totalMat, div+"tau_"+tau+".mat")
        sumd= sumd/50000000000.0
        print("Matrix: %stau_%s.mat\ndivisor: %s\ntau: %s\ntotalPoints: %s\nfiles: %s\n\n"%(div,tau,div,tau,sumd,f))


#get the div=1 outputs
for tau in taus:
        totalMat = np.zeros((bins,bins))
        sumd=0
        f=0
        for i in glob.glob(fileGlob+tau+".mat"):
            addMat = readMatrix(i)
            totalMat += addMat
            sumd+=addMat.sum()
            f=f+1
        writeMatrix(totalMat, "1tau_"+tau+".mat")
        sumd= sumd/50000000000.0
        print("Matrix: 1tau_%s.mat\ndivisor: 1\ntau: %s\ntotalPoints: %s\nfiles: %s\n\n"%(tau,tau,sumd,f))




