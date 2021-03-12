#!/usr/bin/python


import sys
import random
import math
import numpy as np
from scipy import integrate
from scipy.linalg import expm, logm
from scipy import interpolate
import os

from helperFunctions import *
from trajectories import *

from constants import *

numRun=float(sys.argv[1])
print("Num Run: ", numRun)
start=numRun/48.0 #Even spacing between 0,1
print("Start:" ,start)

#production BS
#Purposefully Throw an Error - This has been run.
#traj1 = (produceTrajectory(5, width, potDerBS,"BS.traj", dt=timestep, T=.0001))


######bintraj only works with numtrajectories == 1!!! #######
traj1 = (produceTrajectory(1, width, potDerBS,"test", dt=timestep, T=TotalTime, start=start))
##############################################################
#Div for paper [1,10,50,100,500,1000]
#Tau for paper[1,5,10,50,100,500,1000]
print("traj made")
#printHist(traj1, width)
#print("hist made")
binTraj(traj1, "%dBS"%numRun,width)
print("traj Binned")
