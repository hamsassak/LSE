import math
import numpy as np
import csv
import scipy.linalg
import time
#import mysql.connector
from subprocess import call
from sys import argv
import requests
import math

importTime = time.time()


'''
Creates a hashtable from the constraint file. This allows us to perform constant time look-ups
'''
#init
impedances = [2.99+1j*2.99,2.15+1j*2.22,2.99+1j*2.99,1.8+1j*2.8,2.2+1j*2.2,2.0+1j*2.0,2.35+1j*2.39,1.75+1j*1.75,1.75+1j*1.75]
admittance = np.divide(1,impedances)
tofrom = np.array([[1,4],[1,3],[2,4],[2,3],[4,7],[7,6],[3,6],[3,5],[4,5]])
cov_v = 10**-5
cov_I = 10**-3
nv = 7
nI = 9
#benchmarking variables
sumT=0
count=0
maxT=0
#W matrix----------------------------------------------------------------
Cov_v = np.diag(np.full(nv,cov_v))
Cov_I = np.diag(np.full(nI,cov_I))
W = scipy.linalg.block_diag(Cov_v,Cov_v,Cov_I,Cov_I)
#number of measurements
voltages = np.zeros((7,2))
currents = np.zeros((9,2))
H0 = np.zeros((2,2*nv))

for i in range(nv):
    Hj = np.array(H0)
    Hj[0][2*i] = 1
    Hj[1][2*i+1] = 1
    if(i == 0):
        Hv = np.array(Hj)
    else:
        Hv = np.append(Hv,Hj,axis=0)
    H0 = np.zeros((2,2*nv))
    #vrect[tofrom[i][0]*admittance[i]][1]
    for i in range(nI):
        #takes note of where the currents are coming from and going to
        #considering attaching that to the current matrix, like how it is in the MATLAB
        To = int(tofrom[i][0]-1)
        From = int(tofrom[i][1]-1)
        #real and imaginary parts
        gj = admittance[i].real
        bj = admittance[i].imag
        #reset Hj matrix
        Hj = np.array(H0)

        Hj[0][To*2] = gj 
        Hj[0][To*2+1] = -bj
        Hj[0][From*2] = -gj 
        Hj[0][From*2+1] = bj
        Hj[1][To*2] = bj
        Hj[1][To*2+1] = gj
        Hj[1][From*2] = -bj
        Hj[1][From*2+1] = -gj
        if(i==0):
            HI = np.array(Hj)
        else:
            HI = np.append(HI,Hj,axis=0)
#a function that takes the magnitude and phase of a measurement
#and returns it in cartesian coordinates
def PtoR(radius, angle):
    real = radius*math.cos(math.radians(angle))
    imaginary = radius*math.sin(math.radians(angle))
    full = real + 1j*imaginary
    return full

HZ = 60 #window size
def buildTables(filename):
		constraints = dict()
		allQueues = dict()
		f = open(filename, 'r')
		for line in f:
				sigID, sel, myType, myMin, myMax, limit = line.strip().split(',')
				constraints[sigID] = (sel,myType,myMin,myMax,limit)
				allQueues[sigID] = [[0]*HZ, [""]*HZ, 0, 0, myType] #measurements, timestamps, index, count
		f.close()
		return constraints, allQueues

constraints, allQueues = buildTables('newConstraints_v2.txt')
startTime = time.time()


verbose = False # print stats, etc. Turn off for automation.
timing = False #print timing stats
innerTime = False #print in loop timing stats
debug = False #print the counter stats, must be turned off for actual processing
def LinearSE(v,I):
    #conversion to rectangular
    for i in range(np.size(v,0)):
        volt = PtoR(v[i][0],v[i][1])
        if (i==0):
            vrect = np.array([[i,volt,volt.real,volt.imag]])
        else:
            vrect = np.append(vrect,[[i,volt,volt.real,volt.imag]],axis=0)
    for i in range(np.size(I,0)):
        cur = PtoR(I[i][0],I[i][1])
        if (i==0):
            Irect = np.array([[i,cur,cur.real,cur.imag]])
        else:
            Irect = np.append(Irect,[[i,cur,cur.real,cur.imag]],axis=0)

    #z matrix for voltages ---------------------------------------------------
    z0 = np.zeros((2,1))
    for i in range(nv):
	    zj = np.array(z0)
	    zj[0] = vrect[i][2]
	    zj[1] = vrect[i][3]
	    if(i == 0):
	    	zv = np.array(zj)
	    else:
	    	zv = np.append(zv,zj,axis=0)
    #z matrix for currents ----------------------------------------------------
    z0 = np.zeros((2,1))
    #vrect[tofrom[i][0]*admittance[i]][1]
    for i in range(nI):
        #takes note of where the currents are coming from and going to
        #considering attaching that to the current matrix, like how it is in the MATLAB
        zj = np.array(z0)
        zj[0] = Irect[i][2]
        zj[1] = Irect[i][3]
        if(i==0):
            zI = np.array(zj)
        else:
            zI = np.append(zI,zj,axis=0)

	#combine both measurements----------------------------------------------
    H = np.append(Hv,HI,axis=0)
    z = np.append(zv,zI,axis=0)
    #I need to go back here and implement a way for my code to connect the currents
    #to their admittances, because as it is right now, the code just assumes they are in order

	#G matrix, G = H' * W * H
    G = H.T @ W @ H
    iG = np.linalg.inv(G)
	#here is the estimate    
    x = np.linalg.inv(H.T@W@H)@H.T @ W @ z
    S_d = np.power(np.diag(iG),0.5).reshape(2*nv,1)
    J = (H @ x - z).T@ W @ ( H @ x - z)
    Jz = np.zeros((32,1))
    for a in range(2*nv):
        Jz[a]=(z[a]-H[a]@x)**2*cov_v
    for a in range(2*nI):
        Jz[a+2*nv] = (z[a+nv*2]-H[a+2*nv]@x)**2*cov_I 
    #sum of residuals
    sumJ = 0
    for i in range(np.size(Jz,0)):
        sumJ = sumJ + Jz[i]
    return(x,S_d,J,sumJ)
while True:
    
    detectionTime = time.time()    
    LSEStartTime = time.time()
    f = open('Model_Validation.csv', 'r')
    w = open('sendToServer', 'w')
    total_constraint_time = 0
    constraint_time_counts = 0
    total_temporal_time = 0
    temporal_time_counts = 0
    stoptime = 0.0
    for line in f:
        micc, sigID, timestamp, value = line.strip().split(',')
        value = float(value)
        if (timestamp != stoptime):
            lse,S_d,J,sumJ = LinearSE(voltages,currents)
            LSEStopTime = time.time()
            t = LSEStopTime-LSEStartTime
            LSEStartTime = time.time()
            sumT +=t
            count+=1
            avgt=sumT/count
            maxT = max(t,maxT)
            print("Here is the LSE: \n This is x: \n {0} \n This is S_d: \n {1} \n This is J: \n {2}\n This is Jsum: \n {3} \n That run took {4} seconds.\n The average time was {5}.\n The max time was {6}.\n".format(str(lse),str(S_d),str(J),str(sumJ),str(t),str(avgt),str(maxT)))
            #print(str(voltages))
            #print(str(currents))
        stoptime = timestamp
        if (constraints[sigID][1] == 'voltage') and int(constraints[sigID][0])<8:
            voltages[int(constraints[sigID][0])-1][0] = value
        if (constraints[sigID][1] == 'voltage phase angle') and int(constraints[sigID][0])<8:
            voltages[int(constraints[sigID][0])-1][1] = value
        if (constraints[sigID][1] == 'current') and int(constraints[sigID][0])<10:
            if (int(constraints[sigID][0])==5 or int(constraints[sigID][0])==7):
                currents[int(constraints[sigID][0])-1][0] = -1*value
            else:
                currents[int(constraints[sigID][0])-1][0] = value            
        if (constraints[sigID][1] == 'current phase angle') and int(constraints[sigID][0])<10:
            currents[int(constraints[sigID][0])-1][1] = value
    f.close()
    w.close()

    data = {}
    with open("sendToServer","r") as dataFile:
        count = 0
        for line in dataFile:
            data[count] = line
            count += 1

    if not data:
        data = {"Anoms":"No"}
    else:
        data["Anoms"] = "Yes"
        
    #r = requests.post("http://192.168.20.200:1234/alert",data=data)
    #print(r.text)
        
print("done")