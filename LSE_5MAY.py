import math
import numpy as np
import csv
import scipy.linalg
import numpy.linalg as lin
import time
import mysql
import mysql.connector


'''
Creates a hashtable from the constraint file. This allows us to perform constant time look-ups
'''
impedances = [[2.99+1j*2.99],[2.15+1j*2.22],[2.99+1j*2.99],[1.8+1j*2.8],[2.2+1j*2.2],[2.0+1j*2.0],[2.35+1j*2.39],[1.75+1j*1.75],[1.75+1j*1.75]]
impedances = np.array(impedances)
cov_v = 10**-5
cov_I = 10**-3

voltages = np.zeros((7,2))
currents = np.zeros((9,2))
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

constraints, allQueues = buildTables('c:/Users/Administrator/Desktop/Original Heuristic AD/newConstraints_v2.txt')
startTime = time.time()
from subprocess import call
#import	mysql
from sys import argv
import requests
import math

importTime = time.time()
verbose = True # print stats, etc. Turn off for automation.
timing = False #print timing stats
innerTime = False #print in loop timing stats
debug = False #print the counter stats, must be turned off for actual processing
x = 0








def LinearSE(y,v,I,cov_v,cov_I):
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
    tofrom = np.array([[1,4],[1,3],[2,4],[2,3],[4,7],[7,6],[3,6],[3,5],[4,5]])
    #making sure arrays given are numpy arrays
    y = np.array(y)
    v = np.array(v)
    I = np.array(I)
    #number of measurements
    nv = np.size(v,0)
    nz = np.size(y,0)
    nI = np.size(I,0)

    Cov_v = cov_v*np.ones((nv,1))
    Cov_I = cov_I*np.ones((nI,1))
    admittance = np.divide(1,y)
    for i in range(np.size(admittance,0)):
        if (i==0):
            Y = np.array([admittance[i].real,admittance[i].imag])
        else:
            Y = np.append(Y,[admittance[i].real,admittance[i].imag],axis=0)

    #W matrix----------------------------------------------------------------
    #at the moment the if - else: is commented out because else never happens
    #if (np.size(Cov_v,0 == 1) or np.size(Cov_v,1)==1):
    Cov_v = np.diag(np.full(nv,cov_v))
    iCov_v = np.diag(np.full(nv,1/cov_v))

    #else:
    #	iCov_v = np.linalg.inv(Cov_v)
    #if (np.size(Cov_I,0)==1 or np.size(Cov_I,1)==1):
    Cov_I = np.diag(np.full(nI,cov_I))
    iCov_I = np.diag(np.full(nI,1/cov_I))
    #else:
    #	iCov_I = np.linalg.inv(Cov_I)
    W = scipy.linalg.block_diag(Cov_v,Cov_v,Cov_I,Cov_I)

    #H matrix for voltages ---------------------------------------------------
    H0 = np.zeros((2,2*nv))
    z0 = np.zeros((2,1))
    for i in range(nv):
	    ji = i + nv
	    Hj = np.array(H0)
	    Hj[0][2*i] = 1
	    Hj[1][2*i+1] = 1
	    zj = np.array(z0)
	    zj[0] = vrect[i][2]
	    zj[1] = vrect[i][3]
	    if(i == 0):
	    	Hv = np.array(Hj)
	    	zv = np.array(zj)
	    else:
	    	Hv = np.append(Hv,Hj,axis=0)
	    	zv = np.append(zv,zj,axis=0)


    #H matrix for currents ----------------------------------------------------
    H0 = np.zeros((2,2*nv))
    z0 = np.zeros((2,1))
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

        zj = np.array(z0)
        zj[0] = Irect[i][2]
        zj[1] = Irect[i][3]
        if(i==0):
            HI = np.array(Hj)
            zI = np.array(zj)
        else:
            HI = np.append(HI,Hj,axis=0)
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
    C = np.linalg.inv(H.T @ W @ H)

    Jz = np.zeros((32,1))
    for a in range(2*nv):
        Jz[a]=(z[a]-H[a]@x)**2*cov_v
    for a in range(2*nI):
        Jz[a+2*nv] = (z[a+nv*2]-H[a+2*nv]@x)**2*cov_I 
    #sum of residuals
    sumJ = 0
    for i in range(np.size(Jz,0)):
        sumJ = sumJ + Jz[i]
    #this is code used to create a csv of an array for debugging
    #file = open(r'data1.csv','w+',newline='')
    #with file:
    #   write = csv.writer(file)
    #   write.writerows(C)
    return(x,S_d,J,sumJ)
while True:
    cnx = mysql.connector.connect(user='root', password='sqldba', host='localHost', database='openPDC')
    cursor = cnx.cursor()
    connTime = time.time()

	# determine the most recent record when the script starts
    query = "SELECT tsmID, Timestamp FROM timeseriesmeasurement WHERE tsmID = (SELECT MAX(tsmID) FROM timeseriesmeasurement)"
    cursor.execute(query)
    results = cursor.fetchall()
    cnx.commit() # commit after every query

    lastRecord = str(results[0][0]) # this is the most recent tsmID...
    lastTime = str(results[0][1])	# ... and its timestamp
    initQueryTime = time.time()
    if timing:
    	print("Initial Import Time: " + str(importTime - startTime))
    	print("Connection Time: " + str(connTime - importTime))
    	print("Initial Query Took: " + str(initQueryTime - connTime))
    	print("Initial Setup Took: " + str(initQueryTime - startTime))
    outF = open(argv[2]+'.csv','w')
    if verbose:
    	print("==============================================================")
    	print("Starting at: " + lastTime + "(" + lastRecord + ")")
    	print("==============================================================")
    counter = 0
    startLoopTime = time.time()
    while True:
		# Build the query using the latest lastRecord value (VHPA measurements)
        query ='''Select timeseriesmeasurement.tsmID, timeseriesmeasurement.SignalID, timeseriesmeasurement.Timestamp, timeseriesmeasurement.Value from timeseriesmeasurement 
			where SignalID in ('b8413f53-b191-11e6-a0b8-001c23c619e5','f4d95912-b196-11e6-a0b8-001c23c619e5','289de3cc-36cb-11e9-a69c-001c23c619e5','59a967b5-8e87-11e8-ac5c-001c23c619e5','eb63e8b6-8f42-11e8-a305-001c23c619e5','289be161-36cb-11e9-a69c-001c23c619e5','59a756e7-8e87-11e8-ac5c-001c23c619e5','4a086cb7-b196-11e6-a0b8-001c23c619e5','c4e90889-b196-11e6-a0b8-001c23c619e5','58b6eab2-b197-11e6-a0b8-001c23c619e5','4a0ab63e-b196-11e6-a0b8-001c23c619e5','b84319e0-b191-11e6-a0b8-001c23c619e5','58b8918c-b197-11e6-a0b8-001c23c619e5','f4db6c35-b196-11e6-a0b8-001c23c619e5','59ae35d3-8e87-11e8-ac5c-001c23c619e5','59ac4109-8e87-11e8-ac5c-001c23c619e5','28a1f58b-36cb-11e9-a69c-001c23c619e5','eb68d0c9-8f42-11e8-a305-001c23c619e5','eb66ac7a-8f42-11e8-a305-001c23c619e5','c4eae0a2-b196-11e6-a0b8-001c23c619e5','28a04068-36cb-11e9-a69c-001c23c619e5','4a0b5530-b196-11e6-a0b8-001c23c619e5','b84391dc-b191-11e6-a0b8-001c23c619e5','58b8fcf6-b197-11e6-a0b8-001c23c619e5','f4dd92fc-b196-11e6-a0b8-001c23c619e5','59aec0a3-8e87-11e8-ac5c-001c23c619e5','59acdda6-8e87-11e8-ac5c-001c23c619e5','28a2f038-36cb-11e9-a69c-001c23c619e5','eb697be9-8f42-11e8-a305-001c23c619e5','eb6736f1-8f42-11e8-a305-001c23c619e5','28a0cf5a-36cb-11e9-a69c-001c23c619e5','c4eb4fb0-b196-11e6-a0b8-001c23c619e5','b841bbb7-b191-11e6-a0b8-001c23c619e5','f4da1185-b196-11e6-a0b8-001c23c619e5','289e86f0-36cb-11e9-a69c-001c23c619e5','59aa1eba-8e87-11e8-ac5c-001c23c619e5','eb64a185-8f42-11e8-a305-001c23c619e5','289c7829-36cb-11e9-a69c-001c23c619e5','59a7d0c0-8e87-11e8-ac5c-001c23c619e5','4a09196b-b196-11e6-a0b8-001c23c619e5')
			 AND tsmID > {} ORDER BY timeseriesmeasurement.Timestamp,timeseriesmeasurement.tsmID;'''.format(str(lastRecord))
    	# Execute query
        cursor.execute(query)
    	# For every measurement retrieved from the query...
        outList = []
        for (tsmID, signalID, timestamp, value) in cursor:
        	if verbose:
        		print(tsmID,signalID,timestamp,value)
        	outList.append(str(tsmID) + ',' + signalID + ',' + str(timestamp) + ',' + str(value))
        	lastRecord = str(tsmID)
    
        outString = ""
        for item in outList:
        	outString += item + "\n"
        counter += len(outList)
        if debug:
        	outString += "Counter: " + str(counter) + "\n"
        outF.write(outString)
        cnx.commit() # commit after every query
        if counter >= int(argv[1]):
            outF.close()
            end = time.time()
            if verbose:
                print(end-startTime)
            break
    	#time.sleep(.5) # sleep for 500 ms before fetching the next batch of data   
    endLoopTime = time.time()
    qTime = endLoopTime - startLoopTime
    if innerTime:
        print("Query Data took: " + str(qTime))
#lineCountTime = time.time()
#call('wc -l out.csv',shell=True)
#lcTime = time.time() - lineCountTime
#if innerTime:
#	print("Line count took: " + str(lcTime))
    detectionTime = time.time()    
    LSEStartTime = time.time()
    f = open('out.csv', 'r')
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
            lse,S_d,J,sumJ = LinearSE(impedances,voltages,currents,cov_v,cov_I)
            LSEStopTime = time.time()
            t = LSEStopTime-LSEStartTime
            LSEStartTime = time.time()
            #print("Here is the LSE: \n This is x: \n {0} \n This is S_d: \n {1} \n This is J: \n {2}\n This is Jsum: \n {3} \n That run took {4} seconds.\n".format(str(lse),str(S_d),str(J),str(sumJ),str(t)))
            #print(str(voltages))
            #print(str(currents))
        #LSE
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
