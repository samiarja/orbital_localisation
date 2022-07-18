import loris
import numpy as np
import matplotlib.pyplot as plt 
import os
import fnmatch
from tqdm import tqdm
import numpy as np
import scipy.io as sio

events = []
FILENAME = "projectile_nospeed"
mat = sio.loadmat('./test_files/' + FILENAME + '.mat')

events = mat['e']
print(events[0][0]["x"].shape)

matX  =  events[0][0]["x"]
matY  =  events[0][0]["y"]
matP  =  events[0][0]["p"]
matTs = events[0][0]["t"]
# matvx = events[0][0]["vx"]
# matvy = events[0][0]["vy"]
print(matX)

nEvents = events[0][0]["x"].shape[1]
x  = matX.reshape((nEvents, 1))
y  = matY.reshape((nEvents, 1))
p  = matP.reshape((nEvents, 1))
ts = matTs.reshape((nEvents, 1))

events = np.zeros((nEvents,4))

events = np.concatenate((ts,x, y, p),axis=1).reshape((nEvents,4))

finalArray = np.asarray(events)
print(finalArray)
# finalArray[:,0] -= finalArray[0,0]

ordering = "txyp"
loris.write_events_to_file(finalArray, "test_files/" + FILENAME + ".es",ordering)
print("File: " + FILENAME + "converted to .es")
