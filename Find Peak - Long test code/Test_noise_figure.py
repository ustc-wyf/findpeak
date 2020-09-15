# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 15:50:56 2020

@author: Hoang Long Nguyen
"""
#%% Import stuff
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def RMS(Arr):
    result = 0
    n = Arr.shape[0]
    for i in range(n):
        result = result + Arr[i]**2
    return np.sqrt(result/n)

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w
#%% Get data
filename1 = "Pulse_20200831-110833_data.txt"
with open(filename1,"r") as r:
    r = r.read().split("\n")
    raw1= json.loads(r[1])

filename2 = "scan_20200831-114154_data.txt"
with open(filename2,"r") as r:
    r = r.read().split("\n")
    raw2= json.loads(r[1])

data1 = np.array(raw1)
data2 = np.array(raw2)
#%% Plot data
plt.plot(data2[0], data2[3])
plt.plot(data1[0], data1[3])
#%% Process data
rms1 = RMS(data1[3])
plt.plot(data1[0],data1[3], color = "C0", label = "data")
plt.axhline(rms1, color = "C1", label = "RMS")
plt.legend()
#%%
l = data1[3].shape[0]
w = int(l*0.01)
mavg = moving_average(data1[3],w)
noise = data1[3,:int(l - w + 1)] - mavg
#np.set_printoptions(precision=2)
plt.plot(data1[3,:int(l - w + 1)], label = "data")
plt.plot(mavg, label = "Moving avg")
plt.plot(noise, label = "noise")
plt.axhline(RMS(noise))
plt.legend()
#%%
smooth = savgol_filter(data1[3], 13, 2, mode="nearest")
dy = np.gradient(smooth)
#plt.plot(data1[3])
plt.plot(dy)