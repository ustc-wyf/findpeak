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

def DW_cal(array, smooth):
    n = len(array)
    numerator = 0
    denominator = 0
    for i in range(n):
        if i == 0:
            numerator = numerator + 0
        else:
            numerator = numerator + ((array[i]-smooth[i]) - (array[i-1] - smooth[i-1]))**2
        denominator = denominator + (array[i]-smooth[i])**2
    return numerator/denominator*n/(n-1)

def smooth_al(array):
    wd = 5
    optimize = True
    DW_min = 5
    while optimize == True:
        smooth = savgol_filter(array, wd, 2)
        DW = DW_cal(array, smooth)
        if abs(2 - DW) < DW_min:
            wd = wd + 2
            DW_min = abs(2 - DW)
        else:
            wd = wd - 2
            smooth = savgol_filter(array, wd, 2)
            DW = DW_cal(array, smooth)
            break
    return smooth, wd

def noise_level(array):
    h_size = len(array) - 2
    h_value = np.zeros(h_size)
    for i in range(1, h_size - 1):
        h_value[i] = abs(array[i] - array[i-1]) + abs(array[i] - array[i+1])
    h_value = h_value/2
    return np.median(h_value)

def RMS(array):
    result = 0
    n = len(array)
    for i in range(n):
        result = result + array[i]**2
    return np.sqrt(result/n)

def find_peak_2nd(array, line):
    index = []
    for i in range(len(array)-1):
        check = (array[i+1] - line)*(array[i] - line)
        if check < 0:
            if array[i] > line:
                index.append(i) 
            else:
                index.append(i+1)
    return index
    # Find peak from intersection. If the peak between intersection is negative, the peak is selected
    peak_index = []
    for i in range(len(index) - 1):
        if index[i+1] - index[i] > 1:
            #check_index = int((index[i+1] + index[i])/2)
            check_index = index[i] + np.argmin(array[index[i]:index[i+1]])
            print(index[i], check_index, index[i+1])
            check_value = array[check_index]
            if check_value < line:
                peak_index.append(check_index)

def check_th(array, peak_index ,line):
    peak_index_2nd = []
    for i in range(len(peak_index)):
        if array[peak_index[i]] > line:
            peak_index_2nd.append(peak_index[i])
    return peak_index_2nd

#%% Get data
filename1 = "Pulse_20200913-162824_data.txt"
with open(filename1,"r") as r:
    r = r.read().split("\n")
    raw1= json.loads(r[1])

filename2 = "vna_cont_spec_pulselike.txt"
with open(filename2,"r") as r:
    r = r.read().split("\n")
    raw2= json.loads(r[1])

data = np.array(raw1)
# Process data to move it to 0 and flip the dip into peak
# Since the y-axis is arbitrary unit, this action can be done without losing information
x_data = data[3]
x_data = x_data - max(x_data)
x_data = abs(x_data)
#%% Find good smoothing parameters
data_sm, wd = smooth_al(x_data)
# Plot data and smoothed data
fig1, ax1 = plt.subplots(figsize= (12,7))
ax1.plot(x_data, color = "C0")
ax1.plot(data_sm, color = "C1")
#ax1.plot(data_2nd_dev*100)
#%% Calculate threshold value
# Threshold for 2nd derivative
data_2nd_dev = savgol_filter(x_data, wd, 2, deriv = 2)
th_h = -5*noise_level(data_2nd_dev) # Calculate the h_value in the paper

#%% Find peak
# Find intersection:
peak_index = find_peak_2nd(data_2nd_dev, th_h)
# Plot founded peak
fig2, ax2 = plt.subplots(figsize= (12,7))
ax2.plot(data_2nd_dev, color = "C0")
ax2.axhline(th_h, color = "C1")
for i in range(len(peak_index)):
    ax2.scatter(peak_index[i], data_2nd_dev[peak_index[i]])
#%% Check 2nd threshold
# Threshold for smoothed data
th_1 = RMS(data_sm)/2 # The criterial here can be varied but should keep fixed between codes
peak_index_2nd = check_th(data_sm, peak_index, th_1)
# Plot peak that sastify the 2nd threshold
fig3, ax3 = plt.subplots(figsize = (12,7))
ax3.plot(data_sm, color = "C0")
ax3.axhline(th_1, color = "C1")
for i in range(len(peak_index_2nd)):
    ax3.scatter(peak_index_2nd[i], data_sm[peak_index_2nd[i]])
#%% Check user input threshold - input by the user
th_u = th_1 * 2
peak_index_user = check_th(data_sm, peak_index_2nd, th_u)
# Plot peak that sastify user defined threshold
fig4, ax4 = plt.subplots(figsize = (12,7))
ax4.plot(data_sm, color = "C0")
ax4.axhline(th_u, color = "C1")
for i in range(len(peak_index_user)):
    ax4.scatter(peak_index_user[i], data_sm[peak_index_user[i]])
#%%
        
