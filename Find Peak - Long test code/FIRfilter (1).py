import sys
import csv
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy import signal


def plot_response(fs, w, h, title):
    "Utility function to plot response functions"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(0.5*fs*w/np.pi, 20*np.log10(np.abs(h)))
    ax.set_ylim(-40, 5)
    ax.set_xlim(0, 0.5*fs)
    ax.grid(True)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Gain (dB)')
    ax.set_title(title)

#filename = sys.argv[1]

filename = "Pulse_20200831-110833_data.txt"
with open(filename,"r") as r:
    r = r.read().split("\n")
    raw1= json.loads(r[1])

data1 = np.array(raw1)
freq = data1[0]
amp = data1[3]

#------------------------------------------------
# Artifically change the x-axis to time and pre
#------------------------------------------------

fake_time = freq*1e-9
fake_time_step = fake_time[1]-fake_time[0]
fake_fs = (1/fake_time_step)

#------------------------------------------------
# FFT the raw data and create a "fake" freq axis
#------------------------------------------------

##N = len(fake_time)-1
##f_step = fake_fs/N
##amp_ft = np.fft.fft(amp) # only for plotting
##fake_freq_axis = np.arange(0,fake_fs+f_step,f_step) # only for plotting
##
#------------------------------------------------
# Create a FIR filter and apply it to input signal
#------------------------------------------------

fs = fake_fs        # Sample rate, Hz
cutoff = 1000    # Desired cutoff frequency, Hz
trans_width = 160  # Width of transition from pass band to stop band, Hz
numtaps = 20      # Size of the FIR filter.

# create the taps (coefficients) 
taps = signal.remez(numtaps,
                    [0, cutoff, cutoff+ trans_width, 0.5*fs],
                    [1, 0], Hz=fs)

# Apply the filter
amp_filtered = signal.lfilter(taps, 1.0, amp) # apply the filter to the signal

#------------------------------------------------
# Trim the delay at the head of the filtered signal.
#------------------------------------------------

freq_trim = freq[numtaps+1::] 
amp_filtered_trim = amp_filtered[numtaps+1::]


#------------------------------------------------
# Calculate the freq resonpose and plot it
#------------------------------------------------
w, h = signal.freqz(taps, [1], worN=2000)
plot_response(fs, w, h, "Low-pass Filter") 

#------------------------------------------------
# Plot input and output signal
#------------------------------------------------

plt.figure()
plt.plot(freq,amp)

##plt.figure()
##plt.plot(fake_freq_axis,abs(amp_ft))

plt.figure()
plt.plot(freq_trim,amp_filtered_trim)
##plt.ylim(0.08, 0.3)

plt.show()



                
            






    

