# findpeak
findpeak from data

#mainprogram7&kit7: 1. can fit multiple peaks in one data window

#mainprogram4&kit4: 1. fit the highest peak and minus it, until the remain is covered by noise.

#myplot: to plot several types of curve easier

#lmfittest1: this .py shows why lmfit is not robust for finding peaks. Run this program 10 times you will find the results are different.

#tips:
1.the .txt file of data must be several columns of data, please delete the title. the default parameters are: 1st column: x data, 4th column: y data.

2.different confinements will generate different results:
  r: signal to noise ratio
  w: critical frequency of low pass filter
  (the non-peaks interval)
