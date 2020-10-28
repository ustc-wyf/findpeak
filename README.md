# findpeak
findpeak from data

#mainprogram13:
1. Fit one maximum with one or two peaks, output the result with lowest chi-square-error.
2. Because there are more overlapped peaks, mainprogram13 adds "reduced chi-square-error" output of the whole fitting section and "mean square error" in the output is not significant now. We can see clearly that the Gaussian model is better than the Lorentzian model!

#mainprogram12:
1. Correct output mistakes.
2. Change estimation of initial value of sigma.

#mainprogram11:
1. Add an input wf. Users can change the width of non-peak region by changing this parameter. This parameter is a default number 1 in mainprogram9
2. Add some legends and titles.

#mainprogram9:
1. Can fit multiple peaks in one data window
2. Can find and fit non-peak section automatically
3. Skip the first line of .txt file of data(skip the title), default parameters: 1st column: x data, 4th column: y data, change it in line 451&452 of the codes.
4. Users only need to input the signal to noise ratio(r): which determines what peaks you want to accept.
5. Run without any other .py file

#mainprogram7&kit7: 
1. Can fit multiple peaks in one data window

#mainprogram4&kit4: 
1. Fit the highest peak and minus it, until the remains are covered by noise.

#myplot: To plot several types of curve easier

#lmfittest1: This .py shows why lmfit is not robust for finding peaks. Run this program 10 times and you will find the results are different.

#tips to mainprogram4&7:
1. The .txt file of data must be several columns of data, please delete the title. the default parameters are: 1st column: x data, 4th column: y data.

2. Different confinements will generate different results:
  r: signal to noise ratio
  w: critical frequency of low pass filter
  (the non-peaks interval)
