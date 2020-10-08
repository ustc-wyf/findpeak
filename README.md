# findpeak
findpeak from data

#mainprogram11:
1. add an input wf. user can change the width of non-peak region by change this parameter. This parameter is a default number 1 in mainprogram9
2. add some legends and titles.

#mainprogram9:
1. can fit multiple peaks in one data window
2. can find and fit non-peak section automatically
3. skip the first line of .txt file of data(skip the title), default parameters: 1st column: x data, 4th column: y data, change it in line 451&452 of the codes.
4. you only need to input the signal to noise ratio(r): which determine what peaks you want to accept.
5. run without any other .py file

#mainprogram7&kit7: 
1. can fit multiple peaks in one data window

#mainprogram4&kit4: 
1. fit the highest peak and minus it, until the remain is covered by noise.

#myplot: to plot several types of curve easier

#lmfittest1: this .py shows why lmfit is not robust for finding peaks. Run this program 10 times you will find the results are different.

#tips to mainprogram4&7:
1.the .txt file of data must be several columns of data, please delete the title. the default parameters are: 1st column: x data, 4th column: y data.

2.different confinements will generate different results:
  r: signal to noise ratio
  w: critical frequency of low pass filter
  (the non-peaks interval)
