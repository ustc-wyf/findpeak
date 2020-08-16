import sympy
import numpy
import math
import sys
import matplotlib.pyplot as plt
import myplot
import tkinter.filedialog as tk
from scipy.optimize import curve_fit

def loadData(infile,k):
    """
    loadData will extract data from a txt file to a python list
    infile: the file pointer
    k: the ordinal number of colomn you want
    return dataset
    dataset: the colomn you want
    """
    f=open(infile,'r')
    sourceInline=f.readlines()
    dataset=[]
    for line in sourceInline:
        temp1=line.strip('\n')
        temp2=temp1.split(' ')
        dataset.append(float(temp2[k-1])) #chose the line you want
    return dataset

def dd(data,n):
    """dd can calculate density distribution of an array
    n: the number of parts you want to divide the range in
    return density
    density: density distribution of data"""
    ma=max(data)
    mi=min(data)
    width=(ma-mi)/n
    length=len(data)
    density=numpy.zeros(n)
    for i in range(length):
        if data[i]!=ma:
            temp=int((data[i]-mi)/width)
            density[temp]=density[temp]+1/length
    myplot.p1(density,mi,ma)
    plt.show()
    return density

def findpeakg(data, zp, noise, bottom, top, r):
    """findpeakg can find domain of a peak, then fit it by Guassian Curve
    and give mean squared error.
    zp: zeropoint
    bottom: Minimum value of independent variable
    top: Maximum value of independent variable
    r:signal to noise ratio
    return peak
    peak[i][0]: sigma
    peak[i][1]: miu
    peak[i][2]: amplitude
    peak[i][3]: mean square error
    peak[i][4]: left side of the fitting interval
    peak[i][5]: right side of the fitting interval
    """
    length=len(data)
    width=(top-bottom)/(length-1)
    absdata=[]
    peak=[]
    for i in range(length):
        absdata.append(abs(data[i]-zp[i]))
    inma=absdata.index(max(absdata))
    ma=max(absdata)
    i=0
    while ma>noise:
        tempx=[]
        tempy=[]
        tempz=[]
        j=inma-1
        tempbo=inma
        if j>=0:
            while absdata[j]>0.8*ma or absdata[j]>noise and absdata[j]-absdata[tempbo]<ma/10:
                if absdata[j]<absdata[tempbo]:
                    tempbo=j
                j=j-1
                if j<0:
                    break
        j=inma+1
        tempto=inma
        if j<length:
            while absdata[j]>0.8*ma or absdata[j]>noise and absdata[j]-absdata[tempto]<ma/10:
                if absdata[j]<absdata[tempto]:
                    tempto=j
                j=j+1
                if j>=length:
                    break
        if tempto-tempbo>1:
            for k in range(tempbo,tempto+1):
                tempx.append(bottom+width*k)
                tempy.append(math.log(absdata[k]))
            print(i+1)
            print(tempx)
            print(tempy)
            tempz=numpy.polyfit(tempx,tempy,2)
            a=tempz[0]
            b=tempz[1]
            c=tempz[2]
            sigema=math.sqrt(-1/(2*a))
            miu=-b/(2*a)
            amplitude=math.exp(c-b*b/(4*a))
            sum1=0
            for k in range(tempbo,tempto+1):
                v=abs(a*math.exp(-(bottom+width*k-miu)*(bottom+width*k-miu)/(2*sigema*sigema)))
                sum1=sum1+(v-absdata[k])*(v-absdata[k])
            sum1=sum1/(tempto-tempbo+1)
            if amplitude>noise*r:
                peak.append([sigema,miu,amplitude,sum1,tempbo,tempto])
            for k in range(length):
                absdata[k]=abs(absdata[k]-amplitude*math.exp(-(bottom+width*k-miu)*(bottom+width*k-miu)/(2*sigema*sigema)))
        else:
            absdata[tempto]=0
            absdata[tempbo]=0
        inma=absdata.index(max(absdata))
        ma=max(absdata)
        i=i+1
    return peak

def findpeakl(data, zp, noise, bottom, top, r):
    """
    find lorentzian peak
    """
    length=len(data)
    width=(top-bottom)/(length-1)
    absdata=[]
    peak=[]
    for i in range(length):
        absdata.append(abs(data[i]-zp[i]))
    inma=absdata.index(max(absdata))
    ma=max(absdata)
    i=0
    while ma>noise:
        tempx=[]
        tempy=[]
        tempz=[]
        j=inma-1
        tempbo=inma
        if j>=0:
            while absdata[j]>0.8*ma or absdata[j]>noise and absdata[j]-absdata[tempbo]<ma/10:
                if absdata[j]<absdata[tempbo]:
                    tempbo=j
                j=j-1
                if j<0:
                    break
        j=inma+1
        tempto=inma
        if j<length:
            while absdata[j]>0.8*ma or absdata[j]>noise and absdata[j]-absdata[tempto]<ma/10:
                if absdata[j]<absdata[tempto]:
                    tempto=j
                j=j+1
                if j>=length:
                    break
        if tempto-tempbo>1:
            for k in range(tempbo,tempto+1):
                tempx.append(bottom+width*k)
                tempy.append(1/absdata[k])
            print(i+1)
            print(tempx)
            print(tempy)
            tempz=numpy.polyfit(tempx,tempy,2)
            a=tempz[0]
            b=tempz[1]
            c=tempz[2]
            gama2=(4*a*c-b*b)/(4*a*a)
            miu=-b/(2*a)
            amplitude=1/a
            sum1=0
            for k in range(tempbo,tempto+1):
                v=abs(a/((bottom+width*k-miu)*(bottom+width*k-miu)+gama2))
                sum1=sum1+(v-absdata[k])*(v-absdata[k])
            sum1=sum1/(tempto-tempbo+1)
            if amplitude/gama2>noise*r:
                peak.append([gama2,miu,amplitude,sum1,tempbo,tempto])
            for k in range(length):
                absdata[k]=abs(absdata[k]-amplitude/((bottom+width*k-miu)*(bottom+width*k-miu)+gama2))
        else:
            absdata[tempto]=0
            absdata[tempbo]=0
        inma=absdata.index(max(absdata))
        ma=max(absdata)
        i=i+1
    return peak
