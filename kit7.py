import sympy
import numpy
import math
import matplotlib.pyplot as plt
import myplot
import tkinter.filedialog as tk
from scipy.optimize import curve_fit
from lmfit import Model,Parameters,report_fit
from lmfit.models import GaussianModel,LorentzianModel

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
    i=0
    fsnum=0#fitting section number
    fsec=[]#fitting section
    fdata=[]#signal of fitting section
    fnum=[]#fitting number
    fm=[]#index of max and min points
    while(i<length):
        if absdata[i]>noise*r:
            fsnum=fsnum+1
            fsec.append([])
            fdata.append([])
            tempmax=absdata[i]
            tempmin=absdata[i]
            inma=i
            inmi=i
            fnum.append(0)
            fm.append([])
            direction=1#1:rising,0:descending
            while(absdata[i]>noise*r):
                if direction==1:
                    if absdata[i]>tempmax:
                        tempmax=absdata[i]
                        inma=i
                    elif absdata[i]<tempmax-noise:
                        direction=0
                        fm[fsnum-1].append([inma,inmi])
                        tempmin=absdata[i]
                        inmi=i
                        fnum[fsnum-1]=fnum[fsnum-1]+1
                elif direction==0:
                    if absdata[i]<tempmin:
                        tempmin=absdata[i]
                        inmi=i
                    elif absdata[i]>tempmin+noise:
                        direction=1
                        tempmax=absdata[i]
                        inma=i
                fsec[fsnum-1].append(bottom+width*i)
                fdata[fsnum-1].append(absdata[i])
                i=i+1
                if i>=length:
                    break
            if fm[fsnum-1]==[]:
                del fsec[fsnum-1]
                del fdata[fsnum-1]
                del fnum[fsnum-1]
                del fm[fsnum-1]
                fsnum=fsnum-1
        i=i+1
    for i in range(fsnum):
        pars=Parameters()
        j=0
        mod=GaussianModel(prefix='g1_')
        pars.update(GaussianModel(prefix='g%i_'%(j+1)).make_params())
        sigma0=math.sqrt((width*(fm[i][j][0]-fm[i][j][1]))**2/(2*math.log(absdata[fm[i][j][0]]/absdata[fm[i][j][1]])))
        pars['g%i_center'%(j+1)].set(value=bottom+width*fm[i][j][0],min=fsec[i][0],max=fsec[i][-1])
        pars['g%i_sigma'%(j+1)].set(value=sigma0,min=sigma0/20,max=sigma0*20)
        pars['g%i_amplitude'%(j+1)].set(value=absdata[fm[i][j][0]]/0.3989423*sigma0,min=noise*r,max=absdata[fm[i][j][0]]*20/0.3989423*sigma0)
        for j in range(1,fnum[i]):
            mod=mod+GaussianModel(prefix='g%i_'%(j+1))
            pars.update(GaussianModel(prefix='g%i_'%(j+1)).make_params())
            sigma0=math.sqrt((width*(fm[i][j][0]-fm[i][j][1]))**2/(2*math.log(absdata[fm[i][j][0]]/absdata[fm[i][j][1]])))
            pars['g%i_center'%(j+1)].set(value=bottom+width*fm[i][j][0],min=fsec[i][0],max=fsec[i][-1])
            pars['g%i_sigma'%(j+1)].set(value=sigma0,min=sigma0/20,max=sigma0*20)
            pars['g%i_amplitude'%(j+1)].set(value=absdata[fm[i][j][0]]/0.3989423*sigma0,min=noise*r,max=absdata[fm[i][j][0]]*20/0.3989423*sigma0)
        result=mod.fit(fdata[i],pars,x=fsec[i])
        print(result.fit_report())
        plt.plot(fsec[i],fdata[i],'bo')
        plt.plot(fsec[i],result.best_fit,'r-')
        plt.show()
        tempbo=int((fsec[i][0]-bottom)/width)
        tempto=int((fsec[i][-1]-bottom)/width)
        for k in range(fnum[i]):
            amplitude=pars['g%i_height'%(k+1)].value
            sigma=pars['g%i_sigma'%(k+1)].value
            miu=pars['g%i_center'%(k+1)].value
            sum1=0
            for p in range(tempbo,tempto+1):
                v=abs(amplitude*math.exp(-(bottom+width*p-miu)*(bottom+width*p-miu)/(2*sigma*sigma)))
                sum1=sum1+(v-absdata[k])*(v-absdata[k])
            sum1=sum1/(tempto-tempbo+1)
            peak.append([sigma,miu,amplitude,sum1,tempbo,tempto])
    return peak

def findpeakl(data, zp, noise, bottom, top, r):
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
    i=0
    fsnum=0#fitting section number
    fsec=[]#fitting section
    fdata=[]#signal of fitting section
    fnum=[]#fitting number
    fm=[]#index of max and min points
    while(i<length):
        if absdata[i]>noise*r:
            fsnum=fsnum+1
            fsec.append([])
            fdata.append([])
            tempmax=absdata[i]
            tempmin=absdata[i]
            inma=i
            inmi=i
            fnum.append(0)
            fm.append([])
            direction=1#1:rising,0:descending
            while(absdata[i]>noise*r):
                if direction==1:
                    if absdata[i]>tempmax:
                        tempmax=absdata[i]
                        inma=i
                    elif absdata[i]<tempmax-noise:
                        direction=0
                        fm[fsnum-1].append([inma,inmi])
                        tempmin=absdata[i]
                        inmi=i
                        fnum[fsnum-1]=fnum[fsnum-1]+1
                elif direction==0:
                    if absdata[i]<tempmin:
                        tempmin=absdata[i]
                        inmi=i
                    elif absdata[i]>tempmin+noise:
                        direction=1
                        tempmax=absdata[i]
                        inma=i
                fsec[fsnum-1].append(bottom+width*i)
                fdata[fsnum-1].append(absdata[i])
                i=i+1
                if i>=length:
                    break
            if fm[fsnum-1]==[]:
                del fsec[fsnum-1]
                del fdata[fsnum-1]
                del fnum[fsnum-1]
                del fm[fsnum-1]
                fsnum=fsnum-1
        i=i+1
    for i in range(fsnum):
        pars=Parameters()
        j=0
        mod=LorentzianModel(prefix='l1_')
        pars.update(LorentzianModel(prefix='l%i_'%(j+1)).make_params())
        sigma0=abs(width*(fm[i][j][0]-fm[i][j][1]))*math.sqrt(fm[i][j][0]/fm[i][j][1]-1)
        pars['l%i_center'%(j+1)].set(value=bottom+width*fm[i][j][0],min=fsec[i][0],max=fsec[i][-1])
        pars['l%i_sigma'%(j+1)].set(value=sigma0,min=sigma0/20,max=sigma0*20)
        pars['l%i_amplitude'%(j+1)].set(value=absdata[fm[i][j][0]]*sigma0/0.3183099,min=noise*r,max=absdata[fm[i][j][0]]*20*sigma0/0.3183099)
        for j in range(1,fnum[i]):
            mod=mod+LorentzianModel(prefix='l%i_'%(j+1))
            pars.update(LorentzianModel(prefix='l%i_'%(j+1)).make_params())
            sigma0=math.sqrt((width*(fm[i][j][0]-fm[i][j][1]))**2/(2*math.log(absdata[fm[i][j][0]]/absdata[fm[i][j][1]])))
            pars['l%i_center'%(j+1)].set(value=bottom+width*fm[i][j][0],min=fsec[i][0],max=fsec[i][-1])
            pars['l%i_sigma'%(j+1)].set(value=sigma0,min=sigma0/20,max=sigma0*20)
            pars['l%i_amplitude'%(j+1)].set(value=absdata[fm[i][j][0]]*sigma0/0.3183099,min=noise*r,max=absdata[fm[i][j][0]]*20*sigma0/0.3183099)
        result=mod.fit(fdata[i],pars,x=fsec[i])
        print(result.fit_report())
        plt.plot(fsec[i],fdata[i],'bo')
        plt.plot(fsec[i],result.best_fit,'r-')
        plt.show()
        tempbo=int((fsec[i][0]-bottom)/width)
        tempto=int((fsec[i][-1]-bottom)/width)
        for k in range(fnum[i]):
            gama2=(pars['l%i_sigma'%(k+1)].value)**2
            amplitude=pars['l%i_height'%(k+1)].value*gama2
            miu=pars['l%i_center'%(k+1)].value
            sum1=0
            for p in range(tempbo,tempto+1):
                v=abs(amplitude/((bottom+width*p-miu)*(bottom+width*p-miu)+gama2))
                sum1=sum1+(v-absdata[k])*(v-absdata[k])
            sum1=sum1/(tempto-tempbo+1)
            peak.append([gama2,miu,amplitude,sum1,tempbo,tempto])
    return peak
