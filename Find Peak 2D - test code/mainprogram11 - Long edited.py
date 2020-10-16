#%%import stuff
import json
#from tkinter import filedialog
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.signal import savgol_filter
from lmfit import Parameters
from lmfit.models import GaussianModel,LorentzianModel
#%%defination
def p2(a,miu,s,zp,bottom,top,n,count):
    '''plot a Guassian distribution'''
    width = (top - bottom)/(n - 1)
    x = []
    y = []
    for i in range(n):
        x.append(bottom + i*width)
        ytemp = zp[i] + a*math.exp( - (bottom + width*i - miu)*(bottom + width*i - miu)/(2*s*s))
        y.append(ytemp)
    plt.plot(x,y,label = 'Gaussian %d'%count)

def p3(a,miu,g2,zp,bottom,top,n,count):
    '''plot a Lorentzian distribution'''
    width = (top - bottom)/(n - 1)
    x = []
    y = []
    for i in range(n):
        x.append(bottom + i*width)
        ytemp = zp[i] + a/((bottom + width*i - miu)*(bottom + width*i - miu) + g2)
        y.append(ytemp)
    plt.plot(x,y,label = 'Lorentzian %d'%count)

def loadData(infile,k):
    """
    loadData will extract data from a txt file to a python list
    infile: the file pointer
    k: the ordinal number of row you want
    return dataset
    dataset: the colomn you want
    """
    f = open(infile,'r')
    #f = f.read().split("\n")
    #raw = json.loads(f[1])
    f = f.read()
    raw = json.loads(f)
    data = np.array(raw)
    dataset = data[k]
    return dataset          

def loadData2D(infile):
    '''
    Parameters
    ----------
    infile : string
        File name of the target file.

    Returns
    -------
    x_data : array
        x axis.
    y_data : array
        y axis.
    dataset : ndarray
        2D mao data.

    '''
    f = open(infile,'r')
    #f = f.read().split("\n")
    #raw = json.loads(f[1])
    f = f.read()
    raw = json.loads(f)
    y_data = [*raw.keys()]# Unpack into a list literal, get y data (frequency)
    x_data = np.array(raw[y_data[0]][0]) #Get x data (time)
    image_shape = [len(y_data), len(x_data)]
    dataset = np.zeros(image_shape)
    for i in range(len(y_data)):
        dataset[i,:] = np.array(raw[y_data[i]][1])
    return x_data, y_data, dataset     
 
def DW_cal(data, data_sm):
    """DW_cal can calculate DW value of the original data and smoothed data
    input
    data: array, original data
    data_sm: array, smoothed data
    return
    DW value: double
    """
    n = len(data)
    numerator = 0
    denominator = 0
    for i in range(n):
        if i == 0:
            numerator = numerator + 0
        else:
            numerator = numerator + ((data[i] - data_sm[i]) - (data[i-1] - data_sm[i-1]))**2
        denominator = denominator + (data[i] - data_sm[i])**2
    return numerator/denominator*n/(n - 1)

def smooth_al(data):
    """smooth_al can find an appropriate window size of SG algorithm and return smoothed data
    data: array, original data
    return
    smooth: array, smoothed data
    wd: int, window size
    """
    wd = 5
    optimize = True
    DW_min = 5
    while optimize == True:
        smooth = savgol_filter(data, wd, 2)
        DW = DW_cal(data, smooth)
        if abs(2 - DW) < DW_min:
            wd = wd + 2
            DW_min = abs(2 - DW)
        else:
            wd = wd - 2
            smooth = savgol_filter(data, wd, 2)
            DW = DW_cal(data, smooth)
            break
    return smooth, wd

def noise_level(data):
    """noise_level calculate threshold of a set of array
    data: array, original data
    data_sm: array, smoothed data, must be equal length with data
    return
    90%biggest value of the deviation between data and smoothed data
    """
    length=len(data) - 2
    dev=[]
    for i in range(1,length - 1):
        dev.append((abs(data[i] - data[i-1]) + abs(data[i] - data[i + 1]))/2)
    dev.sort()
    return dev[round(0.9*length)]

def cutter(data,wf):
    length = len(data)
    absdata = []
    for i in range(length):
        absdata.append(abs(data[i]))
    l = max(absdata)*1.01# multiplus 1.01 to avoid index out of range
    step = l/100
    p = []#point number
    pd = []#point density
    for i in range(100):
        p.append(0)
    for i in range(length):
        p[int(absdata[i]/step)] = p[int(absdata[i]/step)] + 1
    for i in range(100):
        sum = 0
        for j in range(i + 1):
            sum = sum + p[j]
        pd.append(sum/(step*(i + 1) + l))#default parameter
# =============================================================================
#     plt.plot(range(len(pd)),pd)
#     plt.title('η(w) - w')
#     plt.xlabel('w')
#     plt.ylabel('η')
#     plt.show()
# =============================================================================
    mindex = pd.index(max(pd))
    return wf*(mindex + 1)*step#weight factor

def RMS(x,y,temp):
    """RMS:root mean square
    """
    result = 0
    n = len(x)
    for i in range(n):
        result = result + (y[i] - np.polyval(temp,x[i]))**2
    return np.sqrt(result/n)

def nrf(data,thr):
    """noise region finding"""
    pr = []#peak region
    flag = 0
    length = len(data)
    for i in range(length):
        if(flag==0 and abs(data[i])>thr):
            flag = 1
            bottom = i
        elif(flag==1 and abs(data[i])<=thr):
            flag = 0
            top = i - 1
            for j in range(max(0,2*bottom - top),min(length,2*top - bottom)):
                pr.append(j)
        elif(flag==1 and i==length - 1):
            top = i
            for j in range(max(0,2*bottom - top),min(length,2*top - bottom)):
                pr.append(j)
    region = list(set(range(len(data))) - set(pr))
    return region

def noise_fitting(x,data,wd,wf):
    length = len(data)
    data_1st_dev = savgol_filter(data, wd, 2, deriv = 1)
    data_trans = []
    avg = sum(data)/len(data)
    for i in range(length):
        data_trans.append(data[i] - avg)
    thr1 = cutter(data_trans,wf)
    thr2 = cutter(data_1st_dev,wf)
    z1 = []
    z2 = []
    z3 = []
    z4 = []
    for i in range(len(data_1st_dev)):
        z1.append(thr1)
        z2.append( - thr1)
        z3.append(thr2)
        z4.append( - thr2)
# =============================================================================
#     plt.plot(range(len(data_trans)),data_trans)
#     plt.plot(range(len(z1)),z1)
#     plt.plot(range(len(z2)),z2)
#     plt.title('data and threshold')
#     plt.show()
#     plt.plot(range(len(data_1st_dev)),data_1st_dev)
#     plt.plot(range(len(z3)),z3)
#     plt.plot(range(len(z4)),z4)
#     plt.title('data_1st_dev and threshold')
#     plt.show()
# =============================================================================
    #print(thr1,thr2)
    nx = []
    ny = []
    region1 = nrf(data_trans,thr1)
    region2 = nrf(data_1st_dev,thr2)
    region = list(set(region1)&set(region2))
    for i in range(len(region)):
        nx.append(x[region[i]])
        ny.append(data[region[i]])
# =============================================================================
#     plt.plot(nx,ny)
#     plt.title('non-peak section')
#     plt.show()
# =============================================================================
    flag = 1
    degree = 3
    while(flag==1):
        temp1 = np.polyfit(nx,ny,degree)
        temp2 = np.polyfit(nx,ny,degree + 1)
        rms1 = RMS(nx,ny,temp1)
        rms2 = RMS(nx,ny,temp2)
        if(rms2/rms1>=0.98):
            flag = 0
        else:
            degree = degree + 1
    """rms=[]
    for i in range(10):
        temp=np.polyfit(nx,ny,i)
        rms.append(RMS(nx,ny,temp))
    plt.plot(range(1,11),rms)
    plt.xlabel('n')
    plt.ylabel('rms(n)')
    plt.show()"""
    #print(degree)
    temp = np.polyfit(nx,ny,degree)
    zeroline = []
    tempdev = []
    for i in range(length):
        zeroline.append(np.polyval(temp,x[i]))
    for i in range(len(nx)):
        tempdev.append(abs(np.polyval(temp,nx[i]) - ny[i]))
    tempdev.sort()
    noise = tempdev[int(0.98*len(nx))]
    return zeroline,noise,degree
    
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
    length = len(data)
    width = (top - bottom)/(length - 1)
    absdata = []
    peak = []
    for i in range(length):
        absdata.append(abs(data[i] - zp[i]))
    i = 0
    fsnum = 0#fitting section number
    fsec = []#fitting section
    fdata = []#signal of fitting section
    fnum = []#fitting number
    fm = []#index of max and min points
    while(i<length):
        if absdata[i]>noise*r:
            fsnum = fsnum + 1
            fsec.append([])
            fdata.append([])
            tempmax = absdata[i]
            tempmin = absdata[i]
            inma = i
            inmi = i
            fnum.append(0)
            fm.append([])
            direction = 1#1:rising,0:descending
            while(absdata[i]>noise*r):
                if direction==1:
                    if absdata[i]>tempmax:
                        tempmax = absdata[i]
                        inma = i
                    elif absdata[i]<tempmax - noise*r:
                        direction = 0
                        fm[fsnum - 1].append([inma,inmi])
                        tempmin = absdata[i]
                        inmi = i
                        fnum[fsnum - 1] = fnum[fsnum - 1] + 1
                elif direction==0:
                    if absdata[i]<tempmin:
                        tempmin = absdata[i]
                        inmi = i
                    elif absdata[i]>tempmin + noise*r:
                        direction = 1
                        tempmax = absdata[i]
                        inma = i
                fsec[fsnum - 1].append(bottom + width*i)
                fdata[fsnum - 1].append(absdata[i])
                i = i + 1
                if i>=length:
                    break
            if fm[fsnum - 1]==[]:
                del fsec[fsnum - 1]
                del fdata[fsnum - 1]
                del fnum[fsnum - 1]
                del fm[fsnum - 1]
                fsnum = fsnum - 1
        i = i + 1
    for i in range(fsnum):
        pars = Parameters()
        j = 0
        mod = GaussianModel(prefix = 'g1_')
        pars.update(GaussianModel(prefix = 'g%i_'%(j + 1)).make_params())
        sigma0 = math.sqrt((width*(fm[i][j][0] - fm[i][j][1]))**2/(2*math.log(absdata[fm[i][j][0]]/absdata[fm[i][j][1]])))
        pars['g%i_center'%(j + 1)].set(value = bottom + width*fm[i][j][0],min = fsec[i][0],max = fsec[i][ - 1])
        pars['g%i_sigma'%(j + 1)].set(value = sigma0,min = sigma0/20,max = sigma0*20)
        pars['g%i_amplitude'%(j + 1)].set(value = absdata[fm[i][j][0]]/0.3989423*sigma0,min = noise*r/0.3989423*sigma0,max = absdata[fm[i][j][0]]*20/0.3989423*sigma0)
        for j in range(1,fnum[i]):
            mod = mod + GaussianModel(prefix = 'g%i_'%(j + 1))
            pars.update(GaussianModel(prefix = 'g%i_'%(j + 1)).make_params())
            sigma0 = math.sqrt((width*(fm[i][j][0] - fm[i][j][1]))**2/(2*math.log(absdata[fm[i][j][0]]/absdata[fm[i][j][1]])))
            pars['g%i_center'%(j + 1)].set(value = bottom + width*fm[i][j][0],min = fsec[i][0],max = fsec[i][-1])
            pars['g%i_sigma'%(j + 1)].set(value = sigma0,min = sigma0/20,max = sigma0*20)
            pars['g%i_amplitude'%(j + 1)].set(value = absdata[fm[i][j][0]]/0.3989423*sigma0,min = noise*r/0.3989423*sigma0,max = absdata[fm[i][j][0]]*20/0.3989423*sigma0)
# =============================================================================
#         result = mod.fit(fdata[i],pars,x = fsec[i])
#         #print(result.fit_report())
#         plt.plot(fsec[i],fdata[i],'bo',label = 'original')
#         plt.plot(fsec[i],result.best_fit,'r-',label = 'fitting')
#         plt.title('Gaussian fitting')
#         plt.show()
# =============================================================================
        tempbo = int((fsec[i][0] - bottom)/width)
        tempto = int((fsec[i][-1] - bottom)/width)
        for k in range(fnum[i]):
            amplitude = pars['g%i_height'%(k + 1)].value
            sigma = pars['g%i_sigma'%(k + 1)].value
            miu = pars['g%i_center'%(k + 1)].value
            sum1 = 0
            for p in range(tempbo,tempto + 1):
                v = abs(amplitude*math.exp( - (bottom + width*p - miu)*(bottom + width*p - miu)/(2*sigma*sigma)))
                sum1 = sum1 + (v - absdata[k])*(v - absdata[k])
            sum1 = sum1/(tempto - tempbo + 1)
            peak.append([sigma,miu,amplitude,sum1,tempbo,tempto])
    return peak

def findpeakl(data, zp, noise, bottom, top, r):
    """findpeakl can find domain of a peak, then fit it by Lorenztian Curve
    and give gama.
    zp: zeropoint
    bottom: Minimum value of independent variable
    top: Maximum value of independent variable
    r:signal to noise ratio
    return peak
    peak[i][0]: gama2
    peak[i][1]: miu
    peak[i][2]: amplitude
    peak[i][3]: mean square error
    peak[i][4]: left side of the fitting interval
    peak[i][5]: right side of the fitting interval
    """
    length = len(data)
    width = (top - bottom)/(length - 1)
    absdata = []
    peak = []
    for i in range(length):
        absdata.append(abs(data[i] - zp[i]))
    i = 0
    fsnum = 0#fitting section number
    fsec = []#fitting section
    fdata = []#signal of fitting section
    fnum = []#fitting number
    fm = []#index of max and min points
    while(i<length):
        if absdata[i]>noise*r:
            fsnum = fsnum + 1
            fsec.append([])
            fdata.append([])
            tempmax = absdata[i]
            tempmin = absdata[i]
            inma = i
            inmi = i
            fnum.append(0)
            fm.append([])
            direction = 1#1:rising,0:descending
            while(absdata[i]>noise*r):
                if direction==1:
                    if absdata[i]>tempmax:
                        tempmax = absdata[i]
                        inma = i
                    elif absdata[i]<tempmax - noise*r:
                        direction = 0
                        fm[fsnum - 1].append([inma,inmi])
                        tempmin = absdata[i]
                        inmi = i
                        fnum[fsnum - 1] = fnum[fsnum - 1] + 1
                elif direction==0:
                    if absdata[i]<tempmin:
                        tempmin = absdata[i]
                        inmi = i
                    elif absdata[i]>tempmin + noise*r:
                        direction = 1
                        tempmax = absdata[i]
                        inma = i
                fsec[fsnum - 1].append(bottom + width*i)
                fdata[fsnum - 1].append(absdata[i])
                i = i + 1
                if i>=length:
                    break
            if fm[fsnum - 1]==[]:
                del fsec[fsnum - 1]
                del fdata[fsnum - 1]
                del fnum[fsnum - 1]
                del fm[fsnum - 1]
                fsnum = fsnum - 1
        i = i + 1
    for i in range(fsnum):
        pars = Parameters()
        j = 0
        mod = LorentzianModel(prefix = 'l1_')
        pars.update(LorentzianModel(prefix = 'l%i_'%(j + 1)).make_params())
        sigma0 = abs(width*(fm[i][j][0] - fm[i][j][1]))/math.sqrt(absdata[fm[i][j][0]]/absdata[fm[i][j][1]] - 1)
        pars['l%i_center'%(j + 1)].set(value = bottom + width*fm[i][j][0],min = fsec[i][0],max = fsec[i][ - 1])
        pars['l%i_sigma'%(j + 1)].set(value = sigma0,min = sigma0/20,max = sigma0*20)
        pars['l%i_amplitude'%(j + 1)].set(value = absdata[fm[i][j][0]]*sigma0/0.3183099,min = noise*r*sigma0/0.3183099,max = absdata[fm[i][j][0]]*20*sigma0/0.3183099)
        for j in range(1,fnum[i]):
            mod = mod + LorentzianModel(prefix = 'l%i_'%(j + 1))
            pars.update(LorentzianModel(prefix = 'l%i_'%(j + 1)).make_params())
            sigma0 = abs(width*(fm[i][j][0] - fm[i][j][1]))/math.sqrt(absdata[fm[i][j][0]]/absdata[fm[i][j][1]] - 1)
            pars['l%i_center'%(j + 1)].set(value = bottom + width*fm[i][j][0],min = fsec[i][0],max = fsec[i][ - 1])
            pars['l%i_sigma'%(j + 1)].set(value = sigma0,min = sigma0/20,max = sigma0*20)
            pars['l%i_amplitude'%(j + 1)].set(value = absdata[fm[i][j][0]]*sigma0/0.3183099,min = noise*r*sigma0/0.3183099,max = absdata[fm[i][j][0]]*20*sigma0/0.3183099)
# =============================================================================
#         result = mod.fit(fdata[i],pars,x = fsec[i])
#         #print(result.fit_report())
#         plt.plot(fsec[i],fdata[i],'bo',label = 'original')
#         plt.plot(fsec[i],result.best_fit,'r-',label = 'fitting')
#         plt.title('Lorentzian fitting')
#         plt.show()
# =============================================================================
        tempbo = int((fsec[i][0] - bottom)/width)
        tempto = int((fsec[i][ - 1] - bottom)/width)
        for k in range(fnum[i]):
            gama2 = (pars['l%i_sigma'%(k + 1)].value)**2
            amplitude = pars['l%i_height'%(k + 1)].value*gama2
            miu = pars['l%i_center'%(k + 1)].value
            sum1 = 0
            for p in range(tempbo,tempto + 1):
                v = abs(amplitude/((bottom + width*p - miu)*(bottom + width*p - miu) + gama2))
                sum1 = sum1 + (v - absdata[k])*(v - absdata[k])
            sum1 = sum1/(tempto - tempbo + 1)
            peak.append([gama2,miu,amplitude,sum1,tempbo,tempto])
    return peak

#%%Get data
#filename = filedialog.askopenfilename(filetypes = [('TXT','txt')])
filename = "vna_cont_spec_flux_wider.txt"
x, y, alldata =  loadData2D(filename)
peak_dict = {}
#for k in range(len(y)-1, -1, -1):
for k in range(len(y)):
    data = alldata[k,:]
    length = len(data)
    top = max(x)
    bottom = min(x)
    widthx = (top - bottom)/(length - 1)
    mi = min(data)
    ma = max(data)
    # =============================================================================
    # r = float(input("input the signal to noise ratio(reference value: 1.5~2)\n"))
    # wf = float(input("input the weight factor of the width of non-peak region(reference value: 1, if some peaks are considered as noise, turn down this factor)\n"))
    # =============================================================================
    r = 1.5
    wf = 1
    #%% Find good smoothing parameters
    data_sm, wd  =  smooth_al(data)
    # =============================================================================
    # plt.plot(x,data,label = 'original data')
    # plt.plot(x,data_sm,label = 'smoothed data')
    # plt.legend()
    # plt.title('original data and smoothed data(SG algorithm)')
    # plt.show()
    # =============================================================================
    # zero point and noise
    zeroline,noise,degree = noise_fitting(x, data_sm, wd, wf)
    #print('noise level is %f'%noise)
    # =============================================================================
    # plt.plot(x,data_sm,label = 'smoothed data')
    # plt.plot(x,zeroline,label = 'zero line')
    # plt.legend()
    # plt.title('smoothed data and zeropoint')
    # plt.show()
    # =============================================================================
    #display result
    peak1 = findpeakg(data,zeroline,noise,bottom,top,r)
    peak2 = findpeakl(data,zeroline,noise,bottom,top,r)
# =============================================================================
#     print('noise level: %f'%noise)
#     print('polynomial degree: %d'%degree)
# =============================================================================
    m = len(peak1)
    n = len(peak2)
    mse_gauss = 0
    mse_lorent = 0
    for i in range(m):
        mse_gauss += peak1[i][3]
    for i in range(n):
        mse_lorent += peak2[i][3]
    # =============================================================================
    # plt.plot(x,data,label = 'data')
    # plt.show()
    # plt.plot(x,data,label = 'data')
    # =============================================================================
    # =============================================================================
    # for i in range(m):
    #     tempbo = peak1[i][4]
    #     tempto = peak1[i][5]
    #     a1 = peak1[i][2]
    #     miu1 = peak1[i][1]
    #     sigema = peak1[i][0]
    #     mse1 = peak1[i][3]
    #     print("the peak%d:\nGaussian: A:%f, miu:%f, sigma:%f, mean square error:%f"%(i + 1,a1,miu1,sigema, mse1))
    #     p2(-a1,miu1,sigema,zeroline,bottom + widthx*tempbo,bottom + widthx*tempto,tempto - tempbo + 1,i + 1)
    #     plt.plot(x[tempbo:tempto + 1],data_sm[tempbo:tempto + 1],label = 'original')
    #     plt.legend()
    #     plt.show()
    # for i in range(n):
    #     tempbo = peak2[i][4]
    #     tempto = peak2[i][5]
    #     a2 = peak2[i][2]
    #     miu2 = peak2[i][1]
    #     mse2 = peak2[i][3]
    #     gama2 = peak2[i][0]
    #     if gama2>0:
    #         print("the peak%d:\nLorentzian: A:%f, miu:%f: gama:%f, mean square error:%f"%(i + 1,a2/gama2,miu2,math.sqrt(gama2),mse2))
    #     else:
    #         print("the peak%d:\nLorentzian: A:%f, miu:%f: gama2:%f, mean square error:%f"%(i + 1,a2/gama2,miu2,gama2,mse2))
    #     p3(-a2,miu2,gama2,zeroline,bottom + widthx*tempbo,bottom + widthx*tempto,tempto - tempbo + 1,i + 1)
    #     plt.plot(x[tempbo:tempto + 1],data_sm[tempbo:tempto + 1],label = 'original')
    #     plt.legend()
    #     plt.show()
    # =============================================================================
    peak_list = []
    if mse_gauss < mse_lorent:
        for i in range(m):
            a1 = peak1[i][2]
            miu1 = peak1[i][1]
            sigema = peak1[i][0]
            peak_list.append(miu1)
            #plt.scatter(miu1, y[k])
            #p2(-a1,miu1,sigema,zeroline,bottom,top,len(data),i + 1)
    else:
        for i in range(n):
            a2 = peak2[i][2]
            miu2 = peak2[i][1]
            gama2 = peak2[i][0]
            peak_list.append(miu2)
            #plt.scatter(miu2, y[k])
            #p3(-a2,miu2,gama2,zeroline,bottom,top,len(data),i + 1)
    #print("{}".format(y[k]))
    peak_dict["{}".format(y[k])] = peak_list
    # =============================================================================
    # plt.legend()
    # plt.show()
    # =============================================================================
with open('{}-peak_dict.txt'.format(filename), 'w') as f:
    print(peak_dict, file=f)