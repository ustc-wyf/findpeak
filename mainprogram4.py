import kit4
import sys
import matplotlib.pyplot as plt
import myplot
import math
import tkinter
import numpy
from scipy import signal

"""input data"""
filename=tkinter.filedialog.askopenfilename(filetypes=[('TXT','txt')])
data=kit4.loadData(filename,4)
x=kit4.loadData(filename,1)
length=len(data)
top=max(x)
bottom=min(x)
widthx=(top-bottom)/(length-1)
mi=min(data)
ma=max(data)
myplot.p1(data,bottom,top)
plt.show()
r=float(input("input the signal to noise ratio(reference value: 1.5~2)\n"))

"""denoise"""
flagd=input("Do you want to denoise?(which may generate distortion of your peak)\nYes: 1, No: 0\n")
if flagd=='1':
    w=float(input("input the cut-off frequency(Tip: 0<w<1, reference value: 0.2~0.3):\n"))
    fourierb,fouriera=signal.butter(4,w,'low')
    sf=signal.filtfilt(fourierb,fouriera,data)
    myplot.p1(data,0,len(data)-1)
    myplot.p1(sf,0,len(data)-1)
    plt.show()
    for i in range(length):
        data[i]=sf[i]
    mi=min(data)
    ma=max(data)

"""find zero line and noise"""
zeroline=[]
flag1=input("Choose your method to find zero line and noise:\n1. By density distribution: input 1\n2. By assign non-peak area: input 2\n3. By assign zero line and noise: input 3\n")
if flag1 == '1':
    density=kit4.dd(data,50)#manually set parameters: 50
    mden=max(density)
    widthden=(ma-mi)/49
    ze=[]
    for i in range(50):
        ze.append(0)
    zp1=kit4.findpeakg(density,ze,mden/10,mi,ma,1)
    zp2=kit4.findpeakl(density,ze,mden/10,mi,ma,1)
    print(zp1,zp2)
    noise=2.58*zp1[0][0]#a point outside zp+-noise is ?% not because of a noise fluctuation
#68.27%--1sigema, 95%--1.96sigema, 99%--2.58sigema
    zp=zp1[0][1]
    for i in range(length):
        zeroline.append(zp)
    print("\n\n\nzeropoint:%f, noise:%f\n\n\n"%(zp,noise))
elif flag1 =='2':
    i=0
    tempx=[]
    tempy=[]
    print("\n###input selected interval without overlap, print '+' to exit###\n")
    while 1:
        i=i+1
        left=input("input the left end of interval %d:\n"%i)
        if left=='+':
            break
        right=input("input the right end of interval %d:\n"%i)
        if right=='+':
            break
        for j in range(int((float(left)-bottom)/widthx),int((float(right)-bottom)/widthx)):
            if j in range(length):
                tempx.append(x[j])
                tempy.append(data[j])
    degree=input("input the degree of polynomial(reference: times of concavity and convexity transitions + 1)\n")
    temp=numpy.polyfit(tempx,tempy,int(degree))
    devdata=[]
    for i in range(length):
        zeroline.append(numpy.polyval(temp,x[i]))
        devdata.append(abs(zeroline[i]-data[i]))
    devdensity=kit4.dd(data,50)
    devmden=max(devdensity)
    devma=max(devdata)
    devmi=min(devdata)
    widthdev=(devma-devmi)/49
    ze=[]
    for i in range(50):
        ze.append(0)
    devzp1=kit4.findpeakg(devdensity,ze,devmden/10,devmi,devma,1)
    devzp2=kit4.findpeakl(devdensity,ze,devmden/10,devmi,devma,1)
    noise=2.58*devzp1[0][0]#a point outside zp+-noise is ?% not because of a noise fluctuation
#68.27%--1sigema, 95%--1.96sigema, 99%--2.58sigema
    myplot.p1(data,bottom,top)
    plt.plot(x,zeroline,label='zero line')
    plt.show()
    print("\n\n\nnosie:%f\n\n\n"%noise)
elif flag1 =='3':
    zp=float(input("input zero point\n"))
    noise=float(input("input noise\n"))
    for i in range(length):
        zeroline.append(zp)
else:
    print("wrong input")
    exit()

"""show the result"""
p1=kit4.findpeakg(data,zeroline,noise,bottom,top,r)
p2=kit4.findpeakl(data,zeroline,noise,bottom,top,r)
m=len(p1)
n=len(p2)
if flag1=='1':
    print("\n\n\nzeropoint:%f, noise:%f\n\n\n"%(zp,noise))
elif flag1=='2':
    print("\n\n\nnosie:%f\n\n\n"%noise)
else:
    print("\n\n\nzeropoint:%f, noise:%f\n\n\n"%(zp,noise))
for i in range(m):
    tempbo=p1[i][4]
    tempto=p1[i][5]
    a1=p1[i][2]
    miu1=p1[i][1]
    sigema=p1[i][0]
    mse1=p1[i][3]
    print("the peak%d:\nGuassian: A:%f, miu:%f, sigema:%f, mean square error:%f\n"%(i+1,a1,miu1,sigema, mse1))
    myplot.p2(-a1,miu1,sigema,zeroline,bottom+widthx*tempbo,bottom+widthx*tempto,tempto-tempbo+1)
    myplot.p1(data[tempbo:tempto+1],bottom+widthx*tempbo,bottom+widthx*tempto)
    plt.show()
for i in range(n):
    tempbo=p2[i][4]
    tempto=p2[i][5]
    a2=p2[i][2]
    miu2=p2[i][1]
    mse2=p2[i][3]
    gama2=p2[i][0]
    if gama2>0:
        print("the peak%d:\nLorentzian: A:%f, miu:%f: gama:%f, mean square error:%f\n"%(i+1,a2/gama2,miu2,math.sqrt(gama2),mse2))
    else:
        print("the peak%d:\nLorentzian: A:%f, miu:%f: gama2:%f, mean square error:%f\n"%(i+1,a2/gama2,miu2,gama2,mse2))
    myplot.p3(-a2,miu2,gama2,zeroline,bottom+widthx*tempbo,bottom+widthx*tempto,tempto-tempbo+1)
    myplot.p1(data[tempbo:tempto+1],bottom+widthx*tempbo,bottom+widthx*tempto)
    plt.show()
    
for i in range(m):
    a1=p1[i][2]
    miu1=p1[i][1]
    sigema=p1[i][0]
    myplot.p2(-a1,miu1,sigema,zeroline,bottom,top,len(data))
for i in range(n):
    a2=p2[i][2]
    miu2=p2[i][1]
    gama2=p2[i][0]
    myplot.p3(-a2,miu2,gama2,zeroline,bottom,top,len(data))
myplot.p1(data,bottom,top)
plt.show()
                    
