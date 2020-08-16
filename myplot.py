import matplotlib.pyplot as plt
import math
import numpy

def p1(data,bottom,top):
    length=len(data)
    width=(top-bottom)/(length-1)
    x=[]
    for i in range(length):
        x.append(bottom+i*width)
    plt.plot(x,data,label='original')
    plt.legend()

def p2(a,miu,s,zp,bottom,top,n):
    '''plot a Guassian distribution'''
    width=(top-bottom)/(n-1)
    x=[]
    y=[]
    for i in range(n):
        x.append(bottom+i*width)
        ytemp=zp[i]+a*math.exp(-(bottom+width*i-miu)*(bottom+width*i-miu)/(2*s*s))
        y.append(ytemp)
    plt.plot(x,y,label='Guassian')
    plt.legend()

def p3(a,miu,g2,zp,bottom,top,n):
    '''plot a Lorentzian distribution'''
    width=(top-bottom)/(n-1)
    x=[]
    y=[]
    for i in range(n):
        x.append(bottom+i*width)
        ytemp=zp[i]+a/((bottom+width*i-miu)*(bottom+width*i-miu)+g2)
        y.append(ytemp)
    plt.plot(x,y,label='Lorentzian')
    plt.legend()
