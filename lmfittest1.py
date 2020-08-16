from numpy import exp,sqrt
from lmfit import Model
import random
import matplotlib.pyplot as plt

def twogaussian(x,a1,miu1,sigma1,a2,miu2,sigma2,a3,miu3,sigma3):
    temp=a1*exp(-(x-miu1)*(x-miu1)/(2*sigma1*sigma1))+a2*exp(-(x-miu2)*(x-miu2)/(2*sigma2*sigma2))+a3*exp(-(x-miu3)*(x-miu3)/(2*sigma3*sigma3))
    return temp
x=[]
y=[]
for i in range(-50,50):
    x.append(i)
    y.append(exp(-(i+25)*(i+25)/8)+2*exp(-(i-3)*(i-3)/18)+3*exp(-(i-24)*(i-24)/8)+random.randint(0,9)/40-0.125)
mod=Model(twogaussian)
pars=mod.make_params(a1=1,miu1=-20,sigma1=1,a2=1,miu2=0,sigma2=1,a3=1,miu3=20,sigma3=1)
result=mod.fit(y,pars,x=x)
print(result.fit_report())
plt.plot(x,y,'bo')
plt.plot(x,result.best_fit,'r-')
plt.show()
