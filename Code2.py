import numpy as np
import math as m
import matplotlib.pyplot as plt
from Equs import calcs_class


##Ephemeris/Эфимириды
DawnData=np.genfromtxt('DawnFull_tb.txt')*1000
EarthData=np.genfromtxt('Earth_tb.txt')*1000
MarsData=np.genfromtxt('Mars_tb.txt')*1000

EarthX=EarthData[:,0]
EarthY=EarthData[:,1]

MarsX=MarsData[:,0]
MarsY=MarsData[:,1]

k = np.sqrt(EarthX[0]**2 + EarthY[0]**2)

EarthX /= k
EarthY /= k


MarsX /= k
MarsY /= k


plt.scatter(0,0, color='yellow', s=100, marker='o')
plt.plot(EarthX,EarthY,color='blue',label='Орбита Земли')
plt.scatter(EarthX[0],EarthY[0], color='blue', s=40, marker='o')
plt.plot(MarsX,MarsY,color='red',label='Орбита Марса')
plt.scatter(MarsX[-1],MarsY[-1], color='red', s=40, marker='o')




calcs = calcs_class(MarsX[-1],MarsY[-1])

input = np.random.uniform(-1,1,4)


inp = calcs.fit(input)

out = calcs.calc(inp)
X,Y, u1,u2 = calcs.get_points()

T = [510*i/10000 for i,item in enumerate(X)]



plt.plot(X,Y,color='black',label='Траектория КА')

plt.axis('equal')
plt.legend(loc = 2)
plt.show()

plt.plot(T,u1,label='Угол')
plt.legend(loc = 2)
plt.show()

plt.plot(T,u2,label='Ускорение')
plt.legend(loc = 2)
plt.show()