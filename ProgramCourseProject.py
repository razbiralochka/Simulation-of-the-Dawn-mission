import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


DawnData=np.genfromtxt('DawnFull_tb.txt')
EarthData=np.genfromtxt('Earth_tb.txt')
MarsData=np.genfromtxt('Mars_tb.txt')
VestaData=np.genfromtxt('Vesta_tb.txt')
CeresData=np.genfromtxt('Ceres_tb.txt')





X=DawnData[:,0]
Y=DawnData[:,1]
Z=DawnData[:,2]

EarthX=EarthData[:,0]
EarthY=EarthData[:,1]
EarthZ=EarthData[:,2]

MarsX=MarsData[:,0]
MarsY=MarsData[:,1]
MarsZ=MarsData[:,2]

VestaX=VestaData[:,0]
VestaY=VestaData[:,1]
VestaZ=VestaData[:,2]

CeresX=CeresData[:,0]
CeresY=CeresData[:,1]
CeresZ=CeresData[:,2]


plt.plot(CeresX,CeresY)
plt.plot(EarthX,EarthY)
plt.plot(VestaX,VestaY)
plt.plot(MarsX,MarsY)
plt.plot(X,Y)

plt.axis('equal')
plt.grid()
plt.show()