import numpy as np
import math as m
import matplotlib.pyplot as plt



DawnData=np.genfromtxt('DawnFull_tb.txt')*1000
EarthData=np.genfromtxt('Earth_tb.txt')*1000
MarsData=np.genfromtxt('Mars_tb.txt')*1000
VestaData=np.genfromtxt('Vesta_tb.txt')*1000
CeresData=np.genfromtxt('Ceres_tb.txt')*1000

##Dawn JPL coodrds
X=DawnData[:,0]
Y=DawnData[:,1]
Z=DawnData[:,2]

CalcX=[]
CalcY=[]

##Earth JPL coords
EarthX=EarthData[:,0]
EarthY=EarthData[:,1]
EarthZ=EarthData[:,2]

##Mars JPL Coords
MarsX=MarsData[:,0]
MarsY=MarsData[:,1]
MarsZ=MarsData[:,2]


##Vesta JPL Coords
VestaX=VestaData[:,0]
VestaY=VestaData[:,1]
VestaZ=VestaData[:,2]


##Ceres JPL Coords
CeresX=CeresData[:,0]
CeresY=CeresData[:,1]
CeresZ=CeresData[:,2]


radius = m.sqrt(X[0]*X[0]+Y[0]*Y[0])
angle = m.atan(Y[0]/X[0])

V0=[-3.267334027266680E+00,3.345804245220296E+01]*1000
Vmars=[2.239956000086868E+01, 1.417958296152212E+01]




phiV=m.sqrt(V0[0]*V0[0]+V0[1]*V0[1])
print(phiV)
rV=0;

m0=747.1+425+45.6
F0=3*50*0.000

c0=26000
a0=F0/m0
consumption=F0/c0
fuel = 425
k = np.zeros(4)
q = np.zeros(4)
day = 1
day=0;
dday = 1
mas=747.1+fuel+45.6
while day < 510:
    CalcX.append(radius*m.cos(angle))
    CalcY.append(radius*m.sin(angle))
    a = F0/mas
    
    A=(phiV*phiV)/radius
    B=1/(radius*radius)
    rV=rV+dday*25*3600*(A-B+0.866*a)
    A=-(phiV*rV)/radius
    phiV = phiV + dday *25*3600* (A+1*a)
    radius = radius + dday *25*3600* rV  
    angle = angle + dday *25*3600* (phiV/radius)
    mas = mas - dday *25*3600* consumption
   
    
    day = day + dday 
    
     

print(mas)

plt.plot(X,Y,color='black');
plt.plot(CeresX,CeresY,color='blue')
plt.scatter(CeresX[2715],CeresY[2715], color='blue', s=40, marker='o')

plt.plot(EarthX,EarthY,color='green')
plt.scatter(X[0],Y[0], color='green', s=40, marker='o')

plt.plot(VestaX,VestaY,color='cyan')
plt.scatter(VestaX[-1],VestaY[-1], color='cyan', s=40, marker='o')

plt.plot(MarsX,MarsY,color='red')
plt.scatter(MarsX[-1],MarsY[-1], color='red', s=40, marker='o')

plt.plot(CalcX,CalcY,color='orange');

plt.axis('equal')





##plt.grid()
plt.show()