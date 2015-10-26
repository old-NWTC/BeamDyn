import matplotlib.pyplot as plt
import numpy as np

dataA = np.loadtxt('ModBD21OCT15.dat')   #ANSYS Data
dataAA= np.loadtxt('ModBD22OCT15.dat')   #ANSYS Data with Acceleration
dataBT= np.loadtxt('QiDisp_Tip_Mod1_BD.out') #BD Tip Displacement
dataBR= np.loadtxt('Qi_Mod1Disp.out') #BD Root (mass) Displacement
dataBA= np.loadtxt('Qi_Mod1Vel.out')  #BD Velocities and acceleartions

plt.figure(1)
plt.xlabel('Time (s)',fontsize=20,fontname='Helvetica')
plt.ylabel('$U_z (m)$',fontsize=20)
plot1 = plt.plot(dataA[:,0],dataA[:,1],'-r',label='ANSYS')
plot2 = plt.plot(dataBR[:,0],dataBR[:,1],'-b',label='BeamDyn')
plt.legend() 

plt.grid()


plt.figure(2)
plt.xlabel('Time (s)',fontsize=20,fontname='Helvetica')
plt.ylabel('$U_z (m)$',fontsize=20)
plot1 = plt.plot(dataA[:,0],dataA[:,2],'-r',label='ANSYS')
plot2 = plt.plot(dataBT[:,0],dataBT[:,3],'-b',label='BeamDyn')
plt.legend() 

plt.grid()

plt.figure(3)
plt.xlabel('Time (s)',fontsize=20,fontname='Helvetica')
plt.ylabel('$A_z (m/s^2)$',fontsize=20)
plot1 = plt.plot(dataAA[:,0],dataAA[:,3],'-r',label='ANSYS')
plot2 = plt.plot(dataBA[:,0],dataBA[:,4],'-b',label='BeamDyn')
plt.legend() 
plt.grid()

plt.show()