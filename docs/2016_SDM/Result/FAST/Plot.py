import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('Static.dat')   


plt.figure(1)
plt.xlabel('Number of Nodes',fontsize=20,fontname='Helvetica')
plt.ylabel('$U_{flap} (m)$',fontsize=20)
plot1 = plt.plot(data[:,0],data[:,1],'-ro',label='Gauss')
plot2 = plt.plot(data[:,0],data[:,2],'-bo',label='Trapezoidal')
plt.legend() 
plt.grid()


plt.figure(2)
plt.xlabel('Number of Nodes',fontsize=20,fontname='Helvetica')
plt.ylabel('$Blade Mass (kg)$',fontsize=20)
plot1 = plt.plot(data[:,0],data[:,3],'-ro',label='Gauss')
plot2 = plt.plot(data[:,0],data[:,4],'-bo',label='Trapezoidal')
plt.legend() 
plt.grid()

plt.show()