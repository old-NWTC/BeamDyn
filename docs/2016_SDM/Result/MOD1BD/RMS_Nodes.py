import matplotlib.pyplot as plt
import numpy as np

def powerfunc(x):
    return x * x


xvalues = [3,5]

ypower = [0.1,10**(-1.6020599)]

#ypower = powerfunc(xvalues)

ts = [3,4,5,6,9]
rms = [0.254332711217978,0.174299308,0.0310857102456719,0.0262158285550141,0.0266300601840303]

#fig = plt.figure(figsize=(13.5, 6.0), dpi=100)

#plt.title('powerction lin-lin scale')
plt.xlabel('Number of Nodes',fontsize=20,fontname='Helvetica')
plt.ylabel(r'$ \varepsilon_{RMS} (U_{z})$',fontsize=20)
plt.text(3.0,0.04,'$O(n^2)$',fontsize=18)
#plt.xscale('log')
plt.yscale('log')
plt.plot(xvalues,ypower,'r--',ts,rms,'-bo',markerfacecolor='None',markersize=8.0)
plt.grid()

plt.show()