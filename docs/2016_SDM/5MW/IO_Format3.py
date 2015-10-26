import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys

d1 = pylab.loadtxt("5MW_Blade_Property.txt",skiprows=2)
d2 = np.zeros(11)
for i in range(0,11):
    d1[i,3]=d1[i,3]*1.04536
    d2[i] = d1[i,9]+d1[i,10]

f=open('5MW_Blade_F3.inp','w')
for i in range(0,11):
    (f.write('%14.5E' % d1[i,0]),
    f.write('%14.5E' % d1[i,3]),f.write('%14.5E' % d1[i,10]),f.write('%14.5E' % d1[i,9]),f.write('%14.5E' % 0),f.write('%14.5E' % d1[i,14]),f.write('%14.5E' % 0),
    f.write('%14.5E' % d1[i,7]),f.write('%14.5E' % d1[i,5]),f.write('%14.5E' % d1[i,4]),f.write('%14.5E' % 0),f.write('%14.5E' % 0),f.write('%14.5E' % 0),
    f.write('%14.5E' % d1[i,6]),f.write('%14.5E' % d1[i,7]),f.write('%14.5E' % d1[i,7]),f.write('%14.5E' % 0),f.write('%14.5E' % 0),f.write('%14.5E' % 0),
    f.write('\n')
    )
f.close()