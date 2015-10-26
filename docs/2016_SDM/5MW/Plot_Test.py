import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys

d1 = pylab.loadtxt("Test19.out")
pylab.plot(d1[:,0],d1[:,1])
pylab.show()