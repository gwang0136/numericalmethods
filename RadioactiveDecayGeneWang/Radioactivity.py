# @author Gene Wang
# @version 13 November 2017
# Defines a method to model radioactive decay
# using an iterative method and differential
# equations to solve for the particle count
# at each point.
#

# Imports python modules to assist with math and numerical calculations
# Numpy: functions to create arrays of specified sizes
# Math: basic mathematical functions and constants (e, Sin, cos)
import math

import numpy as np

# Function models radioactive decay using an iterative 
# method and the differential equation dN(t) = (-N(t)/tao)dt 
# to solve for the particle count at each point in a given 
# number of iterations, an initial particle count value, 
# a given delta t value, and a tao value that represents
# the time between decays. The function writes each set 
# of values (time, particle count) in a given .txt file 
# which can be changed in the method.
# @param init the initial particle count
# @param delt the delta t value used in the calculation.
# @param tao the time between decays.
# @param iterations the number of iterations to be run
#
def radio(init, delt, tao, iterations):
    
    file = open('Radioactivity.txt','w') #change the specific txt file here
    
    N0 = init
    
    file.write(" " + str(0)) #writes the initial t value (0)
    
    file.write("\t" + str(N0) + "\n") #writes the initial partical count
    
    N1 = 0 #"next" value initialized as zero
    
    t = 0 #initial time set as zero
    
    for x in range(1,iterations):
        
        N1 = N0 * (1 - delt/tao)
        
        file.write(" " + str(y)) #writes the current t value
        
        file.write("\t" + str(N1) + "\n") #writes the current particle count
        
        N0 = N1
        
        t += delt
        
    file.close() #closes the file

#def radio(init, delt, tao, iterations):

#runs the method
#change the initial particle count, delta t value,
#tao value, and iterations here in the parameters.
radio(6.022e23, .001, 29, 100000)