# @author Gene Wang
# @version 15 November 2017
# Defines a method to model simple pendulum motion using an iterative 
# method and differential equations to solve for the angular displacement
# and angular velocity at each point in time. Please refer to the 
# following design document for the differential equations and their derivation:
# https://athena.harker.org/course/1152351352/materials/gp/1359313135
#

# Imports python modules to assist with math and numerical calculations
# Numpy: functions to create arrays of specified sizes
# Math: basic mathematical functions and constants (e, Sin, cos)
import math

import numpy as np

# This function models simple pendulum motion using an iterative 
# method and differential equations to solve for the angular displacement
# and angular velocity at each point in time. The equations used for the 
# change in angular velocity and angular displacement are equations 4 and 5
# in the above design document respectively. The function writes each set 
# of values (time, theta(in radians)) in a given .txt file which can be 
# changed in the method.
# @param initw the initial angular velocity
# @param inittheta the initial angular displacement in degrees
# @param g the gravitational constant
# @param l the length of the pendulum
# @param delt the time interval
# @param iterations the number of iterations used
#
def pendulummotion(initw, inittheta, g, l, delt, iterations):
    
    file = open('HarmonicMotion.txt','w') #change the specific txt file here
    
    w0 = initw
    
    t0 = math.radians(inittheta) #changes the inputted inital angle to radians
    
    file.write(" " + str(0)) #writes the initial t value (0)
    
    file.write("\t" + str(t0) + "\n") #writes the initial angular displacement
    
    w1 = 0
    
    t1 = 0
    
    t = delt
    
    for x in range(iterations):
        
        w1 = w0 - (((g/l) * math.sin(t0)) * delt) 
        
        t1 = t0 + (w1 * delt)
        
        file.write(" " + str(t)) #writes the current t value
        
        file.write("\t" + str(t1) + "\n") #writes the current angular displacement
    
        w0 = w1
        
        t0 = t1
        
        t += delt
        
    file.close() #closes the file

#defines start values used in the calculation of the pendulum model
initialw = 0 #initial angular velocity (set at 0)
initialtheta = 60 #initial angular displacement in degrees
g = 9.80665 #value of gravitational constant
length = 10 #length of the pendulum
deltat = .01 #time interval
numiter = 10000 #number of iterations

#runs the pendulum model function using the above initial start values
pendulummotion(initialw, initialtheta, g, length, deltat, numiter)