# @author Gene Wang
# @version 2 January 2018
# Defines classes and functions to model gravitational attraction between an N number of bodies 
# interacting with each other, or the "N-Body" Problem. The differential equations used for these 
# interactions and their derivations are defined in the design document titled "N-Body v1" written 
# for the course  ATCS: Numerical Methods for the first semester of the 2017-18 school year at 
# The Harker School by its instructor, Dr. Eric Nelson.
#

# Imports python modules to assist with math and numerical calculations
# as well as plots
# Numpy: functions to create arrays of specified sizes
# Math: basic mathematical functions and constants (e, Sin, cos)
# Matplotlib: plotting software in python
# Axes3D: allows for 3D plots of the objects.
import math

import numpy as np

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

pi = np.arccos(-1) #establishes a numerical constant value for PI

G = 6.674e-11 #m^3⋅kg−1⋅s−2 #establishes a numerical constant value for G

# This class defines the methods and initializes the variables for a Star object.
# An object of the class Star has variables that describe its position in space
# (as a vector), its mass, and its velocity (as a vector).
#
class Star:
    
    # Constructor initializes all of the instance variables upon
    # instantiation of the object.
    # @param self the bound version of this object
    # @param position a Vector containing the position of the Star (in parsecs)
    # @param mass the mass of the Star (in solar masses)
    # @param velocity a Vector containing the velocity of the Star (in km/s) 
    def __init__(self, position, mass, velocity):
        
        selfx = position.getx()#*3.086e+16 #changes the parsec value to meters
        
        selfy = position.gety()#*3.086e+16 #changes the parsec value to meters
        
        selfz = position.getz()#*3.086e+16 #changes the parsec value to meters
        
        self.pos = Vector(selfx,selfy,selfz)
        
        self.mass = mass * 1.98e30 #changes the solar mass value to kg
        
        self.v = velocity
        
    # This function returns the current position of the Star object
    # @param self the bound version of this object
    # @return a Vector with the current position of the object
    def getpos(self):
            
        return self.pos
    
    # This function sets the current position of the Star object
    # to the Vector parameter given
    # @param self the bound version of this object
    # @param vector the vector position that you want to set the position of this object to
    def setpos(self,vector):
        
        self.pos = vector
        
    # This function returns the current mass of the Star object
    # @param self the bound version of this object
    # @return the current mass of the object
    def getmass(self):
            
        return self.mass
    
    # This function sets the current mass of the Star object
    # to the parameter given
    # @param self the bound version of this object
    # @param val the mass value that you want to set the mass of this object to
    def setmass(self,val):
        
        self.mass = val
        
    # This function returns the current velocity of the Star object
    # @param self the bound version of this object
    # @return the current velocity of the object as a vector
    def getvelocity(self):
            
        return self.v
    
    # This function sets the current velocity of the Star object
    # to the Vector parameter given
    # @param self the bound version of this object
    # @param vector the vector velocity that you want to set the velocity of this object to
    def setvelocity(self,vector):
        
        self.v = vector

# This class defines the methods and initializes the variables for a Vector object.
# An object of the class Vector has variables that describe three different components
# for a three dimensional vector, and contains methods to change and acquire these values.
#
class Vector:
    
    # Constructor initializes all of the instance variables upon
    # instantiation of the object.
    # @param self the bound version of this object
    # @param xposition the x component of the Vector
    # @param yposition the y component of the Vector
    # @param zposition the z component of the Vector
    def __init__(self, xposition, yposition, zposition):
        
        self.x = xposition
        
        self.y = yposition
        
        self.z = zposition
        
        
    # This function returns the x component value of the vector
    # @param self the bound version of this object
    def getx(self):
        
        return self.x
    
    # This function returns the y component value of the vector
    # @param self the bound version of this object
    def gety(self):
        
        return self.y
    
    # This function returns the z component value of the vector
    # @param self the bound version of this object
    def getz(self):
        
        return self.z

# This function takes in two separate vectors and adds them together using vector addition.
# @param vector1 the first vector to add
# @param vector2 the second vector to add
# @return a new vector that represents the sum of the 2 given vectors.
#
def addition(vector1,vector2):
    
    newvector = Vector(vector1.getx()+vector2.getx(), vector1.gety()+vector2.gety(), vector1.getz()+vector2.getz())
    
    return newvector

# This function takes in two separate vectors and subtracts the second vector 
# from the first using vector subtraction.
# @param vector1 the first vector
# @param vector2 the second vector to subtract from the first
# @return a new vector that represents the difference of the 2 given vectors.
#
def difference(vector1,vector2):
    
    newvector = Vector(vector1.getx()-vector2.getx(), vector1.gety()-vector2.gety(), vector1.getz()-vector2.getz())
    
    return newvector

# This function takes in a vector and calculates the magnitude of the vector
# @param vector the vector that you want the magnitude of
# @return the magnitude of the given vector.
#
def magnitude(vector):

    return (math.sqrt(vector.getx()*vector.getx()+vector.gety()*vector.gety()+vector.getz()*vector.getz()))

# This function takes in a list of Star objects, a time period, and a 
# change in time, and allows the Star objects to evolve based on the interactions
# due to gravitational attraction until the given time is up.
# The gravitational forces are modeled by the equations 11 and 12 in the design document.
# @param stars the array with stars 
# @param time the total amount of time in seconds
# @param delt the change in time
#
#
def evolve(stars, time, delt):
    
    counter = 0
    
    while(counter <= time):
    
        for star1 in stars:
        
            totalsumdelx = 0
            
            totalsumdely = 0
            
            totalsumdelz = 0
            
            pos1 = star1.getpos()
            
            pos1x = pos1.getx()
            
            pos1y = pos1.gety()
            
            pos1z = pos1.getz()
        
            for star2 in stars:
            
                if(star1 != star2):
                    
                    pos2 = star2.getpos()
                    
                    pos2x = pos2.getx()
                        
                    pos2y = pos2.gety()
                        
                    pos2z = pos2.getz()
                    
                    if(magnitude(difference(pos1,pos2)) == 0):
                        
                        newvelocity1 = Vector(star1.getmass*pos1x,star1.getmass*pos1y,star1.getmass*pos1z)
                        
                        newvelocity2 = Vector(star2.getmass*pos2x,star2.getmass*pos2y,star2.getmass*pos2z)
                        
                        newv = addition(newvelocity1,newvelocity2)
                        
                        summass = star1.getmass + star2.getmass
                        
                        newv1 = Vector(newv.getx()/summass, newv.gety()/summass, newv.getz()/summass)
                        
                        star1.setvelocity(newv1)
                        
                        star1.setmass(star1.getmass()+star2.getmass())
                        
                        star2.setvelocity(Vector(0,0,0))
                        
                        star2.setmass(0)
                    
                    else:
                        
                        totalsumdelx += (star2.getmass() * G * (pos1x-pos2x) * delt)/(magnitude(difference(pos1, pos2)) ** 3)
                    
                        totalsumdely += (star2.getmass() * G * (pos1y-pos2y) * delt)/(magnitude(difference(pos1, pos2)) ** 3)
                    
                        totalsumdelz += (star2.getmass() * G * (pos1z-pos2z) * delt)/(magnitude(difference(pos1, pos2)) ** 3)
                    
                
            vx0 = star1.getvelocity().getx()
                    
            vy0 = star1.getvelocity().gety()
                    
            vz0 = star1.getvelocity().getz()
            
            #print(totalsumdelx,totalsumdely,totalsumdelz)
                    
            v1 = Vector(vx0 - totalsumdelx, vy0 - totalsumdely, vz0 - totalsumdelz)
                    
            star1.setvelocity(v1)
                    
        for star in stars:
            
            x0 = star.getpos().getx()
            
            y0 = star.getpos().gety()
            
            z0 = star.getpos().getz()
            
            x1 = x0 + star.getvelocity().getx() * delt
            
            y1 = y0 + star.getvelocity().gety() * delt
            
            z1 = z0 + star.getvelocity().getz() * delt
            
            #print(x1,y1,z1)
            
            r1 = Vector(x1,y1,z1)
            
            star.setpos(r1)
            
            #print(star.getpos().getx())
            
        counter += delt
