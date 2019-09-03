# @author Gene Wang
# @version 28 November 2017
# Defines a function to model projectile motion using an iterative 
# method and differential equations to solve for the x and y components of velocity and
# x and y positions at each point in time. The differential equations and their derivation
# are defined in the design document titled "Projectile Motion v2.pdf" written for the course 
# ATCS: Numerical Methods for the first semester of the 2017-18 school year at The Harker School 
# by the instructor, Dr. Eric Nelson. This function includes a user interface that allows 
# the user to input values for constants used in the calculations and to select 

# Imports python modules to assist with math and numerical calculations
# Numpy: functions to create arrays of specified sizes
# Math: basic mathematical functions and constants (e, Sin, cos)
import math

import numpy as np

pi = np.arccos(-1) #establishes a numerical constant value for PI

# This function models projectile motion using an iterative 
# method and differential equations to solve for the x and y components of velocity and
# x and y positions at each point in time. The differential equations and their derivation
# are defined in the cited design document. The function writes the x and y position values
# at each point in time in the given .txt file, which can be changed independently in 
# the function. When the function is executed, the user will be prompted to select whether 
# to use the default constant values in the code or to manually input values for these constants.
# The user will also be prompted to select an air density model from the ones listed in the code
# (no air density, constant sea-level air density, adiabatic density model, and isothermal density model).
# The equations used to calculate these values are located in the cited design document.
# The function runs for a given number of iterations, or until the projectile passes below the y = 0.
# @param delt the time interval
# @param iterations the number of iterations used
#

def projectilemotion(delt, iterations):
    
    file = open('ProjectileMotion.txt','w') #change the specific txt file here
    
    test = True
    
    adiabatic = False
    
    isothermal = False
    
    sealevelro = 1.225
    
    while(test):
        
        print("'default' - uses the default values stored")
        
        print("'manual' - prompts the user to input values")
        
        selection = input("Enter one of the modes from above: ")
        
        if(selection == "default"):
            
            v0 = 827 #initial velocity in meters/second
            
            R = .0775 #radius of the projectile in meters
            
            m = 43.2 #mass of the projectile in kg
            
            x0 = 0 #initial velocity in the x direction
            
            y0 = 0 #initial velocity in the y direction
            
            theta0 = 20 #initial firing angle
            
            theta0 = math.radians(theta0) #sets the inital theta to radians
            
            g = 9.80665 #gravitational constant
            
            test = False 
            
        elif(selection == "manual"): #prompts the user to input each value
            
            v0 = float(input("Enter an initial velocity (m/s): "))
            
            R = float(input("Enter a radius for the projectile (m): "))
            
            m = float(input("Enter a mass for the projectile (kg): "))
            
            x0 = float(input("Enter an initial x position: "))
            
            y0 = float(input("Enter an initial y position: "))
            
            theta0 = math.radians(float(input("Enter an initial firing angle (in degrees): ")))
            
            g = float(input("Enter a value for the gravitational constant: "))
            
            test = False
        
        else: #keeps prompting the user until a valid selection is entered
        
            print("\n" + "Please enter a valid choice" + "\n")
        
    test = True
    
    while(test): #prompts the user to choose a density setting
        
        print("\n")
        
        print("'none' - uses a no-air friction model")
        
        print("'constant'- uses a constant sea level air-density approximation")
        
        print("'adiabatic' - uses the adiabatic density profile")
        
        print("'isothermal' - uses the isothermal density profile")
        
        density = input("Choose a model for air density from the selections above: ")
        
        if(density == "none"):
            
            ro = 0 #sets the air density to 0
            
            test = False
        
        elif(density == "constant"):
            
            ro = sealevelro #sets the air density to a constant sea level value
            
            test = False
            
        elif(density == "adiabatic"):
            
            ro0 = sealevelro #sets the sea level value
            
            ro = ro0
            
            a = 6.5e-3 #fits to observational data
            
            alpha = 2.5 #fits to observational data
            
            T0 = 291.15 #temperature at sea level in Kelvin
            
            adiabatic = True
            
            test = False
        
        elif(density == "isothermal"):
            
            ro0 = sealevelro #sets the sea level value
            
            ro = ro0
            
            y = 1e4 
            
            isothermal = True
            
            test = False
        
        else:
            
            print("\n" + "Please enter a valid density selection" + "\n")
            
    print("Executing...")
        
    file.write(" " + str(x0)) #writes the initial x position
    
    file.write("\t" + str(y0) + "\n") #writes the initial y position
    
    vx0 = math.cos(theta0) * v0 #sets initial x velocity
    
    vy0 = math.sin(theta0) * v0 #sets initial y velocity
    
    while(iterations > 0):
        
        vx1 = vx0 - (1/m)*2*pi*R*R*ro*(v0)*(vx0)*(delt) #Equation 16 from the cited design document
        
        x1 = x0 + vx1*delt #Equation 17 from the cited design document
        
        vy1 = vy0 -(((1/m)*2*pi*R*R*ro*(v0)*(vy0)) + g)*(delt) #Equation 18 from the cited design document
        
        y1 = y0 + vy1*delt #Equation 19 from the cited design document
        
        v1 = math.sqrt(vx1*vx1 + vy1*vy1) #Equation 3 from the cited design document
        
        if(adiabatic): #changes the density if using the adiabatic model
            
            ro = ro0 * ((1 - (a*y1)/T0)**alpha) #Equation 20 from the cited design document
            
        if(isothermal): #changes the density if using the isothermal model
            
            ro = ro0 * (math.e**(-y1/y)) #Equation 21 from the cited design document
        
        file.write(" " + str(x1)) #writes the current x position
        
        file.write("\t" + str(y1) + "\n") #writes the current y position
        
        if(y1 <= 0): #terminates the calculations if the projectile passes y = 0 and goes negative
            
            iterations = 0
        
        vx0 = vx1 
        
        x0 = x1
        
        vy0 = vy1
        
        y0 = y1
        
        v0 = v1
        
        iterations -= 1
        
    file.close()
    
    print("Complete!")

#runs the projectile motion function
projectilemotion(.01, 10000)