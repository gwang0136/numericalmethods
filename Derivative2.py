# @author Gene Wang
# @version 15 September 2017
# Defines functions that calculate numerical estimates for the
# derivative of and defines the Gaussian function e^-x^2 as well as
# sin(x)/x. Also contains functions to calculate and print the 
# percent error for each numerical value versus the closed form value.
# Gaussian: Takes in an array who's length signifies the domain and
#           returns an 2D array with the values in the domain and their
#           corresponding y values based on the function
# Derivative: Takes in an array who's length signifies the domain and 
#             returns an 2D array with the values in the domain and their
#             corresponding y values based on the closed point derivitave of the function
# Parabolic: Takes in an array with the values of the Gaussian curve and 
#            returns the numerical estimate of the derivative of the function
#            using the parabolic fit method.
#

# Imports python modules to assist with math and numerical calculations
# Numpy: functions to create arrays of specified sizes
# Math: basic mathematical functions and constants (e, Sin, cos)
import numpy as np

import math

np.set_printoptions(threshold=np.nan) #prints entire list of values

# This function takes in an array of specified length, a step size, and a
# x domain and returns an array that contains points (x,y) on the given
# domain that are on the curve of the gaussian function e^-x^2.
# @param arr an array of zeroes of the size of all the points in the domain
#        along a given step size.
# @param h the step size between each x value.
# @param dom an array of two integers that specifies the domain of the 
#        values that are requested
# @return the array containing all the points (x,y) on the given domain along using the 
#         given step size of the function e^-x^2
#
def gaussian(arr,h,dom):
    
    i = dom[0]
    
    for l in range(0,len(arr)):
        
        arr[l][0] = i
        
        i += h
        
    for j in range(0,len(arr)):
        
        arr[j][1] = math.e ** (-(arr[j][0]*arr[j][0]))
        
    return arr

# This function takes in an array of specified length, a step size, and a
# x domain and returns an array that contains points (x,y) on the given
# domain that are on the curve of the closed form derivitave of the 
# gaussian function e^-x^2.
# @param arr an array of zeroes of the size of all the points in the domain
#        along a given step size.
# @param h the step size between each x value.
# @param dom an array of two integers that specifies the domain of the 
#        values that are requested
# @return the array containing all the points (x,y) on the given domain along using the 
#         given step size of the closed form derivative of the function e^-x^2
#
def derivativegauss(arr,h,dom):
    
    curval = dom[0]
    
    for l in range(0,len(arr)):
        
        arr[l][0] = curval
    
        curval += h

    for j in range(0,len(arr)):
        
        arr[j][1] = -2*arr[j][0]* (math.e ** (-(arr[j][0]*arr[j][0])))
        
    return arr

# This function takes in an array of (x,y) values and 
# estimates the derivative of the function formed by the
# given points using the parabolic fit method.
# @param arr the array of (x,y) values that you want to estimate
#        the derivative for.
# @return the array containing all the points (x,y) corresponding to
#         the estimation of the derivitave of the inputed values using
#         the parabolic fit method.
#
def parabolic(arr):
    
    parabolic = np.zeros([len(arr),2])
    
    for i in range(0,len(arr)):
        
        parabolic[i][0] = arr[i][0]
    
    for j in range(1,len(arr)-1):
        
        x1 = arr[j-1][0]
        
        x2 = arr[j][0]
        
        x3 = arr[j+1][0]
        
        y1 = arr[j-1][1]
        
        y2 = arr[j][1]
        
        y3 = arr[j+1][1]

        numeratora = x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1)
        
        denominatora = (x1-x2)*(x1-x3)*(x2-x3)
        
        numeratorb = (x3**2)*(y1 - y2) + (x1**2)*(y2 - y3) + (x2**2)*(y3 - y1)
        
        denominatorb = (-x1 + x2) * (x2 - x3) * (-x1 + x3)
        
        a = numeratora/denominatorb
        
        b = numeratorb/denominatorb
        
        parabolic[j][1] = (2 * a * x2) + b
        
        #for j in range(1,len(arr)-1)
        
    parabolic[len(arr)-1][1] = arr[len(arr)-1][1]
    
    parabolic[0][1] = arr[0][1]
    
    return parabolic

# This function calculates the RMS error of the 
# parabolic estimate of an arbitrary function's derivative by taking in
# both the estimate and the closed form derivative and
# comparing the values between the two.
# @param actual an array of (x,y) values of the closed form
#        derivative of an arbitrary function.
# @param estimate an array of (x,y) values of the parabolic estimate
#        of the derivative of an arbitrary function.
# @return the RMS error between the parabolic estimate
#         and the closed form derivative.
#
def parabolicerror(actual,estimate):
    
    parabolicsum = 0
    
    for j in range(0,len(actual)):
        
        parabolicsum += abs(actual[j][1]-estimate[j][1]) ** 2
        
    parabolicrms = math.sqrt(parabolicsum/len(estimate))
    
    return parabolicrms

# This function takes in an array of specified length, a step size, and a
# x domain and returns an array that contains points (x,y) on the given
# domain that are on the curve of the sinc function sin(x)/x.
# @param arr an array of zeroes of the size of all the points in the domain
#        along a given step size.
# @param h the step size between each x value.
# @param dom an array of two integers that specifies the domain of the 
#        values that are requested
# @return the array containing all the points (x,y) on the given domain along using the 
#         given step size of the function sin(x)/x.
#
def sincfunc(arr,h,dom):
    
    currval = dom[0]
    
    for l in range(0,len(arr)):
        
        arr[l][0] = currval
        
        currval += h
        
    for j in range(0,len(arr)):
        
        arr[j][1] = math.sin(arr[j][0])/arr[j][0]
        
    return arr

# This function takes in an array of specified length, a step size, and a
# x domain and returns an array that contains points (x,y) on the given
# domain that are on the curve of the closed form derivitave of the 
# sinc function sin(x)/x.
# @param arr an array of zeroes of the size of all the points in the domain
#        along a given step size.
# @param h the step size between each x value.
# @param dom an array of two integers that specifies the domain of the 
#        values that are requested
# @return the array containing all the points (x,y) on the given domain along using the 
#         given step size of the closed form derivative of the function sin(x)/x.
#
def derivativesinc(arr,h,dom):
    
    i = dom[0]
    
    for l in range(0,len(arr)):
        
        arr[l][0] = i
        
        i += h
        
    for j in range(0,len(arr)):
        
        currx = arr[j][0]
        
        arr[j][1] = (math.cos(currx)/currx)-(math.sin(currx)/(currx ** 2))
        
    return arr

# This function utilizes all of the previously defined functions to store the 
# functions and their derivative estimations in arrays. This function also establishes
# the domain, step size, and zero arrays utilized in the estimation functions.
#
def main():
    
    dom = [-100,100] #requested domain
    
    h = 0.1 #requested step-size
    
    n = (dom[1]-dom[0])*(int(1/h)) #number of values along domain using given step size.
    
    init = np.zeros([n+1,2]) #array of zeroes (x,y) of specified length
    
    init2 = np.zeros([n+1,2]) #array of zeroes (x,y) of specified length
    
    #Gaussian Curve
    
    gauss = gaussian(init,h,dom)
    
    derivgauss = derivativegauss(init2,h,dom)
    
    parabolicgauss = parabolic(gauss)
    
    parabolicerrorgauss = parabolicerror(derivgauss,parabolicgauss)
    
    #Sinc Function
    
    sinc = sincfunc(init,h,dom)
    
    derivsinc = derivativesinc(init2,h,dom)

    parabolicsinc = parabolic(sinc)
    
    parabolicerrorsinc = parabolicerror(derivsinc,parabolicsinc)

#Calls the main function
main()