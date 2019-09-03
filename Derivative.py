# @author Gene Wang
# @version 5 September 2017
# Defines functions that calculate numerical estimates for the
# derivative of and defines the Gaussian function e^-x^2. Also
# contains a function to calculate and print the percent error
# for each numerical value versus the closed form value.
# Gaussian: Takes in an array who's length signifies the domain and
#           returns an 2D array with the values in the domain and their
#           corresponding y values based on the function
# Derivative: Takes in an array who's length signifies the domain and 
#             returns an 2D array with the values in the domain and their
#             corresponding y values based on the closed point derivitave of the function
# Twopoint: Takes in an array with the values of the Gaussian curve and 
#           returns the numerical estimate of the derivative of the function
#           using the twopoint method.
# Threepoint: Takes in an array with the values of the Gaussian curve and 
#             returns the numerical estimate of the derivative of the function
#             using the threepoint method.
# Fivepoint: Takes in an array with the values of the Gaussian curve and 
#             returns the numerical estimate of the derivative of the function
#             using the fivepoint stencil method.
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
# @param dom an array of two integers that specifies the domain of the values that are
#        requested.
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
# domain that are on the curve of the closed form derivitave of the gaussian function e^-x^2.
# @param arr an array of zeroes of the size of all the points in the domain
#        along a given step size.
# @param h the step size between each x value.
# @param dom an array of two integers that specifies the domain of the values that are
#        requested
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
# given points using the twopoint method: taking a midpoint
# between the two points providing the slope of the current
# point and the next point as the derivative estimation
# of the function at the midpoint.
# @param arr the array of (x,y) values that you want to estimate
#        the derivative for.
# @return the array containing all the points (x,y) corresponding to
#         the estimation of the derivitave of the inputed values using
#         the twopoint method.
#
def twopoint(arr):
    
    twopoint = np.zeros([len(arr),2])
    
    for i in range(0,len(arr)-1):
        
        twopoint[i][0] = (arr[i][0]+arr[i+1][0])/2
        
    twopoint[len(arr)-1][0] = arr[len(arr)-1][0]
    
    for j in range(0,len(arr)-1):
        
        pointx = arr[j][0]
        
        pointy = arr[j][1]
        
        nextx = arr[j+1][0]
        
        nexty = arr[j+1][1]
        
        twopoint[j][1] = (nexty-pointy)/(nextx-pointx)
        
    twopoint[len(arr)-1][1] = arr[len(arr)-1][1]
    
    return twopoint

# This function takes in an array of (x,y) values and 
# estimates the derivative of the function formed by the
# given points using the threepoint method: taking the slope between the
# previous and next point at each point as the derivative estimation.
# of the function at the current point.
# @param arr the array of (x,y) values that you want to estimate
#        the derivative for.
# @return the array containing all the points (x,y) corresponding to
#         the estimation of the derivitave of the inputed values using
#         the threepoint method.
#
def threepoint(arr):
    
    threepoint = np.zeros([len(arr),2])
    
    for i in range(0,len(arr)):
        
        threepoint[i][0] = arr[i][0]
        
    for j in range(1,len(arr)-1):
        
        prevx = arr[j-1][0]
        
        prevy = arr[j-1][1]
        
        nextx = arr[j+1][0]
        
        nexty = arr[j+1][1]
        
        threepoint[j][1] = (nexty-prevy)/(nextx-prevx)
        
    threepoint[len(arr)-1][1] = arr[len(arr)-1][1]
    
    threepoint[0][1] = arr[0][1]
    
    return threepoint

# This function takes in an array of (x,y) values and 
# estimates the derivative of the function formed by the
# given points using the fivepoint stencil method.
# @param arr the array of (x,y) values that you want to estimate
#        the derivative for.
# @param h the stepsize between the each given x value.
# @return the array containing all the points (x,y) corresponding to
#         the estimation of the derivitave of the inputed values using
#         the fivepoint stencil method.
#
def fivepoint(arr,h):
    
    fivepoint = np.zeros([len(arr),2])
    
    for i in range(0,len(arr)):
        
        fivepoint[i][0] = arr[i][0]
        
    for j in range(2,len(arr)-2):
        
        fivepoint[j][1] = ((-1*arr[j+2][1])+(8*arr[j+1][1])-(8*arr[j-1][1])+(arr[j-2][1]))/(12*h)
        
    fivepoint[len(arr)-1][1] = arr[len(arr)-1][1]
    
    fivepoint[len(arr)-2][1] = arr[len(arr)-2][1]
    
    fivepoint[0][1] = arr[0][1]
    
    fivepoint[1][1] = arr[1][1]
    
    return fivepoint

# This function calculates the RMS error of the 
# twopoint estimate of an arbitrary function by taking in
# both the estimate and the closed form derivative and
# comparing the values between the two.
# @param actual an array of (x,y) values of the closed form
#        derivative of an arbitrary function.
# @param estimate an array of (x,y) values of the twopoint estimate
#        of the derivative of an arbitrary function.
# @return the RMS error between the twopoint estimate and the closed form derivative.
#
def twopointerror(actual,estimate):
    
    twopointsum = 0
    
    for j in range(0,len(actual)):
        
        twopointsum += abs(actual[j][1]-estimate[j][1]) ** 2 
        
    twopointrms = math.sqrt(twopointsum/len(estimate))
    
    return twopointrms

# This function calculates the RMS error of the 
# threepoint estimate of an arbitrary function by taking in
# both the estimate and the closed form derivative and
# comparing the values between the two.
# @param actual an array of (x,y) values of the closed form
#        derivative of an arbitrary function.
# @param estimate an array of (x,y) values of the threepoint estimate
#        of the derivative of an arbitrary function.
# @return the RMS error between the threepoint estimate and the closed form derivative.
#
def threepointerror(actual,estimate):
    
    threepointsum = 0
    
    for j in range(0,len(actual)):
        
        threepointsum += abs(actual[j][1]-estimate[j][1]) ** 2
        
    threepointrms = math.sqrt(threepointsum/len(estimate))
    
    return threepointrms

# This function calculates the RMS error of the 
# fivepoint estimate of an arbitrary function by taking in
# both the estimate and the closed form derivative and
# comparing the values between the two.
# @param actual an array of (x,y) values of the closed form
#        derivative of an arbitrary function.
# @param estimate an array of (x,y) values of the fivepoint estimate
#        of the derivative of an arbitrary function.
# @return the RMS error between the fivepoint estimate and the closed form derivative.
#
def fivepointerror(actual,estimate):
    
    fivepointsum = 0
    
    for j in range(0,len(actual)):
        
        fivepointsum += abs(actual[j][1]-estimate[j][1]) ** 2
        
    fivepointrms = math.sqrt(fivepointsum/len(estimate))
    
    return fivepointrms

# This function takes in an array of specified length, a step size, and a
# x domain and returns an array that contains points (x,y) on the given
# domain that are on the curve of the sinc function sin(x)/x.
# @param arr an array of zeroes of the size of all the points in the domain
#        along a given step size.
# @param h the step size between each x value.
# @param dom an array of two integers that specifies the domain of the values that are
#        requested
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
# @param dom an array of two integers that specifies the domain of the values that are
#        requested
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
    
    twopointestimategauss = twopoint(gauss)
    
    threepointestimategauss = threepoint(gauss)
    
    fivepointestimategauss = fivepoint(gauss,h)
    
    twopointrmserrorgauss = twopointerror(derivgauss,twopointestimategauss)
    
    threepointrmserrorgauss = threepointerror(derivgauss,threepointestimategauss)
    
    fivepointrmserrorgauss = fivepointerror(derivgauss,fivepointestimategauss)
    
    #Sinc Function
    
    sinc = sincfunc(init,h,dom)
    
    derivsinc = derivativesinc(init2,h,dom)
    
    twopointestimatesinc = twopoint(sinc)
    
    threepointestimatesinc = threepoint(sinc)
    
    fivepointestimatesinc = fivepoint(sinc,h)
    
    twopointrmserrorsinc = twopointerror(derivsinc,twopointestimatesinc)
    
    threepointrmserrorsinc = threepointerror(derivsinc,threepointestimatesinc)
    
    fivepointrmserrorsinc = fivepointerror(derivsinc,fivepointestimatesinc)

    #def main

#Calls the main method.
main()
    