# @author Gene Wang
# @version 23 October 2017
# Defines a method that finds a functional fit for a given set of data
# points from a text file by minimizing the RMS error using the method
# of steepest descent. Also contains helper methods that find partial
# derivatives, calculate the error and partial error, and defines the
# arbitrary function.
# Function: takes in an array that contains an x value along with a 
#           predetermined number of q values and returns the value
#           of a general function at the point x with the given q values.
# Partial Deriv: estimates the partial derivative at the point x
#                with respect to the change in a given q value 
#                using the five-point stencil method.
# Error: calculates the RMS error of the function using the given
#        q values.
# Partial Error: calculates the RMS error of the derivative of the function
#                using the given q values and with respect to the change in
#                the specified q value.
# Least Squares: finds a functional fit for a given set of data points by
#                minimizing the RMS error using the method of steepest descent
#

# Imports python modules to assist with math and numerical calculations
# Numpy: functions to create arrays of specified sizes
# Math: basic mathematical functions and constants (e, Sin, cos)
import numpy as np

import math

np.set_printoptions(threshold=np.nan) #prints entire list of values

# This function takes in an array containing
# an x value along with various q values
# and returns the value of a general function
# at the point x with the given q values.
# @param arr the array with the x and q values.
# @return the value of the selected function
#         at the point x with the given q values.
#
def function(arr):
    
    x = arr[0]
    
    q0 = arr[1]
    
    q1 = arr[2]
    
    q2 = arr[3]
    
    q3 = arr[4]
    
    #q4 = arr[5]
    
    #q5 = arr[6]
    
    #q6 = arr[7]
    
    y = (q0) * (math.e **( - ((x-q1)**2) / ((q2)**2)) ) + q3
    
    #y = q0*x + q1
    
    #y = q0*((x-q1)**2) + q2
    
    #y = q0 + q1*(math.e **(-((q2-x)**2)/(q3)**2)) + q4*(math.e **(-((q5-x)**2)/(q6)**2))
    
    return y

# This function estimates the partial derivative 
# at the point x with respect to the change in 
# a given q value using the five-point stencil method.
# @param vals the array with the x and q values.
# @param h the stepsize.
# @param index the index of the requested 
#        q value in the array vals.
# @return the estimated value of the partial derivative
#         at the point x.
#
def partialderiv(vals, h, index):
    
    arr1 = np.zeros(len(vals))
    
    arr2 = np.zeros(len(vals))
    
    arr3 = np.zeros(len(vals))
    
    arr4 = np.zeros(len(vals))
    
    for x in range(len(vals)):
        
        arr1[x] = vals[x]
        
        arr2[x] = vals[x]
        
        arr3[x] = vals[x]
        
        arr4[x] = vals[x]
        
    arr1[index] += 2*h
    
    arr2[index] += h
    
    arr3[index] -= h
    
    arr4[index] -= 2*h
    
    partial = ((-1*function(arr1))+(8*function(arr2))-(8*function(arr3))+(function(arr4)))/((12*h))

    return partial
    
# This function calculates the RMS error of the 
# function using the given q values.
# @param realvals an array containing the 
#        given data points in the form [x,y].
# @param vals an array with the x and q values.
# @return the RMS error using the given q values.
#
def error(realvals,vals):
    
    totalsumsquare = 0
    
    for x in range(len(realvals)):
        
        vals[0] = x
        
        totalsumsquare += (realvals[x][1]-function(vals))**2
        
    return (.5 * totalsumsquare)

# This function calculates the RMS error of the derivative 
# of the function using the given q values and with respect 
# to the change in the specified q value.
# @param realvals an array containing the
#        given data points in the form [x,y].
# @param vals an array with the x and q values.
# @param index the index of the requested
#        q value in the array vals.
# @param h the stepsize.
# @return the RMS partial error of the derivative function.
#
def partialerror(realvals,vals,index,h):
    
    totalsumpartial = 0
    
    for x in range(len(realvals)):
        
        vals[0] = x

        totalsumpartial += (realvals[x][1]-function(vals))*(partialderiv(vals,h,index))

    return totalsumpartial

# This function finds a functional fit for a given set 
# of data points by minimizing the RMS error using the 
# method of steepest descent.
# @param realvals an array containing the 
#        given data points in the form [x,y].
# @param initialqvals an array containing the 
#        user designated starting q values.
# @param iterations the number of iterations
#        the code runs before terminating.
# @param h the stepsize.
# @param tolerance the tolerance value that
#        the user are willing to accept.
# @return an array with the q values of the
#         function that has the smallest RMS 
#         error through the number of iterations
#         or a value of error that is less than
#         the designated tolerance value.
# 
def leastsquares(realvals, initialqvals, iterations, h, tolerance):
    
    lam = 1
    
    numiter = iterations
    
    currx = realvals[0][0]
    
    totalvals = np.zeros(len(initialqvals)+1)
    
    totalvals[0] = currx
    
    newtotalvals = np.zeros(len(totalvals))

    finalqvals = np.zeros(len(initialqvals))

    for x in range(1,len(totalvals)):
        
        totalvals[x] = initialqvals[x-1]
        
    currerror = error(realvals, totalvals)
    
    print(currerror)

    for i in range(iterations):
        
        if(currerror <= tolerance):
            
            for x in range(1,len(totalvals)):
                
                finalqvals[x-1] = totalvals[x]
            
            print("q-values found within an error tolerance of: " + str(tolerance))
            
            return finalqvals
        
        for j in range(1,len(totalvals)):
            
            delq = lam * partialerror(realvals,totalvals,j,h)
            
            newtotalvals[j] = totalvals[j] + delq

        preverror = currerror
        
        currerror = error(realvals, newtotalvals)
        
        if((abs(currerror - preverror) <= tolerance) and test):
            
            for x in range(1,len(newtotalvals)):
                
                finalqvals[x-1] = newtotalvals[x]
            
            print("Error function changing by a value less than or equal to: " + str(tolerance))
            
            print("(Not changing sufficiently quick enough)")
            
            print("Error value: " + str(currerror))
            
            print("q-values: ")
            
            return finalqvals
        
        if(currerror > preverror):
            
            currerror = preverror
            
            test = False
            
            lam /= 2
            
        if(currerror < preverror):
            
            for x in range(1,len(newtotalvals)):
                
                totalvals[x] = newtotalvals[x]
                
                test = True

    for x in range(1,len(totalvals)):
                
        finalqvals[x-1] = totalvals[x]
    
    print("q-values found after " + str(numiter) + " iterations")
    
    print("Error value: " + str(currerror))
    
    return finalqvals
        
# main method
# initializes the number and value of
# starting q values

n = 4 #number of q values

qvals = np.zeros(n)

#starting points for the q values

qvals[0] = 300000

qvals[1] = 10

qvals[2] = 10

qvals[3] = 500

#qvals[4] = 

#qvals[5] = 

#qvals[6] = 

# Opens the text file with the values and reads the 
# data into an array as [x,y] data points.
with open('GeigerHisto.txt') as data:
    
    counter = 0
    
    for line in data:
        
        counter += 1 
    
readvals = np.zeros([counter,2])
     
with open('GeigerHisto.txt') as data:
    
    index = 0
    
    for line in data:
        
        numbers = line.split()
        
        readvals[index][0] = float(numbers[0])
        
        readvals[index][1] = float(numbers[1])
        
        index +=1

# Calls the least squares method with a user selected
# number of iterations, step size, and tolerance. 
leastsquares(readvals, qvals, 10000, 1e-7, 1e-100)

