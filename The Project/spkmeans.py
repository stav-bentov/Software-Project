import numpy as np
import pandas as pd
import sys
import mykmeanssp

def call_jacobi(K, filename, max_iter):
    epsilon = 0.00001 # epsilon = 1.0 x 10^-5
    with open('filename', 'r') as f:#todo check of course
        A = [[int(num) for num in line.split(',')] for line in f] # maybe todo  read it from C program??

    res = mykmeanssp.fit_jacobi(len(A),K,max_iter,A,epsilon)
    #todo sorting of eigenvalues , vectors increasingly
    printOutput_jacobi(res)


def printOutput_jacobi(returnFromFit): #todo check if jacobi also requires commas!!!!
    output_res=""
    no_of_rows = len(returnFromFit)

    for i in range(no_of_rows):
        len_of_row_i = len(returnFromFit[i])
        for j in range(len_of_row_i):
            output_res+=str('%.4f' %(returnFromFit[i][j])) # returnFromFit[1] = centroid got from kmenassp.c
            if (j!=len_of_row_i-1):
                output_res += ","
        output_res += "\n"
    print(output_res)


def call_spk(K, data_points_array): # same as in ex2
    pass


def call_wam(K, data_points_array):
    pass


def call_ddg(K, data_points_array):
    pass


def call_lnorm(K, data_points_array):
    pass







def checkInput(given_input, argLen): #todo it in main func (check if k is an int + >=0 and that argv[2] is in enum list)
    if (argLen != 4):
        print("Invalid Input!")
        sys.exit()
    #checking if K is a non-negative integer
    isValid = given_input[1].isnumeric()
    if(int(given_input[1]) < 0 or isValid == False):
        print("Invalid Input!")
        sys.exit()

def findK(): # heuristic method
    pass


def getNDataPoints(filename):
    try:
        return pd.read_csv(filename, header=None) # Put data from file to a data array
    except:
        print("An Error Has Occurred")
        sys.exit()


def main(argv):

    '''=========================Checking input========================='''
    argLen = len(argv)
    checkInput(argv,argLen)

    '''=========================Set arguments========================='''
    max_iter = 100 # Default value for jacobi todo
    K = (int)(argv[1])
    if K == 0:
        K = findK() # heuristic method in case that K is 0
    goal = argv[2] #todo enum!

    if goal == "jacobi":
        call_jacobi(K, argv[3]) #todo check whats the matrix to pass : Lnorm or just the matrix we get from user?
    else:
        data_points_array = getNDataPoints(argv[3])
        N = len(data_points_array)
        if (K > N): #todo also check in jacobi
            print("Invalid Input!")
            sys.exit()

        if goal == "spk":
            call_spk(K, data_points_array)
        if goal == "wam":
            call_wam(K, data_points_array)
        if goal == "ddg":
            call_ddg(K, data_points_array)
        if goal == "lnorm":
            call_lnorm(K, data_points_array)
        else: #goal is invalid
            print("Invalid Input!")
            sys.exit()

        #todo dimension=data_points_array.shape[1]-1



if __name__ == '__main__':
    main(sys.argv)
