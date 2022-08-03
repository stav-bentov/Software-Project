import math

import numpy as np
import pandas as pd
import sys
import mykmeanssp

def print_output(returnFromFit,N,goal):  # todo check if jacobi also requires commas and also print as columns or rows!
    output_res = ""
    no_of_rows=N
    if(goal=="jacobi"):
        no_of_rows+=1

    for i in range(no_of_rows):
        for j in range(N):
            output_res += str('%.4f' %(returnFromFit[i][j]))  # todo check with Stav about %
            if (j != N-1):
                output_res += ","
        output_res += "\n"
    print(output_res)

#todo return list
# '''=========================get data from file========================='''
def getNDataPoints(filename, goal):
    try:
        if goal == "jacobi":
            with open('argv[3]', 'r') as f:  # todo check of course
                return [[float(num) for num in line.split(',')] for line in f]
        else: # goal is not jacobi
            return pd.read_csv(filename, header=None)  # Put data from file to a data array
    except:
        print("An Error Has Occurred")
        sys.exit()

def printOutput(returnFromFit,dimension,K):
    output_res=""
    # Printing first row - indices of K randomly chosen centroids (initial centroids)
    for i in range(K):
        output_res += str(int(returnFromFit[0][i])) # returnFromFit[0] = randomly chosen centroids indices
        if (i != K-1):
            output_res += ","
    output_res += "\n"

    # Printing the K calculated final centroids
    for i in range(K):
        for j in range(dimension):
            output_res+=str('%.4f' %(returnFromFit[1][i][j])) # returnFromFit[1] = centroid got from kmenassp.c
            if (j!=dimension-1):
                output_res += ","
        output_res += "\n"
    print(output_res)

# '''=========================checking input========================='''
def checkInput(given_input, argLen): #todo it in main func (check if k is an int + >=0 and that argv[2] is in enum list)
    if (argLen != 4):
        print("Invalid Input!")
        sys.exit()
    #checking if K is a non-negative integer
    isValid = given_input[1].isnumeric()
    if(int(given_input[1]) < 0 or isValid == False):
        print("Invalid Input!")
        sys.exit()

# ''' spk from ex2'''#
def kMeans_init(K, mergeDf):

    Centroids_array = []  # Saves the centroids u1, ... , uK

    mergeDf = mergeDf.sort_values(by=[0])
    data_points_array = mergeDf.to_numpy()  # Saves the data points

    N = len(data_points_array)
    if (K > N):
        print("Invalid Input!")
        sys.exit()

    D_array = np.array([0.0 for i in range(N+1)])
    Pr_array = np.array([0.0 for i in range(N)]) # Pr_array[i] = probability of
    Index_array = np.array([i for i in range(N)])
    np.random.seed(0)
    Centroids_array.append(data_points_array[np.random.choice(Index_array)]) # miu 1
    
    for i in range(1, K): # next k-1 centroids
        find_D(D_array, data_points_array, N, Centroids_array) # calculating D_l
        Pr_array = np.array([(D_array[l] / D_array[N]) for l in range(N)]) # Pr_array[i] = D_l / sum(d_l for each 1<=l<=N)
        index = np.random.choice(Index_array,p=Pr_array) # choosing randomly an index of a datapoint to be centroid
        Centroids_array.append(data_points_array[index])
    return Centroids_array


def find_D(D_array, datapoints_array, N, Centroids_array):
    D_array[N] = 0.0
    for l in range(N):
        # D_array[l] = min(x_l - miu_j)^2 for j s.t. 1<=j<=number of chosen centroids
        D_array[l] = np.min([calc(datapoints_array[l][1:], centroid[1:]) for centroid in Centroids_array])
        D_array[N] = D_array[N] + D_array[l] # last cell in D_array is the sum of all D's


def calc(x, y): # used to calculate norm
    z = np.subtract(x, y)
    return np.sum(np.multiply(z, z))


def isfloat(num):# check if a number is a float, used for epsilon
    try:
        float(num)
        return True
    except ValueError:
        return False


def checkInput(given_input, argLen):
    # Validate that the command line arguments are in correct format
    if (argLen < 5) or (argLen > 6):
        print("Invalid Input!")
        sys.exit()

    isValid = True

    # Check each argument
    for i in range(argLen):
        if argLen == 6: # Argument maxIter has been given by user
            # Checking if first 2 arguments (K, maxIter) are integers
            if(i == 1 or i == 2):
                isValid = given_input[i].isnumeric()
                # Checking max_iter>0 and k>0
                if (isValid):
                    if(int(given_input[i]) <= 0):
                        isValid = False
            # checking if epsilon is float and =>0
            if (i == 3):
                isValid = isfloat(given_input[i])
                if(isValid):
                    if(float(given_input[i])<0):
                        isValid = False
        # arglen == 5, argument maxIter has not been given by user 
        else:
            if(i == 1):
                isValid = given_input[i].isnumeric()
                # Checking k>0
                if(isValid):
                    if(int(given_input[i])<=0):
                        isValid=False
            # Checking epsilon is float and =>0
            if(i == 2):
                isValid = isfloat(given_input[i])
                if(isValid):
                    if(float(given_input[i])<0):
                        isValid=False
        # One of the checks above failed
        if(isValid == False):
            print("Invalid Input!")
            sys.exit()


def printOutput(returnFromFit,dimension,K):
    output_res=""
    # Printing first row - indices of K randomly chosen centroids (initial centroids)
    for i in range(K):
        output_res += str(int(returnFromFit[0][i])) # returnFromFit[0] = randomly chosen centroids indices
        if (i != K-1):
            output_res += ","
    output_res += "\n"

    # Printing the K calculated final centroids
    for i in range(K):
        for j in range(dimension):
            output_res+=str('%.4f' %(returnFromFit[1][i][j])) # returnFromFit[1] = centroid got from kmenassp.c
            if (j!=dimension-1):
                output_res += ","
        output_res += "\n"
    print(output_res)

def callFit(N,K, max_iter, epsilon, merge_data, centroids,dimension): #calling fit function from kmeans.c
    # Get arguments for fit function
    datapoints = merge_data.iloc[: , 1:]
    centroids=np.array(centroids)
    centroids_indices=centroids[:,0]
    centroids=centroids[:,1:]
    datapoints_list=datapoints.values.tolist()
    centroids_list=centroids.tolist()

    # N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    try:
        final_centroids = mykmeanssp.fit(N,K,max_iter,datapoints_list,centroids_list,epsilon,dimension)
        return (centroids_indices,final_centroids)
    except:
        print("An Error Has Occurred")
        sys.exit()

def mergeInputs(input_1,input_2):
    try:
        data1= pd.read_csv(input_1, header=None) # Put data from file_1 to a data1
        data2= pd.read_csv(input_2, header=None) # Put data from file_2 to a data2
        mergeDf=pd.merge(data1, data2, on=0)  # Merge (inner join) the two data_frames into one, based on mutual indices
    except:
        print("An Error Has Occurred")
        sys.exit()
    return mergeDf


def main(argv):
    '''=========================Checking input========================='''
    argLen = len(argv)
    checkInput(argv,argLen)

    '''=========================Set arguments========================='''
    K = (int)(argv[1])
    goal = argv[2] # todo enum!
    data_points_array = getNDataPoints(argv[3], goal)
    N = data_points_array.shape[0]
    D=data_points_array.shape[1]

    if K > N:
        print("Invalid Input!")
        sys.exit()

    if(goal != "jacobi" and goal != "spk" and goal != "wam" and goal != "ddg" and goal != "lnorm"):
        print("Invalid Input!")
        sys.exit()

    try:
        goal_matrix = mykmeanssp.fit(N,K,D,data_points_array.tolist(),goal)
        if(goal!="spk"):
            print_output(goal_matrix,N,goal)
        else:
            
    except:
        print("An Error Has Occurred")
        sys.exit()




if __name__ == '__main__':
    main(sys.argv)
