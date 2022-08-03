import math

import numpy as np
import pandas as pd
import sys
from enum import Enum
import mykmeanssp

class goal(Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5

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
        data = pd.read_csv(filename, header=None)  # Put data from file to a data array
        return data.to_numpy()
    except:
        print("An Error Has Occurred")
        sys.exit()

def print_output_spk(returnFromFit,dimension,K):
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
    # check K: K isn't an integer or K is a negative integer
    is_valid = given_input[1].isnumeric()
    if is_valid:
        is_valid = int(given_input[1]) < 0
    # check enum: one of the 5 given options
    if is_valid:
        is_valid = given_input[2] in [curr_goal.name for curr_goal in goal]

    if not is_valid:
        print("Invalid Input!")
        sys.exit()

''' ========================= spk from ex 2 ========================='''
def kMeans_init(K, data_points_array):

    Centroids_array = []  # Saves the centroids u1, ... , uK

    N = len(data_points_array)
    '''if (K > N):
        print("Invalid Input!")
        sys.exit()''' #todo- needed?

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

def callFit(N, K, dimension, data_points, centroids, goal): #calling fit function from kmeans.c
    # Get arguments for fit function
    centroids=np.array(centroids)
    centroids_indices=centroids[:,0]
    centroids=centroids[:,1:]
    centroids_list=centroids.tolist()

    # N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    try:
        final_centroids = mykmeanssp.fit(N,K,dimension,data_points,goal,centroids_list)
        return (centroids_indices,final_centroids)
    except:
        print("An Error Has Occurred")
        sys.exit()

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

    try:
        goal_matrix = mykmeanssp.fit(N,K,D,data_points_array.tolist(),goal)
        if(goal!="spk"):
            print_output(goal_matrix,N,goal)
        else:
            if(K==0):
                K=len(goal_matrix[0])
            centroids= kMeans_init(K, goal_matrix)
            print_output_spk(callFit(N,K,data_points_array,centroids,K),K,K)
            
    except:
        print("An Error Has Occurred")
        sys.exit()




if __name__ == '__main__':
    main(sys.argv)
