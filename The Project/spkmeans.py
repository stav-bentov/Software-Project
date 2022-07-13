import math

import numpy as np
import pandas as pd
import sys
import mykmeanssp

# exclusively parameter is used to determine if we want to perform just this goal or maybe to use the function's result


# '''=========================JACOBI========================='''
def call_jacobi(N, K, A, max_iter=100, exclusively=False):
    epsilon = 0.00001  # epsilon = 1.0 x 10^-5
    try:
        res = mykmeanssp.fit_jacobi(N, K, max_iter, A, epsilon)
    except:
        print("An Error Has Occurred")
        sys.exit()
    # todo ooooooooooo increasing sorting of eigenvalues , vectors
    if exclusively:  # otherwise, we only want the matrix it performs without printing
        printOutput_jacobi(res)
    return res  # returns eigenValues (1st row) and eigenVectors (2nd row onwards)


def printOutput_jacobi(returnFromFit):  # todo check if jacobi also requires commas and also print as columns or rows!
    output_res = ""
    no_of_rows = len(returnFromFit)

    for i in range(no_of_rows):
        len_of_row_i = len(returnFromFit[i])
        for j in range(len_of_row_i):
            output_res += str('%.4f' % (returnFromFit[i][j]))  # todo check with Stav about %
            if (j != len_of_row_i-1):
                output_res += ","
        output_res += "\n"
    print(output_res)


    # '''=========================SPK========================='''

def call_spk(N, K, data_points_array):
    # steps 1+2
    Lnorm = call_lnorm(N, K, data_points_array)  # todo - call_wam (step 1 in algorithm) will be called by call_lnorm which will use it as W

    # steps 3+4
    U = call_jacobi(N, K, Lnorm)[1:, :K] # K eigenvectors only (first row was eigenvalues and 0 - k-1 columns), todo check that each eigenvector is in his own column |||


    # todo in numpy Tij = Uij / (sum(Uij^2))^1/2
    # step 5
    T = U.apply(lambda Uij: create_T(U, Uij))  # todo oooo create T out of U
    # adding observation index in the left column of T (0 - N-1), like in EX2
    new_col = np.arange(N)
    T.insert(loc=0, column='observation index', value=new_col)

    # step 6
    kMeans(N, K, T, 300) # call kMeans algorithm from EX2 but with the matrix T instead of filesname



def create_T(U, Uij):  # used to create T out of U (step 5)
    pass

    # '''=========================kMEANS from EX2 (part of SPK)========================='''


def kMeans(N, K, T_matrix, max_iter):  # todo check about epsilon : epsilon = 0.00001 or something else?
    np.random.seed(0)
    random_centroids = kMeans_init(N, K, T_matrix) # randomly picks centroids like kmeans_pp in EX2
    printOutput_spk(callFit(N, K, max_iter, epsilon, T_matrix, random_centroids, K), K, K)


def kMeans_init(N, K, mergeDf):
    np.random.seed(0)
    Centroids_array = []  # Saves the centroids u1, ... , uK
    data_points_array = mergeDf.to_numpy()  # Saves the data points

    D_array = np.array([0.0 for i in range(N+1)])
    Pr_array = np.array([0.0 for i in range(N)]) # Pr_array[i] = probability of
    Index_array = np.array([i for i in range(N)])

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


def calc(x, y):  # used to calculate norm
    z = np.subtract(x, y)
    return np.sum(np.multiply(z, z))

def printOutput_spk(returnFromFit,dimension,K):
    output_res=""
    # Printing first row - indices of K randomly chosen centroids (initial centroids)
    for i in range(K):
        output_res += str(int(returnFromFit[0][i]))  # returnFromFit[0] = randomly chosen centroids indices
        if (i != K-1):
            output_res += ","
    output_res += "\n"

    # Printing the K calculated final centroids
    for i in range(K):
        for j in range(dimension):
            output_res+=str('%.4f' %(returnFromFit[1][i][j]))  # returnFromFit[1] = centroid got from kmenassp.c
            if (j!=dimension-1):
                output_res += ","
        output_res += "\n"
    print(output_res)


def callFit(N,K, max_iter, epsilon, merge_data, centroids,dimension):  # calling fit function from kmeans.c
    # Get arguments for fit function
    datapoints = merge_data.iloc[:, 1:]
    centroids=np.array(centroids)
    centroids_indices=centroids[:, 0]
    centroids=centroids[:, 1:]
    datapoints_list=datapoints.values.tolist()
    centroids_list=centroids.tolist()

    # N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    try:
        #important todo: should we just copy kmeans.c (from EX2) to our new C file? or add it as a different file?
        final_centroids = mykmeanssp.fit_spk(N,K,max_iter,datapoints_list,centroids_list,epsilon,dimension)
        return (centroids_indices, final_centroids)
    except:
        print("An Error Has Occurred")
        sys.exit()

# '''=========================END OF SPK========================='''


    # '''=========================WAM========================='''
def call_wam(N, K, data_points_array, exclusively = False):
    return [[5, 2], [24, 5]]


    # '''=========================DDG========================='''
def call_ddg(N, K, data_points_array, exclusively = False):
    return [[5, 2], [24, 5]]


    # '''=========================LNORM========================='''
def call_lnorm(N, K, data_points_array, exclusively = False):  # needs W_MATRIX, returns Lnorm matrix and prints it if exclusively == True
    W = call_wam(N, K, data_points_array)
    # D^-1/2 = (Diagonal Degree Matrix)^-1/2
    return [[5, 2], [24, 5]]


    # '''=========================heuristic method========================='''
def findK(N, data_points_array):
    res = call_jacobi(N, 5,call_lnorm(N, 5, data_points_array))  # K doesn't matter for call_jacobi, input matrix is Lnorm
    eigenValues = res[0] #first row of call_jacobi's returned value represents the eigenValues of the input matrix
    eigenValues = eigenValues[::-1] # change from an increasing order to a decreasing order
    K = -math.inf # todo change - be careful
    for i in range(math.floor(eigenValues/2)):
        delta_i = math.fabs(eigenValues[i]-eigenValues[i+1])
        if (delta_i > K): # if equal (==) we use lowest index
            K = delta_i  # K = argmax(delta_i)
    return K


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


def main(argv):
    '''=========================Checking input========================='''
    argLen = len(argv)
    checkInput(argv,argLen)

    '''=========================Set arguments========================='''
    K = (int)(argv[1])
    goal = argv[2] # todo enum!
    data_points_array = getNDataPoints(argv[3], goal)
    N = len(data_points_array)

    if K == 0:
        K = findK(N, data_points_array) # heuristic method in case that K is 0

    if K > N:
        print("Invalid Input!")
        sys.exit()

    if goal == "jacobi":
        call_jacobi(N, K, data_points_array, 100, True) #todo check whats the matrix to pass : Lnorm or just the matrix we get from user?
    elif goal == "spk":
        call_spk(N, K, data_points_array)
    elif goal == "wam":
        call_wam(N, K, data_points_array, True)
    elif goal == "ddg":
        call_ddg(N, K, data_points_array, True)
    elif goal == "lnorm":
        call_lnorm(N, K, data_points_array, True)
    else: # goal is invalid
        print("Invalid Input!")
        sys.exit()

        # maybe todo dimension=data_points_array.shape[1]-1



if __name__ == '__main__':
    main(sys.argv)
