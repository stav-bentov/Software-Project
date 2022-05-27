
import numpy as np
import pandas as pd
import sys
import mykmeanssp

# TODO add description to func 
# TODO check errors
# TODO remove all the notes for us


def kMeans_init(K, mergeDf):

    Centroids_array = []  # saves the centroids u1, ... , uK

    mergeDf = mergeDf.sort_values(by=[0])
    data_points_array = mergeDf.to_numpy()  # saves the data points

    N = len(data_points_array)
    if (K > N):
        print("Invalid Input!")
        sys.exit()

    D_array = np.array([0.0 for i in range(N+1)])
    Pr_array = np.array([0.0 for i in range(N)]) # Pr_array[i] = probability of
    Index_array = data_points_array[:, 0].astype(int)
    #Index_array = np.array([i for i in range(N)])# Index_array[i] = ith datapoint in the mergedDf #TODO check if the indices are 1-100 or 3,292,999...
    np.random.seed(0)
    Centroids_array.append(data_points_array[np.random.choice(Index_array)]) # miu 1
    
    for i in range(1, K): # next k-1 centroids
        find_D(D_array, data_points_array, N, Centroids_array) # calculating D_l
        if D_array[N] != 0:  #todo check what to do if it happens
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
    # printing first row - indices of K randomly chosen centroids (initial centroids)
    for i in range(K):
        output_res += str(int(returnFromFit[0][i])) # returnFromFit[0] = randomly chosen centroids indices
        if (i != K-1):
            output_res += ","
    output_res += "\n"

    # printing the K calculated final centroids
    for i in range(K):
        for j in range(dimension):
            output_res+=str('%.4f' %(returnFromFit[1][i][j]))
            if (j!=dimension-1):
                output_res += ","
        output_res += "\n"
    print(output_res)

    fileOpener = open("output_filename.txt", "w")
    fileOpener.write(output_res)
    fileOpener.close()

def callFit(N,K, max_iter, epsilon, merge_data, centroids,dimension): #calling fit function from kmeans.c

    # Get arguments for fit function
    datapoints = merge_data.iloc[: , 1:]
    centroids=np.array(centroids)
    centroids_indices=centroids[:,0]
    centroids=centroids[:,1:]
    datapoints_list=datapoints.values.tolist()
    centroids_list=centroids.tolist()

    # TODO use try and except?
    # N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    try:
        final_centroids = mykmeanssp.fit(N,K,max_iter,datapoints_list,centroids_list,epsilon,dimension)
        return (centroids_indices,final_centroids)
    except:
        print("An Error Has Occurred")
        sys.exit()

def mergeInputs(input_1,input_2):
    data1= pd.read_csv(input_1, header=None) # put data from file_1 to a data_frame 1
    data2= pd.read_csv(input_2, header=None) # put data from file_2 to a data_frame 2
    mergeDf=pd.merge(data1, data2, on=0)  # merge (inner join) the two data_frames into one, based on mutual indices
    #mergeDf.columns=[i for i in range(len(mergeDf.columns))] #Todo check if needed (not only for us)
    return mergeDf

def main(argv): #TODO check if it's okay that K is string (f.e: K = "3")
    '''=========================Checking input========================='''
    argLen = len(argv)
    checkInput(argv,argLen)
    '''=========================Set arguments========================='''
    max_iter = 300 # Default value
    K = (int)(argv[1])
    # Merge data with help of the correct values from argv
    if(argLen == 6):
        max_iter = (int)(argv[2])
        epsilon=float(argv[3])
        merge_data=mergeInputs(argv[4], argv[5])
        N = len(merge_data)
        centroids= kMeans_init(K, merge_data)
    else:
        epsilon=float(argv[2])
        merge_data=mergeInputs(argv[3], argv[4])
        N = len(merge_data)
        centroids= kMeans_init(K, merge_data)

    #TODO check cases with vectors length 1 
    dimension=merge_data.shape[1]-1
    '''===================Done checking and getting input=================='''

    '''====================Get final centroids and print==================='''
    # TODO check if good incase of an error- maybe separate to 2 parts
    printOutput(callFit(N,K,max_iter,epsilon,merge_data,centroids,dimension),dimension,K)

if __name__ == '__main__':
    main(sys.argv)


