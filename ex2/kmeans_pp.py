
import numpy as np
import pandas as pd
import sys
import mykmeanssp

# TODO add description to func 
# TODO check errors
# TODO remove all the notes for us
# TODO check every function uses all inputs

def kMeans_init(K, maxIter, mergeDf, epsilon):

    Centroids_array = []  # saves the centroids u1, ... , uK

    mergeDf=mergeDf.sort_values(by=[0])
    data_points_array=mergeDf.to_numpy() # saves the datapoints

    N = len(data_points_array)

    D_array = np.array([0 for i in range(N+1)])
    Pr_array = np.array([0 for i in range(N)])
    Index_array=np.array([i for i in range(N)])
    np.random.seed(0)
    Centroids_array.append(data_points_array[np.random.choice(Index_array)])
    
    for i in range(1,K):
        find_D(D_array, data_points_array, N, Centroids_array)
        Pr_array = np.array([(D_array[l] / D_array[N]) for l in range(N)])
        index=np.random.choice(Index_array,p=Pr_array)
        Centroids_array.append(data_points_array[index])
    return Centroids_array

def find_D(D_array, datapoints_array, N, Centroids_array):
    D_array[N]=0
    for l in range(N):
        D_array[l] = min([calc(datapoints_array[l][1:],centroid[1:]) for centroid in Centroids_array])
        D_array[N] = D_array[N] + D_array[l]

def calc(x,y):
    z=np.subtract(x, y)
    return np.sum(np.multiply(z,z))

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def checkInput(input,argLen):
    # Validate that the command line arguments are in correct format
    if(argLen < 5 or argLen > 6):
        print("Inavlid Input!")
        sys.exit()

    isValid = True

    # Check each argument
    for i in range(argLen):
        #Argument maxIter has been given by user 
        if(argLen == 6):
            #Checking if first 2 arguments (K, maxIter) are integers
            if(i == 1 or i == 2):
                isValid = input[i].isnumeric()
                #Checking max_iter>0 and k>0
                if(isValid):
                    if(int(input[i])<=0):
                        isValid=False
            #checking if epsilon is float and =>0
            if(i == 3):
                isValid = isfloat(input[i])
                if(isValid):
                    if(float(input[i])<0):
                        isValid=False
        # arglen == 5, argument maxIter has not been given by user 
        else:
            if(i == 1):
                isValid = input[i].isnumeric()
                #Checking k>0
                if(isValid):
                    if(int(input[i])<=0):
                        isValid=False
            #Checking epsilon is float and =>0
            if(i == 2):
                isValid = isfloat(input[i])
                if(isValid):
                    if(float(input[i])<0):
                        isValid=False
        #One of the checks above failed
        if(isValid == False):
            print("Inavlid Input!")
            sys.exit()

def printOutput(returnFromFit,dimension,K):
    output_res=""
    for i in range(K):
        output_res+=str(int(returnFromFit[0][i]))
        if (i!=K-1):
            output_res+=","
    output_res+="\n"

    for i in range(K):
        for j in range(dimension):
            # TODO make sure that if there are zeros- ommit them
            output_res+=str('%.4f' %(returnFromFit[1][i][j]))
            if (j!=dimension-1):
                output_res+=","
        output_res+="\n"
    print(output_res)

def callFit(N,K, max_iter, epsilon, merge_data, centroids,dimension):

    # Get arguments for fit function
    datapoints = merge_data.iloc[: , 1:]
    centroids=np.array(centroids)
    centroids_indices=centroids[:,0]
    centroids=centroids[:,1:]
    datapoints_list=datapoints.values.tolist()
    centroids_list=centroids.tolist()

    # TODO use try and except?
    #N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    return (centroids_indices,mykmeanssp.fit(N,K,max_iter,datapoints_list,centroids_list,epsilon,dimension))

def mergeInputs(input_1,input_2):
    data1= pd.read_csv(input_1, header=None)
    data2= pd.read_csv(input_2, header=None)
    mergeDf=pd.merge(data1, data2, on=0)
    mergeDf.columns=[i for i in range(len(mergeDf.columns))]
    return mergeDf

def main(argv):

    '''=========================Checking input========================='''
    argLen = len(argv)
    checkInput(argv,argLen)
    '''=========================Set arguments========================='''
    max_iter = 300 #Default value
    K = (int)(argv[1])
    # Merge data with help of the correct values from argv
    if(argLen == 6):
        max_iter = (int)(argv[2])
        epsilon=float(argv[3])
        merge_data=mergeInputs(argv[4], argv[5])
        centriods= kMeans_init(K, max_iter, merge_data, argv[3])#TODO check!!
    else:
        epsilon=float(argv[2])
        merge_data=mergeInputs(argv[3], argv[4])
        centriods= kMeans_init(K, max_iter, merge_data, argv[2])#TODO check!!
    N=len(merge_data)
    if (K > N):
        print("Inavlid Input!")
        sys.exit()
    #TODO check cases with vectors length 1 
    dimension=merge_data.shape[1]-1
    '''===================Done checking and getting input=================='''

    '''====================Get final centroids and print==================='''
    # TODO check if good incase of an error- maybe seperete to 2 parts
    printOutput(callFit(N,K,max_iter,epsilon,merge_data,centriods,dimension),dimension,K)

if __name__ == '__main__':
    main(sys.argv)
