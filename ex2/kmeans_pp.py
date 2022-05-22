
import imp
import numpy as np
import pandas as pd
import sys
#import mykmeanssp

def kMeans_init(K, maxIter, mergeDf, epsilon):

    Centroids_array = []  # saves the centroids u1, ... , uK

    mergeDf=mergeDf.sort_values(by=[0])
    #print(mergeDf)
    data_points_array=mergeDf.to_numpy() # saves the datapoints
    #print(data_points_array)

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
    print("type(Centroids_array):",type(Centroids_array))
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


def main(argv):
    max_iter = 300 #default value 
    isValid = True
    argLen = len(argv)
    if(argLen < 5 or argLen > 6):
        print("Inavlid Input!")
        sys.exit()

    for i in range(argLen):
        if(argLen == 6):#argument maxIter has been given by user 
            if(i == 1 or i == 2):#checking if first 2 arguments (K, maxIter) are integers
                isValid = argv[i].isnumeric()
            if(i == 3):#checking if epsilon is float
                isValid = isfloat(argv[i])
                

        else:# arglen == 5, argument maxIter has not been given by user 
            if(i == 1):
                isValid = argv[i].isnumeric()
            if(i == 2):#checking if epsilon is float
                isValid = isfloat(argv[i])
        if(isValid == False):#one of the checks above failed
            print("Inavlid Input!")
            sys.exit()

    if(max_iter <= 0):
        print("Inavlid Input!")
        sys.exit()

    K = (int)(argv[1])
    if(argLen == 6):
        max_iter = (int)(argv[2])
        merge_data=mergeInputs(argv[4], argv[5])
        epsilon=argv[3]
        centriods= kMeans_init(K, max_iter, merge_data, argv[3])#TODO check!!
    else:
        merge_data=mergeInputs(argv[3], argv[4])
        centriods= kMeans_init(K, max_iter, merge_data, argv[2])#TODO check!!
        epsilon=argv[2]

    #TODO check cases with vectors length 1 
    dimension=merge_data.shape[1]-1
    returnFromFit=callFit(K,max_iter,epsilon,merge_data,centriods,dimension)
    output_res=""

    for i in range(K):
        output_res+=str(int(returnFromFit[0][i]))
        if (i!=K-1):
            output_res+=","
    output_res+="\n"

    for i in range(K):
        for j in range(dimension):
            output_res+=str(int(returnFromFit[1][i][j]))
            if (i!=K-1):
                output_res+=","
        output_res+="\n"

def callFit(K, max_iter, epsilon, merge_data, centroids,dimension):
    
    datapoints = merge_data.iloc[: , 1:]
    centroids=np.array(centroids)
    centroids_indices=centroids[:,0]
    centroids=centroids[:,1:]
    N=len(merge_data)
    datapoints_list=datapoints.values.tolist()
    centroids_list=centroids.tolist()

    #N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    return (centroids_indices,fit(N,K,max_iter,datapoints_list,centroids_list,epsilon,dimension))

def mergeInputs(input_1,input_2):
    data1= pd.read_csv(input_1, header=None)
    data2= pd.read_csv(input_2, header=None)
    mergeDf=pd.merge(data1, data2, on=0)
    mergeDf.columns=[i for i in range(len(mergeDf.columns))]
    return mergeDf

if __name__ == '__main__':
    main(sys.argv)
