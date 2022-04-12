import math
import sys

from pip import main

EPSILON=0.001

def kMeans(K, max_iter,input_filename, output_filename):

    fileOpener = open(input_filename, "r")
    
    data_points_array = [] # saves the datapoint x1, ... , xN (array DxN)
    kMeans_array = [] # saves the centroids u1, ... , uK
    kMeans_size_array = [0 for i in range(K)] # saves clusters size |S1|, ..., |SK|
    counter = 0

    # build data_points_array
    while True:
        line = fileOpener.readline()
        if line == '':
            break
        data_points_array.append([float(x) for x in line.split(",")])
        if counter < K :
            kMeans_array.append([float(x) for x in line.split(",")])
        counter += 1
    fileOpener.close()
    ### end of while loop

    N=len(data_points_array)

    if(K>N or K <= 0):
        print("Invalid Input!")
        sys.exit()

    dimension=len(data_points_array[0])

    # for each 0<=i<N data_points_cluster[i]=S for s is the cluster of xi
    data_points_cluster = [0 for i in range(N)]
    counter = 0

    # while euclidean norm is larger then x and number of iteration is less then max_iter
    while (check_euclidean_norm(kMeans_array)) and (max_iter > counter) :
        # calculate cluster for each datapoint
        for i in range(N):
            # update each datapoints's cluster
            data_points_cluster[i] = find_cluster(kMeans_array, data_points_array[i]) # update cluster array
            # updtae number of datapoints that belongs to each cluster 
            kMeans_size_array[data_points_cluster[i] - 1] += 1 # update cluster's size

        kMeans_array = [[float(0) for i in range (dimension)] for j in range(K)]
        
        # calculate updated clusters
        for i in range(len(data_points_array)):
            add_vectors(kMeans_array[data_points_cluster[i]-1], data_points_array[i])
        for i in range(len(kMeans_array)):
            div_vectors(kMeans_array[i], kMeans_size_array[i])
        
        kMeans_size_array=[0 for i in range(K)]
        counter+=1

    fileOpener = open(output_filename, "w")
    kMeans_array = [['%.4f' % (kMeans_array[j][i]) for i in range (len(data_points_array[0]))] for j in range(K)]

    for mean in kMeans_array:
        for i in range(len(mean)):
            if i != (len(mean)-1):
                fileOpener.write(str(mean[i]) + ",")
            else:
                fileOpener.write(str(mean[i]) + "\n")

def check_euclidean_norm(kMeans_array):
    for mean in kMeans_array:
        sum = 0
        for i in mean:
            sum += i**2
        # if there exist a centroid that changes more then epsilon- return true 
        # there is no convergence yet
        if math.sqrt(sum) >= EPSILON:
            return True
    return False

# find the corresponding cluster of datapoint
def find_cluster(kMeans_array, dataPoint):
    meanIndex = 0
    minSum = math.inf
    index = 0

    for mean in kMeans_array :
        index += 1
        sum = 0
        
        for i in range(len(dataPoint)):
            sum += (dataPoint[i] - mean[i])**2

        if sum <= minSum:
            minSum = sum
            meanIndex = index
    return meanIndex

def add_vectors(mean, data_point):
    for i in range(len(mean)):
        mean[i] += data_point[i]

def div_vectors(mean, clusterSize):
    if clusterSize != 0:
        for i in range(len(mean)):
            mean[i] = mean[i] / clusterSize


def main(argv):
    max_iter=200
    isValid=True
    argLen=len(argv)
    if(argLen<4 or argLen>5):
        print("Inavlid Input!")
        sys.exit()
    
    for i in range(argLen):
        if(argLen==5):
            if(i==1 or i==2):
                isValid=argv[i].isnumeric()
            if(i==3 or i==4):
                isValid=(argv[i][-4:]==".txt")
        else:
            if(i==1):
                isValid=argv[i].isnumeric()
            if(i==2 or i==3):
                isValid=(argv[i][-4:]==".txt")
        if(isValid==False):
            print("Inavlid Input!")
        sys.exit()
    
    if(max_iter<=0):
        print("Inavlid Input!")
        sys.exit()

    K=(int)(argv[1])
    if(argLen==5):
        max_iter=(int)(argv[2])
        return kMeans(K,max_iter,argv[3],argv[4])
    return kMeans(K,max_iter,argv[2],argv[3])

if __name__ == "__main__":
    main(sys.argv)
