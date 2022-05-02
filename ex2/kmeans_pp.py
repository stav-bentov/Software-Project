
import numpy as np





def kMeans_init(K, maxIter, input_file, epsilon):

    fileOpener = open(input_filename, "r")

    data_points_array = []  # saves the datapoint x1, ... , xN (array DxN)
    Centroids_array = []  # saves the centroids u1, ... , uK
   
    N = 0

    # build data_points_array
    while True:
        line = fileOpener.readline()
        if line == '': #end of file
            break
        data_points_array.append([float(x) for x in line.split(",")]) #each data point is in different line, seperated with commas
        N += 1
    fileOpener.close()
    # end of while loop

    
    D_array = [0 for i in range(N+1)]
    Pr_array = [0 for i in range(N)]
    np.random.seed(0)
    Centroids_array.append(np.random.choice(data_points_array))

    for i in range(K):
        find_D(D_array, data_points_array, N, Centroids_array)
        Pr_array = [(D_array[l] / D_array[N]) for l in range(N)]
        Centroids_array.append(np.random.choice(data_points_array, Pr_array))


    return Centroids_array


def find_D(D_array, datapoints_array, N, Centroids_array):
    for l in range(N):
        D_array[l] = min([np.linalg.norm(datapoints_array[l] - centroid) for centroid in Centroids_array])
        D_array[N] = D_array[N] + D_array[l]
    

    

    

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
                isValid = argv[i].isfloat()
                

        else:# arglen == 5, argument maxIter has not been given by user 
            if(i == 1):
                isValid = argv[i].isnumeric()
            if(i == 2):#checking if epsilon is float
                isValid = argv[i].isfloat()
        if(isValid == False):#one of the checks above failed
            print("Inavlid Input!")
            sys.exit()

    if(max_iter <= 0):
        print("Inavlid Input!")
        sys.exit()

    K = (int)(argv[1])
    if(argLen == 6):
        max_iter = (int)(argv[2])
        return kMeans(K, max_iter, argv[3], argv[4])#todo check!!
    return kMeans(K, max_iter, argv[2], argv[3])#todo check!!


if __name__ == '__main__':
    main(sys.argv)
