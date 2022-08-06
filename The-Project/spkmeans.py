import math
import numpy as np
import pandas as pd
import sys
from enum import Enum
import my_spkmeans


class Goal(Enum):
    wam = 1
    ddg = 2
    lnorm = 3
    jacobi = 4
    spk = 5
    spk_ex2 = 6

# each goal's output (except spk) is a matrix N*N (wam,ddg,lnorm) or (N+1)*N (jacobi)
def print_output(returnFromFit, N, goal):
    output_res = ""
    num_rows = N
    if(goal == Goal.jacobi):
        num_rows += 1

    for i in range(num_rows):
        for j in range(N):
            # todo check with Stav about %
            output_res += str('%.4f' % (returnFromFit[i][j]))
            if (j != N-1):
                output_res += ","
        output_res += "\n"
    print(output_res)

# return data's matrix, N (=number of points/ number of rows), D (=point's dimension/ number of cols)
def get_goal_input(filename):
    try:
        data = pd.read_csv(filename, header=None)
        return (data.to_numpy(), data.shape[0], data.shape[1])
    except:
        print("An Error Has Occurred")
        sys.exit()

# spk's output is like in ex2- K indexs and centroids
def print_output_spk(returnFromFit, K, D, centroids_index):
    output_res = ""
    # Printing first row - indices of K randomly chosen centroids (initial centroids)
    for i in range(K):
        # returnFromFit[0] = randomly chosen centroids indices
        output_res += str(int(centroids_index[i]))
        if (i != K-1):
            output_res += ","
    output_res += "\n"

    # Printing the K calculated final centroids
    for i in range(K):
        for j in range(D):
            # returnFromFit[1] = centroid got from kmenassp.c
            output_res += str('%.4f' % (returnFromFit[i][j]))
            if (j != D-1):
                output_res += ","
        output_res += "\n"
    print(output_res)

# gets the input, exits if: arglen!=4 or k is'nt a number, or k is negative/0/1. return K,goal
def check_input(given_input, argLen):
    if (argLen != 4):
        print("Invalid Input!")
        sys.exit()
    # check K: K isn't an integer or K is a negative integer (K needs to be greater then 1!)
    is_valid = given_input[1].isnumeric()
    if is_valid:
        is_valid = (int(given_input[1]) > 1 or int(given_input[1]) == 0)
    # check enum: one of the 5 given options
    if is_valid:
        is_valid = given_input[2] in [curr_goal.name for curr_goal in Goal]

    if not is_valid:
        print("Invalid Input!")
        sys.exit()

    return ((int)(given_input[1]), Goal[given_input[2]])


''' ========================= spk from ex 2 ========================='''
def kMeans_init(K, data_points_array):

    Centroids_array = []  # Saves the centroids u1, ... , uK
    Centroids__index_array = []

    N = len(data_points_array)
    '''if (K > N):
        print("Invalid Input!")
        sys.exit()'''  # todo- needed?

    D_array = np.array([0.0 for i in range(N+1)])
    # Pr_array[i] = probability of
    Pr_array = np.array([0.0 for i in range(N)])
    Index_array = np.array([i for i in range(N)])
    np.random.seed(0)

    #added diffrently from ex2
    index=np.random.choice(Index_array)
    Centroids__index_array.append(index)

    Centroids_array.append(data_points_array[index])  # miu 1

    for i in range(1, K):  # next k-1 centroids
        find_D(D_array, data_points_array, N,
               Centroids_array)  # calculating D_l
        # Pr_array[i] = D_l / sum(d_l for each 1<=l<=N)
        Pr_array = np.array([(D_array[l] / D_array[N]) for l in range(N)])
        # choosing randomly an index of a datapoint to be centroid
        index = np.random.choice(Index_array, p=Pr_array)
        Centroids__index_array.append(index)
        Centroids_array.append(data_points_array[index])
    return Centroids__index_array,Centroids_array


def find_D(D_array, datapoints_array, N, Centroids_array):
    D_array[N] = 0.0
    for l in range(N):
        # D_array[l] = min(x_l - miu_j)^2 for j s.t. 1<=j<=number of chosen centroids
        D_array[l] = np.min([calc(datapoints_array[l][1:], centroid[1:])
                            for centroid in Centroids_array])
        # last cell in D_array is the sum of all D's
        D_array[N] = D_array[N] + D_array[l]


def calc(x, y):  # used to calculate norm
    z = np.subtract(x, y)
    return np.sum(np.multiply(z, z))


# calling fit function from kmeans.c
def call_fit_ex2(N, K, dimension, data_points, centroids, goal):
    # N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension
    try:
        final_centroids = my_spkmeans.fit(N, K, dimension, data_points, goal.value, centroids)
        return final_centroids
    except:
        print("An Error Has Occurred")
        sys.exit()


''' ========================= DONE spk from ex 2 ========================='''


def main(argv):
    '''=========================Check input and Set arguments========================='''
    argLen = len(argv)
    K, goal = check_input(argv, argLen)
    data_points_array, N, D = get_goal_input(argv[3])

    if K > N:
        print("Invalid Input!")
        sys.exit()

    '''=========================run goal========================='''
    try:
        goal_matrix = my_spkmeans.fit(N, K, D, data_points_array.tolist(), goal.value, [])
        if(goal != Goal.spk):
            print_output(goal_matrix, N, goal)
        else:
            if(K == 0):
                K = len(goal_matrix[0])
            goal=Goal.spk_ex2
            centroids_index, centroids = kMeans_init(K, goal_matrix)
            D=K
            #goal_matrix= N*K, dimension=D=K
            print_output_spk(call_fit_ex2(N, K, D,goal_matrix, centroids, goal), K, D, centroids_index)

    except:
        print("An Error Has Occurred")
        sys.exit()


if __name__ == '__main__':
    main(sys.argv)
