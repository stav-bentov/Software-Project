#include "spkmeans.h"
/*#include <Python.h>*/

int kMeans(int N, int K, double **Datapoints, double **Centroids, int dimension)
{
    /*
    i, j, counter= counters for loop iterations.
    oldCentroids= saves all centroid's vectors (before change).
    */
    int i, j, counter;
    double **oldCentroids;
    int cluster;

    counter = 0;

    oldCentroids = malloc((sizeof(double *)) * K);
    if (oldCentroids == NULL)
    {
        return FAIL;
    }
    for(i=0;i<K;i++)
    {
        oldCentroids[i] = malloc((sizeof(double)) * (dimension));
        if (oldCentroids[i] == NULL)
        {
            return FAIL;
        }
        for(j=0; j<dimension; j++)
        {
            oldCentroids[i][j]=Centroids[i][j];
        }
    }

    /* Calculate k-means:
    stop iteration when number of iteration is more then man_iter
    or when all of the centroids have changed less then epsilon*/
    while (MAX_ITER_KMEANS > counter)
    {
        for (i = 0; i < N; i++)
        {
            cluster = find_cluster(Centroids, Datapoints[i], dimension, K);

            /* Update for each datapoint- what number of cluster it belongs
            and for each centroids update the number of datapoints that belong to it*/
            Datapoints[i][dimension] = cluster;
            Centroids[cluster - 1][dimension] += 1;
        }

        /* Set all centroids to zero so that updated centroids can be calculated next*/
        for (i = 0; i < K; i++)
        {
            for (j = 0; j < dimension; j++)
            {
                Centroids[i][j] = 0;
            }
        }

        /* Update centroids according to the calculations*/
        for (i = 0; i < N; i++)
        {
            cluster = Datapoints[i][dimension];
            for (j = 0; j < dimension; j++)
            {
                Centroids[cluster - 1][j] += Datapoints[i][j];
            }
        }
        for (i = 0; i < K; i++)
        {
            for (j = 0; j < dimension; j++)
            {
                /*Centroids[i][j]= sum of datapoints that belong to cluster i+1
                Centroids[i][dimension]= number of datapoints that belong to cluster i*/
                Centroids[i][j] = Centroids[i][j] / Centroids[i][dimension];
            }
            Centroids[i][dimension] = 0;
        }

        /* if all centroids changed less then epsilon -done, else- countinue*/
        if (check_euclidean_norm(Centroids, oldCentroids, dimension, K))
        {
            break;
        }

        /* make oldcentroids be the new ones for next iteration*/
        updateOldCentroid(Centroids, oldCentroids, dimension, K);

        counter++;
    }

    free_memory(oldCentroids, K);
    return SUCCESS;
}

/* Gets the new and old centroids, return 1 if all of the centroids didn't change more then epsilon,else-0*/
int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K)
{
    int i, j;
    double e_norm;

    /*Calculate euclidean norm for each centroid*/
    for (i = 0; i < K; i++)
    {
        e_norm = calc_euclidean_norm(newCentroids[i],oldCentroids[i],dimension);
        /* One centroid changed more then epsilon*/
        if (e_norm >= EPSILON_KMEANS)
            return 0;
    }
    /* Every centroids changed less then epsilon */
    return 1;
}

/* Gets the centroids and one datapoint and return datapoint's cluster*/
int find_cluster(double **Centroids, double *Datapoint, int dimension, int K)
{
    int i, j, cluster;
    double sum, minSum;

    cluster=0; /*Default*/

    minSum = DBL_MAX;
    for (i = 0; i < K; i++)
    {
        sum = 0;
        for (j = 0; j < dimension; j++)
        {
            sum += pow((Datapoint[j] - Centroids[i][j]), 2);
        }
        if (minSum >= sum)
        {
            minSum = sum;

            /* Cluster number i+1 because it represented by index cell i*/
            cluster = i + 1;
        }
    }

    return cluster;
}

/* Gets the updated centroids and old ones- update the old centroids*/
void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K)
{
    int i, j;
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            oldCentroids[i][j] = newCentroids[i][j];
        }
    }
}