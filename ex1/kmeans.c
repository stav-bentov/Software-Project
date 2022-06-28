#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#define INVALID "Invalid Input!"
#define ERROR "An Error Has Occurred"

int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K);
int find_cluster(double **Centroids, double *Datapoint, int dimension, int K);
void free_memory(double **Centroids, double **Datapoints, double **oldCentroids, int K, int N);
void kMeans(int K, int max_iter,  char *input_filename,  char *output_filename);
void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K);
int is_number(const char *argument);

static const double EPSILON = 0.001;

int main(int argc, char *argv[])
{

    int i, isValid, K, max_iter;

    /* default of maximum iterations*/
    max_iter = 200;

    /* if the input is valid-1, else-0*/
    isValid = 1;

    /* invalid number of arguments*/
    if (argc < 4 || argc > 5)
    {
        printf(INVALID);
        exit(1);
    }

    /* in case of 5 arguments- there is max_iter,
        in case of 4- no max_iter, so check properly*/
    for (i = 1; i < argc; i++)
    {
        if (argc == 5)
        {
            if (i == 1 || i == 2)
            {
                isValid = is_number(argv[i]);
            }
        }
        else
        {
            if (i == 1)
            {
                isValid = is_number(argv[i]);
            }
        }

        /* if one of the arguments isn't a number*/
        if (!isValid)
        {
            printf(INVALID);
            exit(1);
        }
    }

    K = atoi(argv[1]);
    if (argc == 5)
    {
        max_iter = atoi(argv[2]);
        if (max_iter <= 0)
        {
            printf(INVALID);
            exit(1);
        }
        kMeans(K, max_iter, argv[3], argv[4]);
        exit(0);
    }

    kMeans(K, max_iter, argv[2], argv[3]);
    exit(0);
}

void kMeans(int K, int max_iter, char *input_filename, char *output_filename)
{
    /*
    ifp= file's pointer.
    c= pointer to char while scanning the input\output size.
    N= number of vectors.
    dimension= vector's dimension.
    firstLine= helps knowing when we are scanning the first line and only then calculate the dimension.
    i, j, counter= counters for loop iterations.
    corr= saves each string we scan and helps converting to double.
    Datapoints= saves all datapoint's vectors.
    Centroids= saves all centroid's vectors (updated).
    oldCentroids= saves all centroid's vectors (before change).
    */
    FILE *ifp = NULL;
    char c;
    int N, dimension;
    int firstLine;
    int i, j, counter;
    double corr;
    double **Datapoints;
    double **Centroids;
    double **oldCentroids;
    int cluster;

    i = 0;
    j = 0;
    dimension = 1;
    N = 0;
    firstLine = 1;
    counter = 0;

    ifp = fopen(input_filename, "r");
    if (ifp == NULL)
    {
        printf(ERROR);
        exit(1);
    }

    /* count number of rows= number of vectors, number of ',' in first lise= dimension*/
    while ((c = fgetc(ifp)) != EOF)
    {
        if (firstLine)
            if (c == ',')
                dimension++;
        if (c == '\n')
        {
            N++;
            if (firstLine)
            {
                firstLine = 0;
            }
        }
    }

    if (K > N || K <= 0)
    {
        printf(INVALID);
        exit(1);
    }

    rewind(ifp);
    Datapoints = malloc((sizeof(double *)) * N);
    Centroids = malloc((sizeof(double *)) * K);
    oldCentroids = malloc((sizeof(double *)) * K);

    if (Datapoints == NULL || Centroids == NULL || oldCentroids == NULL)
    {
        printf(ERROR);
        exit(1);
    }

    /* set arrays and then insert datapoints and centroids*/
    for (i = 0; i < N; i++)
    {
        /* last cell contains the datapoint's cluster*/
        Datapoints[i] = malloc((sizeof(double)) * (dimension + 1));
        if (Datapoints[i] == NULL)
        {
            printf(ERROR);
            exit(1);
        }
        if (i < K)
        {
            /* last call contains number of datapoints assigned to centroid's cluster*/
            Centroids[i] = malloc((sizeof(double)) * (dimension + 1));
            oldCentroids[i] = malloc((sizeof(double)) * (dimension));
            if (Centroids[i] == NULL || oldCentroids[i] == NULL)
            {
                printf(ERROR);
                exit(1);
            }
        }

        /* insert vectors coordinate*/
        for (j = 0; j < dimension; j++)
        {
            if (fscanf(ifp, "%lf", &corr) == 1)
            {
                Datapoints[i][j] = (double)corr;
                if (i < K)
                {
                    Centroids[i][j] = (double)corr;
                    oldCentroids[i][j] = (double)corr;
                }
                fgetc(ifp);
            }
        }

        /* zero in last cell [dimension+1]*/
        Datapoints[i][j] = 0;
        if (i < K)
        {
            Centroids[i][j] = 0;
        }
    }

    /* done reading*/
    fclose(ifp);

    /* calculate k-means:
    stop iteration when number of iteration is more then man_iter
    or when all of the centroids have changed less then epsilon*/
    while (max_iter > counter)
    {
        for (i = 0; i < N; i++)
        {
            cluster = find_cluster(Centroids, Datapoints[i], dimension, K);

            /* update for each datapoint- what number of cluster it belongs
            and for each centroids update the number of datapoints that belong to it*/
            Datapoints[i][dimension] = cluster;
            Centroids[cluster - 1][dimension] += 1;
        }

        /* set all centroids to zero so that updated centroid can be calculated next*/
        for (i = 0; i < K; i++)
        {
            for (j = 0; j < dimension; j++)
            {
                Centroids[i][j] = 0;
            }
        }

        /* update centroids according to the calculations*/
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

    ifp = fopen(output_filename, "w");
    if (ifp == NULL)
    {
        printf(ERROR);
        exit(1);
    }

    /* writes k-means output file*/
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            if (j == dimension - 1)
            {
                fprintf(ifp, "%.4f", Centroids[i][j]);
            }
            else
                fprintf(ifp, "%.4f,", Centroids[i][j]);
        }
        fprintf(ifp, "\n");
    }

    fclose(ifp);

    free_memory(Centroids, Datapoints, oldCentroids, K, N);
}

/* gets the new and old centroids, return 1 if all of the centroids didn't change more then epsilon,else-0*/
int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K)
{
    int i, j;
    double sum;

    /*calculate euclidean norm for each centroid*/
    for (i = 0; i < K; i++)
    {
        sum = 0;
        for (j = 0; j < dimension; j++)
        {
            sum += pow(newCentroids[i][j] - oldCentroids[i][j], 2);
        }

        /* one centroid changed more then epsilon*/
        if (sqrt(sum) >= EPSILON)
            return 0;
    }
    /* every centroids changed less then epsilon */
    return 1;
}

/* gets the centroids and one datapoint and return datapoint's cluster*/
int find_cluster(double **Centroids, double *Datapoint, int dimension, int K)
{
    int i, j, cluster;
    double sum, minSum;

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

            /* cluster number i+1 because it represented by index cell i*/
            cluster = i + 1;
        }
    }

    return cluster;
}

/* get all 3 arrays and free all memory that was in use*/
void free_memory(double **Centroids, double **Datapoints, double **oldCentroids, int K, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        if (i < K)
        {
            free(Centroids[i]);
            free(oldCentroids[i]);
        }
        free(Datapoints[i]);
    }
    free(Datapoints);
    free(Centroids);
    free(oldCentroids);
}

/* gets the updated centroids and old ones- update the old centroids*/
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

/* checks if the string that has been given in input is a valid number*/
int is_number(const char *argument)
{
    int i;
    char c;
    i = 0;
    c = argument[0];
    while (c != '\0')
    {
        i++;
        if (!isdigit(c))
            return 0;
        c = argument[i];
    }
    return 1;
}