#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K);
int find_cluster(double **Centroids, double *Datapoint, int dimension, int K);
void free_memory(double **Centroids, double **Datapoints, int K, int N);
int kMeans(int K, int max_iter, const char *input_filename, const char *output_filename);
void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K);
int is_number(const char *argument);
int is_txt(const char *argument);

static const double EPSILON = 0.001;

int main(int argc, char const *argv[])
{

    int i, isValid, K, max_iter;
    max_iter = 200;
    isValid = 1;

    /* invalid number of arguments*/
    if (argc < 4 || argc > 5)
    {
        printf("Invalid Input!");
        return 1;
    }

    for (i = 1; i < argc; i++)
    {
        if (argc == 5)
        {
            if (i == 1 || i == 2)
            {
                isValid = is_number(argv[i]);
            }
            if (i == 3 || i == 4)
            {
                isValid = is_txt(argv[i]);
            }
        }
        else
        {
            if (i == 1)
            {
                isValid = is_number(argv[i]);
            }
            if (i == 2 || i == 3)
            {
                isValid = is_txt(argv[i]);
            }
        }

        if (!isValid)
        {
            printf("Invalid Input!");
            return 1;
        }
    }

    K = atoi(argv[1]);
    if (argc == 5)
    {
        max_iter = atoi(argv[2]);
        return kMeans(K, max_iter, argv[3], argv[4]);
    }

    return kMeans(K, max_iter, argv[2], argv[3]);
}

int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K)
{
    int i, j;
    double sum;
    for (i = 0; i < K; i++)
    {
        sum = 0;
        for (j = 0; j < dimension; j++)
        {
            sum += pow(newCentroids[i][j] - oldCentroids[i][j], 2);
        }
        if (sqrt(sum) >= EPSILON)
            return 0;
    }
    /* every centroids changed doesn't changed more then epsilon */
    return 1;
}

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

            /* cluster number i+1 because in cell number i*/
            cluster = i + 1;
        }
    }

    return cluster;
}

void free_memory(double **Centroids, double **Datapoints, int K, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        if (i < K)
        {
            free(Centroids[i]);
        }
        free(Datapoints[i]);
    }
    free(Datapoints);
    free(Centroids);
}

int kMeans(int K, int max_iter, const char *input_filename, const char *output_filename)
{
    /*
    c= pointer to char while scanning the input\ output size
    N= number of vectors
    dimension= vector's dimension
    firstLine= helps know when we are scanning the first line and only then calculate the dimension
    i, j, counter= counters for loop iteretions
    corr= saves each string we scan and helps converting to double
    Datapoints= saves all datapoint's vectors
    Centroids= saves all centroid's vectors (updated)
    oldCentroids= saves all centroid's vectors (before change)
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

    if (max_iter <= 0)
    {
        printf("Invalid Input!");
        return 1;
    }

    i = 0;
    j = 0;
    dimension = 1;
    N = 0;
    firstLine = 1;
    counter = 0;
    ifp = fopen(input_filename, "r");

    if (ifp == NULL)
    {
        printf("An Error Has Occurred");
        return 1;
    }

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
        printf("Invalid Input!");
        return 1;
    }

    rewind(ifp);
    Datapoints = malloc((sizeof(double *)) * N);
    Centroids = malloc((sizeof(double *)) * K);
    oldCentroids = malloc((sizeof(double *)) * K);

    if (Datapoints == NULL || Centroids == NULL)
    {
        printf("An Error Has Occurred");
        return 1;
    }

    for (i = 0; i < N; i++)
    {
        /* last cell contains the datapoint's cluster*/
        Datapoints[i] = malloc((sizeof(double)) * (dimension + 1));
        if (Datapoints[i] == NULL)
        {
            printf("An Error Has Occurred");
            return 1;
        }
        if (i < K)
        {
            /* last call contains number of datapoints assigned to centroid's cluster*/
            Centroids[i] = malloc((sizeof(double)) * (dimension + 1));
            oldCentroids[i] = malloc((sizeof(double)) * (dimension));
            if (Centroids[i] == NULL)
            {
                printf("An Error Has Occurred");
                return 1;
            }
        }

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

        Datapoints[i][j] = 0;
        if (i < K)
        {
            Centroids[i][j] = 0;
            /*Centroids[i][j+1]=0;*/
        }
    }

    fclose(ifp);

    while (max_iter > counter)
    {
        for (i = 0; i < N; i++)
        {
            cluster = find_cluster(Centroids, Datapoints[i], dimension, K);
            Datapoints[i][dimension] = cluster;
            Centroids[cluster - 1][dimension] += 1;
        }

        for (i = 0; i < K; i++)
        {
            for (j = 0; j < dimension; j++)
            {
                Centroids[i][j] = 0;
            }
        }

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
                Centroids[i][j] = Centroids[i][j] / Centroids[i][dimension];
            }
            Centroids[i][dimension] = 0;
        }

        if (check_euclidean_norm(Centroids, oldCentroids, dimension, K))
        {
            break;
        }

        updateOldCentroid(Centroids, oldCentroids, dimension, K);

        counter++;
    }

    ifp = fopen(output_filename, "w");
    if (ifp == NULL)
    {
        printf("An Error Has Occurred");
        return 1;
    }

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

    free_memory(Centroids, Datapoints, K, N);

    return 0;
}

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

int is_txt(const char *argument)
{
    int i, j;
    char c;
    char strArg[4];

    char str[] = ".txt";
    i = 0;
    j = 0;
    c = argument[0];

    while (c != '\0')
    {
        i++;
        c = argument[i];
    }
    if (i < 5)
    {
        return 0;
    }

    i = i - 4;

    while (j <= 4)
    {
        c = argument[i + j];
        strArg[j] = c;
        j++;
    }

    if (strcmp(strArg, str))
    {
        return 0;
    }

    return 1;
}