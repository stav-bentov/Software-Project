#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/*#define EPSILON = 0.001*/
static const double EPSILON=0.001;

int check_euclidean_norm(double **Centroids ,int dimension, int K)
{
    int i,j;
    double sum;
    for(i=0;i<K;i++)
    {
        sum=0;
        for(j=0;j<dimension;j++)
        {
            sum+=pow(Centroids[i][j],2);
        }
        if(sqrt(sum)>=EPSILON)
            return 1;
    }
    return 0;
}

int find_cluster(double **Centroids,double *Datapoint,int dimension,int K)
{
    int i,j,cluster;
    double sum,minSum;

    minSum=DBL_MAX;
    for(i=0;i<K;i++)
    {
        sum=0;
        for(j=0;j<dimension;j++)
        {
            sum+=pow((Datapoint[j]-Centroids[i][j]),2);
        }
        if(minSum>=sum)
        {
            minSum=sum;
            
            /* cluster number i+1 because in cell number i*/
            cluster=i+1;
        }
    }
    
    return cluster;
}

/* make default max_iter=200*/
int kMeans(int K, char *input_filename, char *output_filename, int max_iter)
{
    FILE *ifp=NULL;
    char c;
    /*char *line;*/
    int N,dimension;
    int firstLine;
    int i,j,counter;
    double corr;
    /*double *ptr;*/
    double **Datapoints,**Centroids;
    int cluster;

    i=0;
    j=0;
    dimension=1;
    N=0;
    firstLine=1;
    ifp= fopen(input_filename,"r");

    if (ifp == NULL)
    {
        printf("Error opening the file %s", input_filename);
        return -1;
    }

    while((c=fgetc(ifp))!=EOF)
    {
        if(firstLine)
            if(c==',')
                dimension++;
        if(c=='\n')
        {
            N++;
            if(firstLine)
            {
                printf("N= %d",N);
                firstLine=0;
            }
        }
    }
    printf("N= %d",N);
    printf("d= %d",dimension);

    rewind(ifp);
    Datapoints=malloc((sizeof(double*))*N);
    Centroids=malloc((sizeof(double*))*K);
    for(i=0;i<N;i++)
    {
        /* last cell contains the datapoint's cluster*/
        Datapoints[i]=malloc((sizeof(double))*(dimension+1));
        if(i<K)
        {
            /* last call contains number of datapoints assigned to centroid's cluster*/
            Centroids[i]=malloc((sizeof(double))*(dimension+1));
        }
        
        for(j=0;j<dimension;j++)
        {
            if(fscanf(ifp,"%lf",&corr)==1)
            {
                Datapoints[i][j]=(double)corr;
                if(i<K)
                {
                    Centroids[i][j]=(double)corr;
                }
                fgetc(ifp);
            }
        }

        Datapoints[i][j]=0;
        if(i<K)
        {
            Centroids[i][j]=0;
            Centroids[i][j+1]=0;
        }
    }

    fclose(ifp);

    while((check_euclidean_norm(Centroids,dimension,K)) && (max_iter>counter))
    {
        for(i=0;i<N;i++)
        {
            cluster=find_cluster(Centroids,Datapoints[i],dimension,K);
            Datapoints[i][dimension]=cluster;
            Centroids[cluster-1][dimension]+=1;
        }

        for(i=0;i<K;i++)
        {
            for(j=0;j<dimension;j++)
            {
                Centroids[i][j]=0;
            }
        }

        for(i=0;i<N;i++)
        {
            cluster=Datapoints[i][dimension];
            for(j=0;j<dimension;j++)
            {
                Centroids[cluster-1][j]+=Datapoints[i][j];
            }
        }
    
        for(i=0;i<K;i++)
        {
            for(j=0;j<dimension;j++)
            {
                Centroids[i][j]=Centroids[i][j]/Centroids[i][dimension];
            }
            Centroids[i][dimension]=0;
        }
        counter++;
    }

    ifp= fopen(output_filename,"w");
    if (ifp == NULL)
    {
        printf("Error opening the file %s", output_filename);
        exit(-1);
    }

    for(i=0;i<K;i++)/*kmeans*/
    {
        for(j=0;j<dimension;j++)
        {
            if(j==dimension-1)
            {
                fprintf(ifp, "%.4f", Centroids[i][j]);
            }
            else
                fprintf(ifp, "%.4f,", Centroids[i][j]);
        }
        fprintf(ifp,"\n");
    }

    fclose(ifp);
    return 1;
}


