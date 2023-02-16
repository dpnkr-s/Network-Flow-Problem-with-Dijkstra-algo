#include <stdio.h>
#include <stdlib.h>
#include "random.h"
#include "sorting.h"
#include "struct.h"
#include "shortest_path.h"

#define N 16
#define A 0.5 // lower bound of the traffic flows
#define B 1.5 // upper bound of the traffic flows

int check_possible_links();

/**
generated_traffic[i]=sum t(i,j) for each i,j
sort generated_traffic[i]                              //the quicksort put values in increasing order!
do
    if l'i-th allows another link outcoming from it
        take k nodes toward it sends more traffic
        select one of them and add it to the topology + subtract the traffic of the added link from generated_traffic[i] + add the updated generated_traffic[i]
        sort
    otherwhise try with con (i+1)-th
while i can add a link
**/

int delta=4;
int delta_tmp_inflow[N]= {};
int delta_tmp_outflow[N]= {};

int main()
{
    long seed=4;
    int i,j,t;
    double traffic_matrix[N][N];
    Struct** traffic_generated_structure;
    int b[N][N]= {};
    int source_index,destination_index;
    double f[N][N]= {};
    double** f_tmp;
    double fmax=0;
    FILE* fp=fopen("Document_Our_Greedy_Algorithm.txt","a");
    traffic_generated_structure=malloc(N*sizeof(Struct*));

    // initialization
    //    - of the traffic matrix with flow between A and B --> traffic_matrix
    //    - of the total traffic injected by each node vector --> traffic_generated_vector
    for (i=0; i<N; i++)
    {
        traffic_generated_structure[i] = newStruct(0.0,i);
        for (j=0; j<N; j++)
        {
            if (i!=j)
            {
                traffic_matrix[i][j] = uniform((float)A,(float)B,&seed);
                addTrafficNode(traffic_generated_structure[i],traffic_matrix[i][j]);
            }
            else
            {
                // in order to avoid loops, we consider the traffic toward the node itself equal to a negative value
                traffic_matrix[i][j]=-1;
            }
        }
    }

//    displayStruct(traffic_generated_structure,N);
//    printf("\n");
    quickSort(traffic_generated_structure,0,N-1);
//    displayStruct(traffic_generated_structure,N);
//    printf("fine primo quicksort\n");


    do
    {
//        printf("choose_source function start\n");
        // find first node which admits a new out-link
        source_index = choose_source(traffic_generated_structure,N);
//        printf("choose_source function end\n");
//        printf("source_index: %d",source_index);
        if (source_index==-1)
        {
            printf("There are no nodes which admit an out-link\n");
            break;
        }
        // choose the destination
        destination_index = choose_destination(traffic_matrix[source_index],N,source_index,b);
        if (destination_index==-1)
        {
            printf("There are no nodes which admit an in-link from that destination\n");
            break;
        }
        // add link
        if (source_index!=destination_index )
        {
            b[source_index][destination_index]++;
            traffic_generated_structure[source_index]->traffic_amount-=traffic_matrix[source_index][destination_index];
//            printf("\ndestination:%d\t source:%d\n",delta_tmp_inflow[destination_index],delta_tmp_outflow[source_index]);
            delta_tmp_inflow[destination_index]++;
            delta_tmp_outflow[source_index]++;
//            printf("destination:%d\t source:%d\n",delta_tmp_inflow[destination_index],delta_tmp_outflow[source_index]);
        }
        else
        {
            break;
        }

        quickSort(traffic_generated_structure,0,N-1);
        // printf("Fine quicksort n: %d\n",t++);

    }
    while(check_possible_links());


    for(j=0; j<N; j++)
    {
        printf("%d ",delta_tmp_inflow[j]);
        printf("%d ",delta_tmp_outflow[j]);
        printf("\n");
    }
    // computing all f_tmp matrices
    for (t=0; t<N; t++)
    {

        f_tmp=dijikstra(b,N,t,traffic_matrix);

        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                f[i][j]+=f_tmp[i][j];
                fprintf(fp,"%.2f/%.2f \t",f_tmp[i][j],f[i][j]);
            }
            fprintf(fp,"\n");

        }
//        fprintf(fp,"***************************************************************************************\n");

    }

    // searching the fmax
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
            if(fmax<f[i][j])
                fmax=f[i][j];
    }

    // output
//    for(i=0; i<N; i++)
//    {
//        for(j=0; j<N; j++)
//            fprintf(fp,"%.2f \t",f[i][j]);
//        fprintf(fp,"\n\n");
//    }
//    fprintf(fp,"******************************************************************************************************");
    fprintf(fp,"Number of nodes: %d\n",N);
    fprintf(fp,"Uniform distribution of traffic between %f and %f \n", (float)A,(float)B);
    fprintf(fp,"Delta: %d\n",delta);
    fprintf(fp,"Fmax: %f\n",fmax);
    fprintf(fp,"******************************************************************************************************\n");
    printf("%f\n",fmax);
    return 0;
}



int check_possible_links()
{

    int i;
    int possible_in=0, possible_out=0;

    for (i=0; i<N; i++)
    {
        if (delta_tmp_inflow[i]<delta)
            possible_in++;
    }

    for (i=0; i<N; i++)
    {
        if (delta_tmp_outflow[i]<delta)
            possible_out++;
    }

    if (possible_in>0 && possible_out>0)
        return 1;

    return 0; //false
}
