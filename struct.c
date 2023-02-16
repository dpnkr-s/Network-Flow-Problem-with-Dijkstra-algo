#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "sorting.h"

#define K 4

extern int delta;
extern int delta_tmp_outflow[];
extern int delta_tmp_inflow[];


void swapStruct(Struct* struct1, Struct* struct2)
{
    Struct* tmp=malloc(sizeof(Struct));
    *tmp=*struct1;
    *struct1=*struct2;
    *struct2=*tmp;

}

void displaySingleStruct(Struct* structure)
{
    printf("Structure amount: %f \t node: %d \n", getAmount(structure),(getNode(structure)+1));
    return;
}

Struct* newStruct(double traffic_amount, int node)
{

    Struct* tmp = malloc(sizeof(Struct));
    tmp->node=node;
    tmp->traffic_amount=traffic_amount;

    return tmp;
}

int choose_source(Struct** traffic_generated_structure, int N)
{
    int index=-1;
    int found=0,i;
    int possible_index;

    for (i=N-1; i>=0 && found==0; i--)
    {
        possible_index=traffic_generated_structure[i]->node;
        if (delta_tmp_outflow[possible_index]<delta)
        {
            index=possible_index;
            found=1;
        }
    }

    return index;
}

int choose_destination(double *traffic_row,int N,int source_index,int b[][16])
{
    int index=-1;
    int i,t;
    int possible_index;
    int possible_links=0;
    Struct* row_struct[N];
    int candidate_nodes[K];

    for (i=0; i<N; i++)
    {
        row_struct[i] = newStruct(traffic_row[i],i);
    }

    quickSort(row_struct,0,N-1);

    for (t=N-1; t>0 && possible_links<K; t--)
    {
        if(delta_tmp_inflow[row_struct[t]->node]<delta)
        {
            candidate_nodes[possible_links]=row_struct[t]->node;
            possible_links++;
        }
    }
    if (possible_links==0)
        return -1;

    possible_index=rand()%possible_links;
    index=candidate_nodes[possible_index];

    for(i=0; i<N; i++)
        free(row_struct[i]);
    free(row_struct);

    return index;
}


double getAmount(Struct* structure)
{
    return structure->traffic_amount;
}


int getNode(Struct* structure)
{
    return structure->node;
}

void addTrafficNode(Struct *structure,double amount)
{
    structure->traffic_amount+=amount;
    return;
}

void displayStruct(Struct** structure,int N)
{

    int i;
    for (i=0; i<N; i++)
        printf("Amount: %f \t Node:%d\n",structure[i]->traffic_amount,(structure[i]->node)+1);
    return;
}



//void display(double list[],const int n)
//{
//    int i;
//    for(i=0; i<n; i++)
//        printf("%f\t",list[i]);
//}
//
