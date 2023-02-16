#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "shortest_path.h"

#define INFINITY 4000

double** dijikstra(int G[][16], int N, int startnode, double traffic_matrix[16][16])
{
	int cost[N][N], distance[N], pred[N];
	int visited[N], count, mindistance, nextnode, i,j;
	int destination_node;
	double** f;

	f=malloc(N*sizeof(double));
	for (i=0;i<N;i++)
    {
        f[i]=malloc(N*sizeof(double));
    }

	// initialization of the cost matrix = topology
	for(i=0;i < N;i++)
		for(j=0;j < N;j++)
			if(G[i][j]==0)
				cost[i][j]=INFINITY;
			else
				cost[i][j]=G[i][j];

    //initialization matrix of visited nodes + distance vector of each node from the startnode(=source_node)
	for(i=0;i< N;i++)
	{
		distance[i]=cost[startnode][i];
		pred[i]=startnode;
		visited[i]=0;
	}
	distance[startnode]=0;
	visited[startnode]=1;
	count=1;

	while(count < N-1){
		mindistance=INFINITY;
		for(i=0;i < N;i++)
			if(distance[i] < mindistance && !visited[i])
			{
				mindistance=distance[i];
				nextnode=i;
			}
		visited[nextnode]=1;
		for(i=0;i < N;i++)
			if(!visited[i])
				if(mindistance+cost[nextnode][i] < distance[i])
				{
					distance[i]=mindistance+cost[nextnode][i];
					pred[i]=nextnode;
				}
			count++;
	}

//     output
	for(i=0; i < N;i++)
		if(i!=startnode)
		{
			j=i;
			destination_node=i;
			do
			{
			    f[pred[j]][j]+=traffic_matrix[startnode][destination_node];
				j=pred[j];
			}
			while(j!=startnode);
		}
    return f;
}
