#include <stdio.h>
#include <stdlib.h>
#include "sorting.h"

int choose_pivot(int i,int j )
{
    return((i+j) /2);
}

void quickSort(Struct** a, int l, int r)
{
    int j;
//    printf("l: %d\t r: %d\n",l,r);
    if( l < r )
    {

        //printf("l: %d\t r: %d\n",l,r);
        //printf("\n");
        // divide and conquer
        j = partition( a, l, r);
        quickSort( a, l, j-1);
        quickSort( a, j+1, r);
    }

}

int partition(Struct** a, int l, int r)
{
    int i, j;
//printf("l: %d\t r: %d\n", l,r);
    Struct* pivot;
    pivot = a[l];
    i = l;
    j = r+1;

    while(1)
    {
//        printf("ok0\n");
        do{
            ++i;
            //printf("i:%d\ta[i]->traffic_amount:%f\tpivot->traffic_amount: %f\tr: %d\n",i,a[i]->traffic_amount,pivot->traffic_amount,r);
        }while( a[i]->traffic_amount <= pivot->traffic_amount && i < r );  /**ABBIAMO TOLTO L'UGUALE A I<=R PER FAR FUNZIONARE IL PRIMO QUICKSORT**/

//        printf("ok1\n");
        do
            --j;
        while( a[j]->traffic_amount > pivot->traffic_amount );

//        printf("ok2\n");
//        printf("i: %d\t j: %d\n", i,j);
        if( i >= j )
            break;

//        printf("ok3\n");
//        printf("%f\t%f\n", a[i]->traffic_amount,a[j]->traffic_amount);
//        printf("i: %d\t j: %d\n", i,j);
        swapStruct(a[i],a[j]);
        //printf("%f\t%f\n", a[i]->traffic_amount,a[j]->traffic_amount);
    }
//    printf("l: %d\t j: %d\n", l,j);
    swapStruct(a[l],a[j]);
    return j;
}







//
//void quicksort(Struct** list,int m,int n)
//{
//    int i,j,k;
//    int t=0;
//    double key;
//    if( m < n)
//    {
//        k = choose_pivot(m,n);
//        swapStruct(list[m],list[k]);
//        key = list[m]->traffic_amount;
//        i = m+1;
//        j = n;
//        while(i <= j)
//        {
//            printf("%d\n", t++);
//            while((i <= n) && (list[i]->traffic_amount <= key))
//                i++;
//            while((j >= m) && (list[j]->traffic_amount > key))
//                j--;
//            if( i < j)
//                swapStruct(list[i],list[j]);
//        }
//        /* swap two elements */
//        swapStruct(list[m],list[j]);
//
//        /* recursively sort the lesser list */
//        quicksort(list,m,j-1);
//        quicksort(list,j+1,n);
//    }
//}
//
//void display(double list[],const int n)
//{
//    int i;
//    for(i=0; i<n; i++)
//        printf("%f\t",list[i]);
//}
//
//
//void swap(double *x,double *y)
//{
//    double temp;
//    temp = *x;
//    *x = *y;
//    *y = temp;
//}
//
