#ifndef PROVA
#define PROVA 1

#include <stdio.h>
#include <stdlib.h>
#include "struct.h"


void swap(double *x,double *y);
int choose_pivot(int i,int j );
void quickSort(Struct** list,int m,int n);
int partition(Struct** a, int l, int r);
//void display(double list[],const int n);

#endif
