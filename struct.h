#ifndef STRUCT_H_INCLUDED
#define STRUCT_H_INCLUDED


typedef struct {
    double traffic_amount;
    int node;
} Struct;

void swapStruct(Struct* struct1, Struct* struct2);
Struct* newStruct(double traffic_amount, int node);
double getAmount(Struct* structure);
int getNode(Struct* structure);
void addTrafficNode(Struct *structure,double amount);
void displayStruct(Struct** structure,int N);
void displaySingleStruct(Struct* structure);
int choose_source(Struct** traffic_generated_structure,int N);
int choose_destination(double *traffic_row,int N,int source_index,int b[][16]);

#endif // STRUCT_H_INCLUDED
