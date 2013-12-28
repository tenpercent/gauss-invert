//#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include "tools.h"
#include "gauss.h"

int main(int argc, char** argv){
  double* a;
  double* b;
  int initID = 0;
  int n=0;
  int m=0;
  int res;
  char* name;
  clock_t t;
  double time=0.;
  FILE* input;
    
  if (argc==3){
    if ( !(n = atoi(argv[1])) || !(m = atoi(argv[2]))){
      printf("This is not a moon!\n");
      return -4;
    }
    //n = atoi(argv[1]);
    //m = atoi(argv[2]);
    a = new double[n*n];
    if (a==NULL){
      printf("Failed to allocate memory for matrix A\n");
      return -2;
    }
    defineMatrixWithFunction(a, n);
    initID=1;
  }
  else if (argc==4){
    if ( !(n = atoi(argv[1])) || !(m = atoi(argv[2]))){
      printf("This is not a moon!\n");
      return -4;
    }
    name = argv[1];
    //n = atoi(argv[2]);
    //m = atoi(argv[3]);
    printf("n=%d m=%d\n", n, m);
    int size = n*n;
    a = new double[size];
    if (a==NULL){
      printf("Failed to allocate memory for matrix A\n");
      return -2;
    }
    input = fopen(name, "r");
    readMatrixFromFile(input, a, n);
    fclose(input);
    initID=2;
  }
  else{
    printf("Usage:\n %s n m\n %s file_name n m\n", argv[0], argv[0]);
    return -1;
  }
  if (initID==0){
    printf("Failed to initialize matrix\n");
    delete[] a;
    return -1;
  }
  b = new double[n*n];
  if (b==NULL){
    printf("Failed to allocate memory for matrix B\n");
    delete[] a;
    return -2;
  }
  makeBlockMatrix(a, n, m);
  printUpperLeftBlock(a, n, m);
  
  idMatrix(b, n);
  makeBlockMatrix(b, n, m);
  
  t = clock();
  res = gaussInvert(a, b, n ,m);
  if (res!=0){
    printf("Method failed!\n");
  }
  time = (clock()-t)*1./CLOCKS_PER_SEC;
  printf("and here we have a record: %.2lf\n", time);
  
  printUpperLeftBlock(b, n, m);

  if (res==0){
    if (initID==2){
      input = fopen(argv[1], "r");
      readMatrixFromFile(input, a, n);
      fclose(input);
    }
    else if (initID==1){
      defineMatrixWithFunction(a, n);
    }
    calculationDeviation(b, a, n, m);
  }
  delete[] a;
  delete[] b;
  return 0;
}
