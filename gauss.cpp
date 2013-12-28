#include "gauss.h"
#include "tools.h"

int gaussInvert(double* const a, double* const b, int n, int m){//n - размер матрицы, m - размер блока
  //предполагаем, что в матрице А нумерация элементов уже "хитрая"
  //константы
  if ((n<1)||(m<1)){
    printf("bad matrix size or block size");
    return -3;
  }
  const int s = n/m;
  const int r = (n-s*m);
  //int matrix_size = n*n;
  const int string_size = n*m;
  const int block_size = m*m;
  const int small_block_size = r*m;
  const int smallest_block_size = r*r;

  int i, j, k;

  int labelMain;
  double temp=0.;
  double minNorm=0.;
  int res=0;

  //double* temp_block_pointer;

  //double* pointer1;
  //double* pointer2;

  int p=0;
  int q=0;//A_pq - блоки
  //double* c;//надо как-то будет заново прочитать исходную матрицу, но это всё потом
  //printf("Let's begin!\n");
  //res = idMatrix(b, n);
  //res = makeBlockMatrix(b, n, m);//а вообще, оно может работать невменяемо долго
  //ладно, единичную матрицу как-нибудь завели.

  double* const temp_block = new double[block_size];//не забыть потом фрей наставить
  double* const temp_block_storage = new double[block_size];
  double* const blocks_order = new double[s];

  if (temp_block==NULL){
    printf("not enough memory for temp_block;\n");
    return -2;
  }
  if (temp_block_storage==NULL){//нужен только для главного блока
    printf("not enough memory for temp_block;\n");
    return -2;
  }
  if (blocks_order==NULL){
    printf("not enough memory for temp_block;\n");
    return -2;
  }

  if (r!=0){
    for (i=0; i<s; i++){
      labelMain=0;//главного элемента не найдено - метка
      p=0; 
      q=0;
      temp = 0.;
      minNorm = 1e+100;//вроде хватит
      //zeroMatrix(temp_block_storage, block_size);
      for (j=i; j<s; j++){
	for (k=i; k<s; k++){
	  res = simpleInvert(a + j*string_size + k*block_size, temp_block, m);
	  if (res==0) {
	    temp=matrixNorm(temp_block, m);
	    if (temp<minNorm){
	      minNorm = temp;
	      p=j;
	      q=k;
	      //и матрицу надо запомнить
	      res = copyMatrix(temp_block, temp_block_storage, block_size);
	      labelMain=1;
	    }
	  }
	}
      }
      
      if (labelMain){
	for (j=0; j<s; j++){
	  res = swapMatrix(a + j*string_size + i*block_size, a + j*string_size + q*block_size, block_size);
	}
	res = swapMatrix(a + s*string_size + i*small_block_size, a + s*string_size + q*small_block_size, small_block_size);
	for (j=i; j<s; j++){
	  res = swapMatrix(a + p*string_size + j*block_size, a + i*string_size + j*block_size, block_size);
	}
	res = swapMatrix(a + p*string_size + s*block_size, a + i*string_size + s*block_size, small_block_size);
	for (j=0; j<s; j++){
	  res = swapMatrix(b + p*string_size + j*block_size, b + i*string_size + j*block_size, block_size);
	}
	res = swapMatrix(b + p*string_size + s*block_size, b + i*string_size + s*block_size, small_block_size);  
	blocks_order[i] = q;
	
	idMatrix(a+i*string_size+i*block_size, m);
	for (j=i+1; j<s; j++){
	  res = simpleMatrixMultiply(temp_block_storage, a + i*string_size+j*block_size, temp_block, m, m, m);
	  res = copyMatrix (temp_block, a + i*string_size+j*block_size, block_size);
	}
	res = simpleMatrixMultiply(temp_block_storage, a + i*string_size+s*block_size, temp_block, m, m, r);
	res = copyMatrix(temp_block, a + i*string_size+s*block_size, small_block_size);
	for (j=0; j<s; j++){
	  res = simpleMatrixMultiply(temp_block_storage, b + i*string_size+j*block_size, temp_block, m, m, m);
	  res = copyMatrix(temp_block, b + i*string_size + j*block_size, block_size);
	}
	res = simpleMatrixMultiply(temp_block_storage, b + i*string_size+s*block_size, temp_block, m, m, r);
	res = copyMatrix(temp_block, b + i*string_size+s*block_size, small_block_size);
		  
	for (j=i+1; j<s; j++){//j-я блочная строчка 
	  for (k=i+1; k<s; k++){//k-й блочный столбец
	    //можно переписать с pointer1 и pointer2
	    res = simpleMatrixMultiply(a+j*string_size+i*block_size, a+i*string_size+k*block_size, temp_block, m, m, m);
	    res = subtractFromMatrix(a+j*string_size+k*block_size, temp_block, block_size);
	  }
	  res = simpleMatrixMultiply(a+j*string_size+i*block_size, a+i*string_size+s*block_size, temp_block, m, m, r);
	  res = subtractFromMatrix(a+j*string_size+s*block_size, temp_block, small_block_size);
	}
	for (j=i+1; j<s; j++){
	  for (k=0; k<s; k++){
	    res = simpleMatrixMultiply(a+j*string_size+i*block_size, b+i*string_size+k*block_size, temp_block, m, m, m);
	    res = subtractFromMatrix(b+j*string_size+k*block_size, temp_block, block_size);
	  }
	  res = simpleMatrixMultiply(a+j*string_size+i*block_size, b+i*string_size+s*block_size, temp_block, m, m, r);
	  res = subtractFromMatrix(b+j*string_size+s*block_size, temp_block, small_block_size);
	}
	for (k = i+1; k<s; k++){
	  res = simpleMatrixMultiply(a+s*string_size+i*small_block_size, a+i*string_size+k*block_size, temp_block, r, m, m);
	  res = subtractFromMatrix(a+s*string_size+k*small_block_size, temp_block, small_block_size);
	}
	res = simpleMatrixMultiply(a+s*string_size+i*small_block_size, a+i*string_size+s*block_size, temp_block, r, m, r);
	res = subtractFromMatrix(a+s*string_size+s*small_block_size, temp_block, smallest_block_size);
	for (k = 0; k<s; k++){
	  res = simpleMatrixMultiply(a+s*string_size+i*small_block_size, b+i*string_size+k*block_size, temp_block, r, m, m);
	  res = subtractFromMatrix(b+s*string_size+k*small_block_size, temp_block, small_block_size);
	}
	res = simpleMatrixMultiply(a+s*string_size+i*small_block_size, b+i*string_size+s*block_size, temp_block, r, m, r);
	res = subtractFromMatrix(b+s*string_size+s*small_block_size, temp_block, smallest_block_size);
	for (j=i+1; j<s; j++){
	  res = zeroMatrix(a+j*string_size+i*block_size, block_size);
	} 
	res = zeroMatrix(a+s*string_size+i*small_block_size, small_block_size);
      }
      else{
	printf("This method failed. Just as planned\n\t --gaussInvert\n");
	delete[] temp_block;
	delete[] temp_block_storage;
	delete[] blocks_order;
	return -1;
      }
    }
    res = simpleInvert(a+s*string_size+s*small_block_size, temp_block, r);
    if (res==0){
      for (i=0; i<s; i++){
	res = simpleMatrixMultiply(temp_block, b+s*string_size+i*small_block_size, temp_block_storage, r, r, m);
	res = copyMatrix(temp_block_storage, b+s*string_size+i*small_block_size, small_block_size);
      }
      res = simpleMatrixMultiply(temp_block, b+s*string_size+s*small_block_size, temp_block_storage, r, r, r);
      res = copyMatrix(temp_block_storage, b+s*string_size+s*small_block_size, smallest_block_size);
      //just in case
      //res = simpleMatrixMultiply(temp_block, a+s*string_size+s*small_block_size, temp_block_storage, r, r, r);
      //res = copyMatrix(temp_block_storage, a+s*string_size+s*small_block_size, smallest_block_size);
      //тоже не забыть заменить этот кусок
      idMatrix(a+s*string_size+s*small_block_size, r);
    }
    else{
      printf("The method failed. So close...\n");
      delete[] temp_block;
      delete[] temp_block_storage;
      delete[] blocks_order;
      return -1;
    }
    for (i=s-1; i>=0; i--){
      res = simpleMatrixMultiply(a+i*string_size+s*block_size, b+s*string_size+s*small_block_size, temp_block, m, r, r);
      res = subtractFromMatrix(b+i*string_size+s*block_size, temp_block, small_block_size);
      for (k=s-1; k>i; k--){
	res = simpleMatrixMultiply(a+i*string_size+k*block_size, b+k*string_size+s*block_size, temp_block, m, m, r);
	res = subtractFromMatrix(b+i*string_size+s*block_size, temp_block, small_block_size);
      }
      for (j=0; j<s; j++){
	res = simpleMatrixMultiply(a+i*string_size+s*block_size, b+s*string_size+j*small_block_size, temp_block, m, r, m);
	res = subtractFromMatrix(b+i*string_size+j*block_size, temp_block, block_size);
	for (k=s-1; k>i; k--){
	  res = simpleMatrixMultiply(a+i*string_size+k*block_size, b+k*string_size+j*block_size, temp_block, m, m, m);
	  res = subtractFromMatrix(b+i*string_size+j*block_size, temp_block, block_size);
	}
      }
    }
    for (i=s-1; i>=0; i--){
      q = blocks_order[i];
      for(j=0; j<s; j++){
	swapMatrix(b+i*string_size+j*block_size, b+q*string_size+j*block_size, block_size);
      }
      swapMatrix(b+i*string_size+s*block_size, b+q*string_size+s*block_size, small_block_size);
    }
  }
  else if (r==0){
    printf("String size: %d\n", string_size);
    printf("Block size: %d\n", block_size);
    for (i=0; i<s; i++){
      labelMain=0;//главного элемента не найдено - метка
      p=0; 
      q=0;
      temp = 0.;
      minNorm = 1e+111;//вроде хватит
      for (j=i; j<s; j++){
	for(k=i; k<s; k++){
	  res = simpleInvert(a+j*string_size+k*block_size, temp_block, m);
	  if (res==0){
	    temp = matrixNorm(temp_block, m);
	    if (temp<minNorm){
	      minNorm = temp;
	      p=j;
	      q=k;
	      copyMatrix(temp_block, temp_block_storage, block_size);
	      labelMain = 1;
	    }
	  }
	}
      }
      //labelMain=1;
      if (labelMain){
	//p=i; q=i;
	//res = simpleInvert(a+i*string_size+i*block_size, temp_block_storage, m);
        for (j=i; j<s; j++){
	  swapMatrix(a+p*string_size+j*block_size, a+i*string_size+j*block_size, block_size); 
	}//свап строчек в левой матрице
	for (j=0; j<s; j++){
	  swapMatrix(b+p*string_size+j*block_size, b+i*string_size+j*block_size, block_size); 
	}
	for (j=0; j<s; j++){
	  swapMatrix(a+j*string_size+i*block_size, a+j*string_size+q*block_size, block_size); 
	}
	blocks_order[i]=q;

	idMatrix(a+i*string_size+i*block_size, m);//это крайне неочевидно, но генерация единичной матрицы
	//работает медленнее, чем умножение
	for (j=i+1; j<s; j++){
	  res = simpleMatrixMultiply(temp_block_storage, a + i*string_size+j*block_size, temp_block, m, m, m);
	  res = copyMatrix (temp_block, a + i*string_size+j*block_size, block_size);
	}
	for (j=0; j<s; j++){
	  res = simpleMatrixMultiply(temp_block_storage, b + i*string_size+j*block_size, temp_block, m, m, m);
	  res = copyMatrix(temp_block, b + i*string_size + j*block_size, block_size);
	} 

	for (j=i+1; j<s; j++){//j-я блочная строчка 
	  for (k=i+1; k<s; k++){//k-й блочный столбец
	    res = simpleMatrixMultiply(a+j*string_size+i*block_size, a+i*string_size+k*block_size, temp_block, m, m, m);
	    res = subtractFromMatrix(a+j*string_size+k*block_size, temp_block, block_size);
	  }
	}
	for (j=i+1; j<s; j++){
	  for (k=0; k<s; k++){
	    res = simpleMatrixMultiply(a+j*string_size+i*block_size, b+i*string_size+k*block_size, temp_block, m, m, m);
	    res = subtractFromMatrix(b+j*string_size+k*block_size, temp_block, block_size);
	  }
	}
	/*
	printf("After subtraction:\n");
	printf("left:\n");
	printBlockMatrix(a, n, m);
	printf("right:\n");
	printBlockMatrix(b, n, m);
	*/
	for (j=i+1; j<s; j++){
	  res = zeroMatrix(a+j*string_size+i*block_size, block_size);
	}	
      }
      else{
	printf("This method failed. Just as planned\n\t --gaussInvert\n");
	delete[] temp_block;
	delete[] temp_block_storage;
	delete[] blocks_order;
	return -1;
      }	
    }
    for (i=s-1; i>=0; i--){
      for (j=0; j<s; j++){
	for (k=s-1; k>i; k--){
	  res = simpleMatrixMultiply(a+i*string_size+k*block_size, b+k*string_size+j*block_size, temp_block, m, m, m);
	  res = subtractFromMatrix(b+i*string_size+j*block_size, temp_block, block_size);
	}
      }
    }
    for(i=s-1; i>=0; i--){
      for (j=0; j<s; j++){
	q=blocks_order[i];
	swapMatrix(b+i*string_size+j*block_size, b+q*string_size+j*block_size, block_size);
      }
    }
  }
  else{
    printf("Something weird is going on...\n");
  }
  delete[] temp_block;
  delete[] temp_block_storage;
  delete[] blocks_order;
  return 0;
}
