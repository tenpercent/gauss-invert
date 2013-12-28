#include "tools.h"

int simpleInvert(double* const a, double* const b, int m){
  //ищем обратный к блоку размера m*m
  double eps = 1e-8;
  if ( !m ) {
    return -1;
  }

  const int size = m * m;
  double temp;
  int point1 = 0;
  int point2 = 0;
  int i, j, k;
  int p = 0;
  double tryMain = 0.;

  double* const c = new double[size];//бэкап данной матрицы
  idMatrix(b, m);
  copyMatrix(a, c, size);//забекапились

  for (i = 0; i < m; i++) {
    temp = 0.;
    tryMain = 0.;
    p = -1;
    for (j = i; j < m; j++){
      tryMain = fabs(a[j * m + i]);
      if (tryMain > temp + eps){
      	temp = tryMain;
      	p = j;
      }
    }
    if(p > -1){
      temp = 0.;
      for (j = i; j < m; j++){
      	point1 = i * m + j;
      	point2 = p * m + j;

      	temp      = a[point1];
      	a[point1] = a[point2];
      	a[point2] = temp;
      }
      for (j = 0; j < m; j++){
      	point1 = i * m + j;
      	point2 = p * m + j;

      	temp      = b[point1];
      	b[point1] = b[point2];
      	b[point2] = temp;
      }
    }
    else{
      //мы не нашли главного элемента, a это значит, что матрица вырождена
      copyMatrix(c, a, size);
      delete[] c;
      return -1;
    }
    
    temp = 1. / a[i * m + i];

    for (j = i + 1; j < m; j++){
      a[i * m + j] *= temp;
    }
    for (j = 0; j < m; j++){
      b[i * m + j] *= temp;
    }
    a[i * m + i] = 1.;
    //теперь в левом верхнем углу - единичный элемент
    for (j = i + 1; j < m; j++){//из j-й строчки вычитаем i-ю      
      temp = a[j * m + i];
      for (k = 0; k < m; k++){
      	b[j * m + k] -= temp * b[i * m + k];
      }
      for (k = i + 1; k < m; k++){
      	a[j * m + k] -= temp * a[i * m + k];
      }
      a[j * m + i] = 0.;
    }
  }//конец прямого хода алгоритма Гаусса

  for (i = m - 1; i > 0; i--){//строчку с номером i вычитаем из
    for (j = i - 1; j >= 0; j--){//первых (i-1) строчек
      for(k = 0; k < m; k++){
      	b[j * m + k] -= b[i * m + k] * a[j * m + i];
      }
    }
  }
  copyMatrix(c, a, size);
  delete[] c;
  return 0;
}
