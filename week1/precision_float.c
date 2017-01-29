#include <stdio.h>
#include <stdlib.h>


int main(void){
  unsigned int id;
  int min;
  float x;
  float x1;
  int i;
  float *interv;
  
  interv=malloc(24*sizeof(float));

  x = 1E-12;
  
  for (i=0;i<25;i++){
	id = (*(unsigned int*)&x);
	id++;
    x1 = (*(float*)&id);
    printf("%.25e\n", x1-x);
    interv[i]=x1-x;
    x=x*10;
  }

  return 0;
}
