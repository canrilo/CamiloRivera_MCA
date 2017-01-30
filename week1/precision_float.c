#include <stdio.h>
#include <stdlib.h>


int main(void){
  unsigned int id;
  int min;
  float x;
  float x1;
  int i;
  float *interv;
  FILE *exp;
  char nombrefof[30] = "Precisiones_float.txt";
  
  exp=fopen(nombrefof,"w");
	if(!exp){
	printf("Hubo problemas abriendo el archivo %s\n",nombrefof);
	}
  
  interv=malloc(24*sizeof(float));

  x = 1E-12;
  
  for (i=0;i<25;i++){
	id = (*(unsigned int*)&x);
	id++;
    x1 = (*(float*)&id);
    //printf("%.25e\n", x1-x);
    interv[i]=x1-x;
    fprintf(exp,"%e\n",interv[i]);
    x=x*10;
  }
  
  fclose(exp);

  return 0;
}
