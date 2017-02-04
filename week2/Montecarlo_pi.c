#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char** argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: Montecarlo_pi num_steps\n");
    exit(1);
  }
  
  int num_steps = atoi(argv[1]);
  // Seed the random number generator to get different results each time
  
  int rank;
  int size;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  srand(rank);
  
  int count =0;
  int i;
  float x,y,z;
  for (i=0;i<num_steps;i++){
	  x=drand48();
	  y=drand48();
	  z = pow(x,2)+pow(y,2);
	  if (z < 1.f){
		  count++;
	  }
  }
  
  
	int *conteos = NULL;
	if (rank==0){
		conteos = malloc(size*sizeof(int));
	}
	MPI_Gather(&count, 1, MPI_INT, conteos, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (rank==0){  
		int total_count=0;
		float my_pi;
		for(i=0;i<size;i++){
			total_count+=conteos[i];
		}
		my_pi=total_count/(num_steps*size*0.25);
		
		printf("El numero pi calculado fue %e\n",my_pi);
	}
  MPI_Finalize();
}
