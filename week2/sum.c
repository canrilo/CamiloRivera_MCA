// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// MPI_Send, MPI_Recv example. Communicates the number -1 from process 0
// to process 1.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int number;
  int number1;
  int i;
  number = rand()*world_rank % 20;
  printf("My name is %d and my number is %d\n",world_rank,number);
  if (world_rank > 0) {

    MPI_Send(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

  }
  else if (world_rank ==0){
    // If we are rank 0, set the number to -1 and send it to process 1
    int sum = number;
    for(i=1;i<world_size;i++){
        MPI_Recv(&number1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        sum+=number1;
    }
    printf("The sum of the numbers is %d\n",sum);

  }
  MPI_Finalize();
}
