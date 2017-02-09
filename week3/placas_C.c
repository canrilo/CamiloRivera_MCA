#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int L = 5, l = 2, d = 1, V0 = 100, m, N, num, rank, size, r;
double h = 0.01953125;

int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);

int main(void)
{
	int up, down, left, right, x0, x1, y0, y1, i=1, j=1, n=0, ii;
	double average;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	MPI_Request send_request,send_request2, recv_request,recv_request2;
	MPI_Status status;
	
	if(!ispowerof2(size))
	{
		printf("Number of processors must be power of 2\n");
		exit(0);
	}
	
	m = 256;   //L/h;
	N = 2*m*m;

	x0 = (int) m/2 - l/(h*2) - 1;
	x1 = (int) m/2 + l/(h*2) - 1;
	y0 = (int) m/2 - d/(h*2) - 1;
	y1 = (int) m/2 + d/(h*2) - 1;
	
	num = m/size + 2; //Number of rows per processor
	if (rank == 0 || rank == (size-1))
	{
		num--;
	}
	r = !(rank==0); 
	//r = 1; 			//Number for index adjustment if rank 0
	//if (rank==0)
	//{
	//	r=0;
	//}
	
	double *V = malloc((m*num)*sizeof(double));
	double *V_new = malloc(m*num*sizeof(double));
	
	V = init(x0, x1, y0, y1, V);
	V_new = init(x0, x1, y0, y1, V);
	
	while (n < N)
	{		
		for(i=1;i < num-1; i++)
		{
			ii = i + rank*m/size - r;
			for(j=1;j < m-1; j++)
			{
				up = transformer(i-1, j);
				down = transformer(i+1, j);
				left = transformer(i, j-1);
				right = transformer(i, j+1);
				if (!(j >= x0 && j <= x1 && ii == y0) && !(j >= x0 && j <= x1 && ii == y1))
				{
					average = (V[up] + V[down] + V[left] + V[right])/4;
					V_new[transformer(i,j)] = average;
				}
			}
		}
		for(i=1;i < num-1; i++)
	    {
			for(j=1;j < m-1; j++)
		  {
			V[transformer(i,j)] = V_new[transformer(i,j)] ;
		  }
		}
		if(!rank==0)
		{
			MPI_Irecv(&V[m], m, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &recv_request);
			MPI_Isend(&V[0], m, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &send_request);
		}
		else if (rank !=(size-1))
		{
			MPI_Irecv(&V[m*(num-2)], m, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_request);
			MPI_Isend(&V[m*(num-1)], m, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &send_request);
		}
		MPI_Wait(&send_request, &status);
		MPI_Wait(&recv_request, &status);
		MPI_Wait(&send_request2, &status);
		MPI_Wait(&recv_request2, &status);
		n += 1;
	}
	free(V_new);
	
	double *V_TOTAL,*V_Reduced;
	V_Reduced = malloc(m*m/size*sizeof(double));
	if (rank==0)
	{
		V_Reduced = V;
	}
	else
	{
		V_Reduced = &V[m];
	}
	free(V);
	if(rank==0)
	{
		V_TOTAL=malloc(m*m*sizeof(double));
	}
	MPI_Gather(V_Reduced,(int)m*m/size, MPI_DOUBLE, V_TOTAL, (int)m*m/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank==0)
	{	
		for(i=0;i < num; i++)
		{
			for(j=0;j < m; j++)
			{
				printf("%f\n", V_TOTAL[transformer(i,j)]);
			}
		}
	}
	MPI_Finalize();
	return 0;
}

int transformer(int i, int j)
{	
	return i*m + j;
}

double *init(int x0, int x1, int y0, int y1, double *array)
{	
	int ii;
	for(i=0;i < num; i++)
	{
		ii = i + rank*m/size - r;
		for(j=0;j < m; j++)
		{
			
			if (ii==y0 && j>=x0 && j<= x1)
			{
				array[transformer(i, j)] = V0/2;
			}
			else if (ii==y1 && j>=x0 && j<= x1)
			{
				array[transformer(i, j)] = -V0/2;
			}
			else
			{
				array[transformer(i,j)] = 0;
			}
		}
	}
}


int ispowerof2(unsigned int x) {
   return x && !(x & (x - 1));
 }
