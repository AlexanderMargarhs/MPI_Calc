#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int main (int argc, char **argv) 
{
    MPI_Status status;
    MPI_Init(&argc, &argv);
	// Initialization. Start
    int Sum = 0, Process, My_Rank, n, Load, i, j, k, Ranking, Leftover, Answer = 1, *Max, Maximum, Position;
	int *A, *B, *C, *D, *Sum_C_D, *Prod_B_C, *Prod_A_B, *Prod_C_D, *A_loc, *B_loc, *C_loc, *D_loc, *Sum_C_D_loc, *Prod_B_C_loc, *Prod_A_B_loc, *Prod_C_D_loc, *Temp;
	
    MPI_Comm_size(MPI_COMM_WORLD, &Process); 
    MPI_Comm_rank(MPI_COMM_WORLD, &My_Rank); 

	while(Answer == 1)
	{
		if(My_Rank == 0)
		{
			printf("Please type the size of the array.\n"); // Get the data of the matrices and initialize the arrays for process 0.
			scanf("%d", &n);
			Load = n / Process;
			Leftover = n % Process;
			if(Leftover == 0)
			{
				A = (int *)malloc(sizeof(int) * n);
				D = (int *)malloc(sizeof(int) * n * n);
				Prod_A_B = (int *)malloc(sizeof(int) * n);
				Prod_C_D = (int *)malloc(sizeof(int) * n * n);
				Sum_C_D = (int *)malloc(sizeof(int) * n * n);
			}
			Max = (int *)malloc(sizeof(int) * Load + Leftover);
			C = (int *)malloc(sizeof(int) * n * n);
			Prod_B_C = (int *)malloc(sizeof(int) * n);
		}
		
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the data that the other processes are going to use.
		MPI_Bcast(&Load, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&Leftover, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		B = (int *)malloc(sizeof(int) * n); // Initialize the matrices that the processes are going to use.
		if(Leftover == 0)
		{
			A_loc = (int *)malloc(sizeof(int) * Load);
			B_loc = (int *)malloc(sizeof(int) * Load);
			C_loc = (int *)malloc(sizeof(int) * n * Load );
			D_loc = (int *)malloc(sizeof(int) * n * Load );
			Temp = (int *)malloc(sizeof(int) * n);		
			Prod_B_C_loc = (int *)calloc(sizeof(int), Load);
			Prod_A_B_loc = (int *)calloc(sizeof(int), Load);
			Prod_C_D_loc = (int *)calloc(sizeof(int), n);
		}		
		else
		{
			C_loc = (int *)malloc(sizeof(int) * n * Load);
			Prod_B_C_loc = (int *)calloc(sizeof(int), Load);
		}
		
		Sum_C_D_loc = (int *)malloc(sizeof(int) * n * Load);

		if(My_Rank == 0)
		{
			if(Leftover == 0) // Check if the n % p is not 0.
			{
				printf("Type the data for the matrix A.\n");	
				for(i = 0; i < n; i ++)
				{
					scanf("%d", &A[i]);
				}
				printf("Type the data for the matrix B.\n");	
				for(i = 0; i < n; i ++)
				{
					scanf("%d", &B[i]);
				}
				printf("Type the data for the matrix C.\n");	
				for(i = 0; i < n * n; i ++)
				{		
					scanf("%d", &C[i]);
				}
				printf("Type the data for the matrix D.\n");	
				for(i = 0; i < n * n; i ++)
				{		
					scanf("%d", &D[i]);
				}	
			}
			else
			{
				printf("Type the data for the matrix B.\n");	
				for(i = 0; i < n; i ++)
				{
					scanf("%d", &B[i]);
				}
				printf("Type the data for the matrix C.\n");	
				for(i = 0; i < n * n; i ++)
				{		
					scanf("%d", &C[i]);
				}
				printf("The matrices A and D are not needed.\n");
			}
		}
		// Initialization. End

		if(Leftover == 0)
		{
			// I. Start
			MPI_Scatter(C, n * Load, MPI_INT, C_loc, n * Load, MPI_INT, 0, MPI_COMM_WORLD); // Send the data of C.
			MPI_Scatter(D, n * Load, MPI_INT, D_loc, n * Load, MPI_INT, 0, MPI_COMM_WORLD); // Send the data of D.
			
			for(i = 0; i < Load * n; i ++)
			{
				Sum_C_D_loc[i] = C_loc[i] + D_loc[i]; // Addition.
			}
			
			MPI_Gather(Sum_C_D_loc, n * Load, MPI_INT, Sum_C_D, n * Load, MPI_INT, 0, MPI_COMM_WORLD); // Gather the data to the Sum_C_D array.
			// I. End
			
			// II. Start

			MPI_Bcast(B, n, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the data of B.
					
			for(i = 0; i < Load; i ++)
			{
				for(j = 0; j < n; j ++)
				{
					Prod_B_C_loc[i] += B[j] * C_loc[i * n + j]; // Find the product of the row C and the array B.
				}	
			}
			MPI_Reduce(Prod_B_C_loc, Max, Load, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

			MPI_Gather(Prod_B_C_loc, Load, MPI_INT, Prod_B_C, Load, MPI_INT, 0, MPI_COMM_WORLD); // Gather the data to the array.
			if(My_Rank != 0)
			{
				free(B); // Free the B array.
			}
			// II. End

			// III. Start 
			MPI_Scatter(A, Load, MPI_INT, A_loc, Load, MPI_INT, 0, MPI_COMM_WORLD); // Scatter the data of A.
			MPI_Scatter(B, Load, MPI_INT, B_loc, Load, MPI_INT, 0, MPI_COMM_WORLD); // Scatter the data of B.

			for(i = 0; i < Load; i ++) // Sum_C_D of the products of A and B matrixs.
			{
				Prod_A_B_loc[i] += (A_loc[i] * B_loc[i]); // Find the Sum of products of the arrays A and B.
			}

			MPI_Gather(Prod_A_B_loc, Load, MPI_INT, Prod_A_B, Load, MPI_INT, 0, MPI_COMM_WORLD); // Process 0 gathers the Sum of products from the other processes.  
			
			if(My_Rank == 0) 
			{
				for(i = 0; i < n; i ++ )
				{
					Sum += Prod_A_B[i]; //Find the final Sum.
				}
			}
			// III. End

			// IV. Start
			if(n == Process) // Check if we can do the calculations.
			{
				if (Process == 1)
				{
					Prod_C_D_loc[0] = C_loc[0] * D_loc[0]; 
				}
				else
				{
					for(Ranking = My_Rank, i = 0; i < Process; i ++, Ranking ++)
					{
						if (Ranking >= n)
						{
							Ranking = 0;
						}
						for(j = 0; j < n; j ++)
						{
							Prod_C_D_loc[j] += (C_loc[Ranking] * D_loc[j]); // Calculate it to the matrix.
						}

						if(My_Rank == 0) 
						{                    
							for(k = 0; k < n; k ++)
							{
								Temp[k] = D_loc[k]; // Save the matrix so we won't lose the data.
							}          
							MPI_Send(&Temp[0], n, MPI_INT, Process - 1, 1, MPI_COMM_WORLD); // Send to the next process and receive data from the last process. 
							MPI_Recv(&D_loc[0], n, MPI_INT, (My_Rank + 1) % n, 1, MPI_COMM_WORLD, &status);
						} 
						else 
						{              
							for(k = 0; k < n; k ++)
							{
								Temp[k] = D_loc[k]; // Save the matrix so we won't lose the data.
							}                                  
							MPI_Recv(&D_loc[0], n, MPI_INT, (My_Rank + 1) % n, 1, MPI_COMM_WORLD, &status); // Send to the next process and receive data from the last process. 
							MPI_Send(&Temp[0], n, MPI_INT, My_Rank - 1, 1, MPI_COMM_WORLD);
						}
					}
				}
			}
			else
			{
				if(My_Rank == 0)
				{
					printf("The size of the arrays is different from the size of the processes ,so, the multiplication of the matrices is not possible.\n\n"); // Message if the calculation is not possible.
				}
			}

			free(Temp); // Free the B array.
			MPI_Gather(Prod_C_D_loc, n, MPI_INT, Prod_C_D, n, MPI_INT, 0, MPI_COMM_WORLD); // Gather the data to the array.
			// IV. End
			
			// Display the data and free the arrays for the Process 0.
			if(My_Rank == 0)
			{
				for(i = 0; i < n * n; i ++)
				{
					printf("Sum_C_D of C and D is: %d\n", Sum_C_D[i]);
				}
				printf("\n");
				// MPI_Reduce(Prod_B_C, Max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

				Maximum = Max[0];
				for(i = 1; i < Load; i++) // Find the maximum number.
				{
					if(Max[i] > Maximum) 
					{
						Maximum = Max[i];
					}
				}

				for(i = 0; i < n; i++) // Find the position.
				{
					if(Prod_B_C[i] == Maximum)
					{
						Position = i; 
					}
				}

				for(j = 0; j < n; j ++)
				{
					printf("Product of C and B is: %d\n", Prod_B_C[j]); // Print the result and maximum.
				}
				printf("The MAX element is: %d \nIt is in position: %d\n", Prod_B_C[Position], Position + 1);

				printf("\nProduct of A matrix and B matrix is: %d\n\n", Sum); // Print the result.
				if(n == Process)
				{
					for(i = 0; i < n * n; i ++)
					{
						printf("Product of C and D is: %d\n", Prod_C_D[i]);	// Print the result.
					}
				}
				// Free the matrices that exist only in the process 0.
				free(A);
				free(B);
				free(C);
				free(D);
				free(Prod_B_C);	
				free(Prod_A_B);
				free(Prod_C_D);
				free(Sum_C_D);
				free(Max);
			}

			// Free the matrices from all the processes.
			free(A_loc);	
			free(B_loc);
			free(C_loc);
			free(D_loc);
			free(Prod_B_C_loc);
			free(Prod_A_B_loc);
			free(Prod_C_D_loc);
			free(Sum_C_D_loc);
		}
		else
		{
			// Check if the n % p is not 0.
			int *sendcount, *displs; // Initialize the matrices.
			sendcount = (int*) malloc(Process * sizeof(int)); 
			displs = (int*) malloc(Process * sizeof(int));
			displs[0] = 0;
			
			for(i = 0; i < Process; i ++)
			{
				sendcount[i] = Load * n; // The count that each process is going to get.
			}
			if(0 < Leftover)
			{ 
				sendcount[0] = (Load + 1) * n; // How much is left for each process to get.
			}
			for(i = 1; i < Process; i ++)
			{
				if(i < Leftover) 
				{
					sendcount[i] = (Load + 1) * n; 
				}
				displs[i] = displs[i - 1] + sendcount[i - 1]; // Find from which point each process is going to get the data.
			}

			if(My_Rank < Leftover)
			{ 
				Load = Load + 1;
			}
			Prod_B_C_loc = (int *)calloc(sizeof(int), Load); 
			C_loc = (int *)malloc(n * Load * sizeof(int));
			MPI_Bcast(B, n, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the data of B.
			MPI_Scatterv(C, sendcount, displs, MPI_INT, C_loc, n * Load, MPI_INT, 0, MPI_COMM_WORLD); // Scatter the matrix.

		
			for(i = 0; i < Load; i ++)
			{
				for(j = 0; j < n; j ++)
				{
					Prod_B_C_loc[i] += B[j] * C_loc[i * n + j]; // Find the product of the row C and the array B.
				}	
			}
			for(i = 0; i < Process; i ++)
			{
				sendcount[i] = sendcount[i] / n; // The actual size we are going to get back from.
				displs[i] = displs[i] / n; // The actual posotion we are going to get back from.
			}
			MPI_Reduce(Prod_B_C_loc, Max, Load, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); // Do the calculations.
			MPI_Gatherv(Prod_B_C_loc, Load, MPI_INT, Prod_B_C, sendcount, displs, MPI_INT, 0, MPI_COMM_WORLD); // Gather the matrices into the process 0.
			
			if(My_Rank == 0)
			{
				printf("The addition of matrices C and D is not possible.\n"); // Print tha the calculation is not possible.
				printf("The multiplication of matrices A and B is not possible.\n");
				printf("The multiplication of matrices C and D is not possible.\n\n");
				
				Maximum = Max[0];
				for(i = 1; i < Load; i ++) // Find maximum.
				{
                    printf("Max[%d] = %d\n", i, Max[i] );
					if(Max[i] > Maximum)
					{
						Maximum = Max[i];
					}
				}

				for(i = 0; i < n; i ++) // Find the position.
				{
					if(Prod_B_C[i] == Maximum)
					{
						Position = i;
					}
				}

				for(j = 0; j < n; j ++) // Print the result.
				{
					printf("Product of C and B is: %d\n", Prod_B_C[j]);
				}
				printf("The MAX element is: %d \nIt is in position: %d\n", Prod_B_C[Position], Position + 1);
				free(C); // Free the matrices from the process 0.
				free(Prod_B_C);
				free(Max);
			}
			// Free the matrices from all the processes.
			free(B);
			free(sendcount);
			free(displs);
			free(C_loc);
			free(Prod_B_C_loc);
		}
		if(My_Rank == 0) // Ask user if he wants to repeat the process with new numbers and the same number of processes.
		{
			printf("1.Continue\n2.Exit\n");
			scanf("%d", &Answer);
			while(Answer != 1 && Answer != 2)
			{
				printf("Wrong, answer please press '1' to Continue or '2' to Exit\n");
				scanf("%d", &Answer);
			}
		}
		MPI_Bcast(&Answer, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the data of the answer.		
	}
	MPI_Finalize();
	return 0;
}