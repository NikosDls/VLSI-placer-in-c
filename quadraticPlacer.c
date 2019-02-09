#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "parser.h"
#include "cg\cg.h"

#include "cg\cg.c"

float solveQP(nodes *nodes, nets nets){
	int i, j, w;		// counters for the loops
	double *A;			// A array. to solve Ax = Bx, same for the By
	double *Bx, *By;	// Bx and By arrays
	double *x, *y;		// the optimal x and y coordinates
	int n;
	
	// initialize the size of arrays
	n = nodes->numberOfNodes - nodes->numberOfTerminals;
	
	// create the array a. map 2D array into 1D
	A = malloc(n * n * sizeof(double *)); 
	
	// create the vectors Bx and By
    Bx = malloc(n * sizeof(double));
    By = malloc(n * sizeof(double));
	
	// the two dimensional logical array A (nxn)
    // can be thought of as a vector of n*n entries, starting with
    // the n entries in the column 1, then the n entries in column 2
    //the entry A(I,J) is then stored in vector location I + J*N
    // and so on
    for(i = 0; i < n*n; i++){
		// initialize all elements of array A and C to zero
		A[i] = 0;
    }
    
    // start counting the execution clocks of tetris algorithm
	clock_t start = clock();
    
    for(i = 0; i < nets.numberOfNets; i++){
		for(j = 0; j < nets.array[i].netDegree - 1; j++){
			//printf("%s - > ", nodes.array[nets.array[i].netNodes[j]].name);
			
			for(w = j + 1; w < nets.array[i].netDegree; w++){
				//printf("%s  ", nodes.array[nets.array[i].netNodes[w]].name);
				// when sparse array is empty, we add the first connection and continue to the next loop
				
				if(nodes->array[nets.array[i].netNodes[j]].terminal == 1){	// if node j its terminal
					if(nodes->array[nets.array[i].netNodes[w]].terminal == 1){	// if nodes j and w are terminals
						// we just move to the next one
						continue;
					}

					// else if node j its terminal and w isnt terminal
					// increase the diagonal (w,w) element weight
					A[(nets.array[i].netNodes[w] * n) + nets.array[i].netNodes[w]]++;

					// set to the Bx and By the x and y coordinates value of the pad
					// to do later for the node i: Bx[i] = weight i * xi, correspondingly for the By
					Bx[nets.array[i].netNodes[w]] = nodes->array[nets.array[i].netNodes[j]].x;
					By[nets.array[i].netNodes[w]] = nodes->array[nets.array[i].netNodes[j]].y;
					
				}else{	// if node j isnt terminal
					if(nodes->array[nets.array[i].netNodes[w]].terminal == 1){	// if node j isnt terminal and node w its terminal
						// increase the diagonal (j,j) element weight
						A[(nets.array[i].netNodes[j] * n) + nets.array[i].netNodes[j]]++;

						// set to the Bx and By the x and y coordinates value of the pad
						// to do later for the node i: Bx[i] = weight i * xi, correspondingly for the By
						Bx[nets.array[i].netNodes[j]] = nodes->array[nets.array[i].netNodes[w]].x;
						By[nets.array[i].netNodes[j]] = nodes->array[nets.array[i].netNodes[w]].y;
						
					}else{	// if nodes j and w arent terminals
						// node j connect with the node w also the w with j
						// A[i][j] not in the array diagonal are A[i][j] = -C[i][j]
						A[(nets.array[i].netNodes[j] * n) + nets.array[i].netNodes[w]]--;
						A[(nets.array[i].netNodes[w] * n) + nets.array[i].netNodes[j]]--;
					}
				}
			}
			//printf("\n");
		}
		//printf("\n");
	}

	// complete the vectors Bx and By
	for(i = 0; i < n; i++){
		Bx[i] = Bx[i] * A[(i * n) + i];
		By[i] = By[i] * A[(i * n) + i];
	}
	
	// complete the array A
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			// elements A[i][j] on the array diagonal are A[i][j] = SUM(all elements on C[i] row) + (weight of any pad wire) 
			if(i == j){
				for(w = 0; w < n; w++){
					if(w == j){	// diagonal element on array A, has already the weight of pad wire, so we dont have to add anything else
						continue;
					}
					
					// else we add to A all other C[i] row elements 
					A[(j * n) + i] += -A[(w * n) + i];
				}
			}	
		}
	}
	/*
	printf("\nA:\n");
	
    // print A column by column
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("%lf ", A[(j * n) + i]);
		}
		printf("\n");
    }
    
    printf("\nBx:\n");
    
    // print Bx column by column
	for(i = 0; i < n; i++){
		printf("%lf ", Bx[i]);
    }
    
    printf("\n\nBy:\n");
    
    // print By column by column
	for(i = 0; i < n; i++){
		printf("%lf ", By[i]);
    }
    */
printf("solve");
	// create the results vectors x and y
	x = malloc(n * sizeof(double));
	y = malloc(n * sizeof(double));
	
	// random initial solutions
	for(i = 0; i < n; i++){
		x[i] = y[i] = 1;
	}

	// compute optimal x
	r8ge_cg(n, A, Bx, x);
	
	// compute optimal y
	r8ge_cg(n, A, By, y);
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
    //printf("\nx:\n");
    
    // save x solutions
	for(i = 0; i < n; i++){
		nodes->array[i].x = x[i];
		//printf("%lf\n", x[i]);
    }
    
    //printf("\n\ny:\n");
    
    // save y solutions
	for(i = 0; i < n; i++){
		nodes->array[i].x = x[i];
		//printf("%lf\n", y[i]);
    }
	
	// return the execution time in seconds of quadratic gp algorithm
	return seconds;	// successful return of solveQP 
}
