#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "parser.h"

void createAandB(nodes nodes, nets nets, float **Ax, float **Ay, float *Bx, float *By){
	int i, j, w;	// counter for the loops
	int n;			// A array size which is the same with array B
	int **C;
	
	// initialize the size of arrays
	n = nodes.numberOfNodes - nodes.numberOfTerminals;
	
	// create the number of columns for array C
	C = malloc(n * sizeof(int *)); 
	
    for(i = 0; i < n; i++){
		// create a new row for array C
		C[i] = malloc(n * sizeof(int));
		
		for(j = 0; j < n; j++){
			C[i][j] = 0;
		}
    }

    // create C array based on each net
	for(i = 0; i < nets.numberOfNets; i++){
		for(j = 0; j < nets.array[i].netDegree; j++){
			// check if node its terminal or not
			// if isnt terminal and node i connect to node j, we add one to C[i][j]
			if(nodes.array[nets.array[i].netNodes[j]].terminal == 0){	// node isnt terminal
				//printf("NODE\n");
				// create the C array based on each net
				if(nets.array[i].netNodesAttributes[j] == 'O'){
					for(w = 0; w < nets.array[i].netDegree; w++){
						if((w == j) || (nets.array[i].netNodesAttributes[w] == 'O')){	// node doesnt connect with herself or with nodes which have the same attribute, so we continue to the next one
							continue;	
						}
						
						// check if node w that connects with node j its terminal
						if(nodes.array[nets.array[i].netNodes[w]].terminal == 0){	// node isnt terminal
							// node j connect with the node w also the w with j
							C[nets.array[i].netNodes[j]][nets.array[i].netNodes[w]]++;
						}else if(nodes.array[nets.array[i].netNodes[w]].terminal == 1){	// node its terminal
							// add on the diagonal element of array A
							Ax[nets.array[i].netNodes[j]][nets.array[i].netNodes[j]]++;
							Ay[nets.array[i].netNodes[j]][nets.array[i].netNodes[j]]++;
							
							// set to the Bx and By the x and y coordinates value of the pad
							// to do later for the node i: Bx[i] = weight i * xi, correspondingly for the By
							Bx[nets.array[i].netNodes[j]] = nodes.array[nets.array[i].netNodes[w]].x;
							By[nets.array[i].netNodes[j]] = nodes.array[nets.array[i].netNodes[w]].y;						
						}
						
					}
				}else if(nets.array[i].netNodesAttributes[j] == 'B'){
					for(w = 0; w < nets.array[i].netDegree; w++){
						if(w == j){	// node doesnt connect with herself, so we continue to the next one
							continue;	
						}
						//printf("%d -> %d\n", j, w);
						//printf("%d -> %d\n", nets.array[i].netNodes[j], nets.array[i].netNodes[w]);
						//getch();
						
						// check if node w that connects with node j its terminal
						if(nodes.array[nets.array[i].netNodes[w]].terminal == 0){	// node isnt terminal
							// node j connect with the node w also the w with j
							C[nets.array[i].netNodes[j]][nets.array[i].netNodes[w]]++;
						}else if(nodes.array[nets.array[i].netNodes[w]].terminal == 1){	// node its terminal
							// add on the diagonal element of array A
							Ax[nets.array[i].netNodes[j]][nets.array[i].netNodes[j]]++;
							Ay[nets.array[i].netNodes[j]][nets.array[i].netNodes[j]]++;
							
							// set to the Bx and By the x and y coordinates value of the pad
							// to do later for the node i: Bx[i] = weight i * xi, correspondingly for the By
							Bx[nets.array[i].netNodes[j]] = nodes.array[nets.array[i].netNodes[w]].x;
							By[nets.array[i].netNodes[j]] = nodes.array[nets.array[i].netNodes[w]].y;								
						}
					}
				}else if(nets.array[i].netNodesAttributes[j] == 'I'){
					for(w = 0; w < nets.array[i].netDegree; w++){
						if((w == j) || (nets.array[i].netNodesAttributes[w] == 'I')){	// node doesnt connect with herself or with nodes which have the same attribute, so we continue to the next one
							continue;	
						}
						
						// check if node w that connects with node j its terminal
						if(nodes.array[nets.array[i].netNodes[w]].terminal == 0){	// node isnt terminal
							// node j connect with the node w also the w with j
							C[nets.array[i].netNodes[j]][nets.array[i].netNodes[w]]++;
						}else if(nodes.array[nets.array[i].netNodes[w]].terminal == 1){	// node its terminal
							// add on the diagonal element of array A
							Ax[nets.array[i].netNodes[j]][nets.array[i].netNodes[j]]++;	
							Ay[nets.array[i].netNodes[j]][nets.array[i].netNodes[j]]++;
							
							// set to the Bx and By the x and y coordinates value of the pad
							// to do later for the node i: Bx[i] = weight i * xi, correspondingly for the By
							Bx[nets.array[i].netNodes[j]] = nodes.array[nets.array[i].netNodes[w]].x;
							By[nets.array[i].netNodes[j]] = nodes.array[nets.array[i].netNodes[w]].y;								
						}
					}	
				}
			}
		}
	}
/*	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("%d\t", C[i][j]);	
		}
		printf("\n");
	}
	printf("\n\n");
*/	
	// complete the arrays Bx and By
	for(i = 0; i < n; i++){
		Bx[i] = Bx[i] * Ax[i][i];
		By[i] = By[i] * Ay[i][i];
	}
	
	// complete the arrays Ax and Ay
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			// elements A[i][j] on the array diagonal are A[i][j] = SUM(all elements on C[i] row) + (weight of any pad wire) 
			if(i == j){
				for(w = 0; w < n; w++){
					if(w == j){	// diagonal element on array A, has already the weight of pad wire, so we dont have to add anything else
						continue;
					}
					
					// else we add to A all other C[i] row elements 
					Ax[i][j] += C[i][w];
					Ay[i][j] += C[i][w];
				}
			}else{	// else all elements A[i][j] not in the array diagonal are A[i][j] = -C[i][j]
				Ax[i][j] = -C[i][j];
				Ay[i][j] = -C[i][j];
			}	
		}
	}
/*	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("%f\t", Ax[i][j]);	
		}
		printf("\t%f\n", Bx[i]);
	}

	printf("\n");
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("%f\t", Ay[i][j]);	
		}
		printf("\t%f\n", By[i]);
	}
*/
	// create the full matrix A+B to do gauss jordan elimination method
	for(i = 0; i < n; i++){
		Ax[i][n] = Bx[i];
		Ay[i][n] = By[i];
	}
	
	return;
}

void solveQP(nodes nodes, nets nets){
	int i, j, k;		// counters for the loops
	float **Ax, **Ay;	// Ax and Ay arrays. Ax its the A array concatenated with Bx (to solve Ax = Bx). Same for the y
	float *Bx, *By;		// Bx and By arrays
	float *x, *y;		// the optimal x and y coordinates
	float c, sum = 0.0;	// temporary variables to do the calculations
	int n;
	
	// initialize the size of arrays
	n = nodes.numberOfNodes - nodes.numberOfTerminals;
	
	// create the number of columns for array A
	Ax = malloc(n * sizeof(float *));
	Ay = malloc(n * sizeof(float *));
	
	// create the array Bx and By
    Bx = malloc(n * sizeof(float));
    By = malloc(n * sizeof(float));
	
    for(i = 0; i < n; i++){
    	// initialize all elements of array B to zero
    	Bx[i] = 0;
    	By[i] = 0;
    	
		// create a new row for array A
    	// the 1 more element in each row its for the Bx or By value
		Ax[i] = malloc((n + 1) * sizeof(float));
		Ay[i] = malloc((n + 1) * sizeof(float));
		
		for(j = 0; j < n; j++){
			// initialize all elements of array A and C to zero
			Ax[i][j] = 0;
			Ay[i][j] = 0;
		}
    }
    
	// create the array x and y
	x = malloc(n * sizeof(float));
	y = malloc(n * sizeof(float));
	
	// create the A and B arrays
	createAandB(nodes, nets, Ax, Ay, Bx, By);

	// first we calculate all x coordinates
	// then we do the same process for the y coordinates
	
	// loop for the generation of upper triangular matrix
	for(j = 0; j < n; j++){
		for(i = 0; i < n; i++){
			if(i > j){
				c = Ax[i][j] / Ax[j][j];
				for(k = 0; k < (n + 1); k++){
					Ax[i][k] = Ax[i][k] - c * Ax[j][k];
				}
			}
		}
	}
    
	// find the Xn
	x[n-1] = Ax[n-1][n] / Ax[n-1][n-1];

	// this loop is for backward substitution
	for(i = n - 2; i >= 0; i--){
		sum = 0.0;
		for(j = i; j < n; j++){
			sum += Ax[i][j] * x[j];
		}
		// find the Xn-1 ... X1
		x[i] = (Ax[i][n] - sum) / Ax[i][i];
	}
	
	printf("\nX coordinates for each node:");
	for(i = 0; i < n; i++){
		printf("\nX%d = %.4f", i+1, x[i]);
	}
	
	// loop for the generation of upper triangular matrix
	for(j = 0; j < n; j++){
		for(i = 0; i < n; i++){
			if(i > j){
				c = Ay[i][j] / Ay[j][j];
				for(k = 0; k < (n + 1); k++){
					Ay[i][k] = Ay[i][k] - c * Ay[j][k];
				}
			}
		}
	}
    
	// find the Yn
	y[n-1] = Ay[n-1][n] / Ay[n-1][n-1];

	// this loop is for backward substitution
	for(i = n - 2; i >= 0; i--){
		sum = 0.0;
		for(j = i; j < n; j++){
			sum += Ay[i][j] * y[j];
		}
		// find the Xn-1 ... X1
		y[i] = (Ay[i][n] - sum) / Ay[i][i];
	}

	printf("\n\nY coordinates for each node:");
	for(i = 0; i < n; i++){
		printf("\nY%d = %.4f", i+1, y[i]);
	}
	
	return;
}
