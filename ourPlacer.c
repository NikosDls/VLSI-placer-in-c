#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "parser.h"
#include "quadraticPlacer.h"
#include "ourPlacer.h"

#define CG_DESCENT_IMPLEMENTATION
#include "cg_nl/cg_descent.h"

double averageWidth;

void createH0(hypergraph *H, int numberOfNodes){
	// create the first level of hypergraph
	H->array = malloc(sizeof(level));
		
	// initially the first level of graph have all nodes
	H->array[0].blockNumber = numberOfNodes;
	
	// set the range of the first level
	H->array[0].range[0] = 0;
	H->array[0].range[1] = numberOfNodes - 1;

	// set the number of levels to one
	H->numberOfLevels = 1;
	
	return;	// successful return of createH0
}

void firstChoiceClustering(hypergraph *H, int level){
	int nodes; // number of nodes for the next level
	
	// set he number of nodes for the next level
	nodes = H->array[level - 1].blockNumber / 5;

	// increase the number of hypergraph levels
	H->array = realloc(H->array, (level + 1) * sizeof(struct level));
	
	// set the range for the new level 
	H->array[level].range[0] = H->array[level - 1].range[0] + nodes + 1;
	H->array[level].range[1] = H->array[level - 1].range[1];

	// set the number of nodes for the new level
	H->array[level].blockNumber = H->array[level].range[1] - H->array[level].range[0];
	
	// update the number of nodes for the previous level
	H->array[level - 1].blockNumber = nodes;

	// update the range for the previous level
	H->array[level - 1].range[1] = H->array[level - 1].range[0] + nodes;

	// save the number of levels
	H->numberOfLevels = level + 1;
	
	//printf("Hlevel   %d\tblock number %d\trange %d-%d\nHlevel-1 %d\tblock number %d\trange %d-%d\n\n", level, H->array[level].blockNumber, H->array[level].range[0],H->array[level].range[1], level-1, H->array[level-1].blockNumber, H->array[level-1].range[0],H->array[level-1].range[1]);
	//getch();
	
	return;	// successful return of firstChoiceClustering
}

void swapC(int *a, int *b){
	int temp = *a; 
	*a = *b; 
	*b = temp; 
}

// the presence of equal elements requires special care in quicksort
// so we use Hoare's partition method
int partitionC(nodes nodes, int arr[], int low, int high){ 
	double pivot = nodes.array[arr[high]].connectivity;	// pivot
	int i = low - 1, j = high + 1;
	
	while(1){ 
		// find leftmost element greater than 
		// or equal to pivot 
		do{ 
			i++; 
		}while(nodes.array[arr[i]].connectivity < pivot); 
   
		// find rightmost element smaller than 
		// or equal to pivot 
		do{ 
			j--; 
        }while(nodes.array[arr[j]].connectivity > pivot); 

		// if two pointers met 
		if(i >= j){
			return j; 
		}
		
		// swap i and j element
		swapC(&arr[i], &arr[j]); 
    }
}

/*
// Lomuto's method is simple and easier to implement, but should not be used for implementing a library sorting method
int partitionC(nodes nodes, int arr[], int low, int high){
	double pivot = nodes.array[arr[high]].connectivity;	// pivot 
	int i = (low - 1);	// index of smaller element 
	int j;
	
	for (j = low; j <= (high - 1); j++){ 
		// if current element is smaller than or equal to pivot
		if(nodes.array[arr[j]].connectivity <= pivot){ 
			i++;    // increment index of smaller element 
			
			// swap i and j element
			swapC(&arr[i], &arr[j]); 
		}
	}

	// swap i + 1 and high element
	swapC(&arr[i + 1], &arr[high]); 
	
	return (i + 1);
}
*/

void quickSortC(nodes nodes, int arr[], int low, int high){ 
	if(low < high){ 
		// pi is partitioning index, arr[p] is now 
		// at right place
		int pi = partitionC(nodes, arr, low, high); 
		
		// separately sort elements before 
		// partition and after partition
		quickSortC(nodes, arr, low, pi - 1); 
		quickSortC(nodes, arr, pi + 1, high); 
	}
}

connectivitySortedNodes sortNodesByConnectivity(nodes nodes){	// sort based on connectivity
	int i, j;	// counters for the loops
	int temp;	// temporary variable for the swaps
	connectivitySortedNodes sortedNodes;	// sorted nodes
	
	// set number of sorted nodes
	sortedNodes.numberOfNodes = nodes.numberOfNodes - nodes.numberOfTerminals;
	
	// creating the sorted array
	sortedNodes.array = malloc(sortedNodes.numberOfNodes * sizeof(int));
	
	// initialze the sorted nodes ids
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		// nodes isnt sorted yet
		// so the ids are in ascending order (0...n), as they are in the file
		sortedNodes.array[i] = i;
	}
	
	/*
	// by using Hoare's partitioning method, we dont have to care for equal elements 
	// variables for the small random number
	double randomNum;
	srand((unsigned int)time(NULL));
	
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		// add a very small random number to each node connectivity
		// so the quicksort do not fail while right-most pivoting
		randomNum = ((double) rand()) / RAND_MAX;
        randomNum /= 70000.0;
        
        nodes.array[i].connectivity += randomNum;
        
    	//printf("%d %s %lf\n", i, nodes.array[sortedNodes.array[i]].name, nodes.array[sortedNodes.array[i]].connectivity);
	}
	*/
	
	// sort the nodes based on their connectivity
	quickSortC(nodes, sortedNodes.array, 0, sortedNodes.numberOfNodes - 1);

	/*
	// sort the nodes for the clustering process
	for(i = sortedNodes.numberOfNodes - 1; i >= 0; i--){
		for(j = 0; j < i; j++){
			if(nodes.array[sortedNodes.array[i]].connectivity < nodes.array[sortedNodes.array[j]].connectivity){	// values must swap
				// swap the values
            	temp = sortedNodes.array[i];
            	sortedNodes.array[i] = sortedNodes.array[j];
            	sortedNodes.array[j] = temp;
			}
		}
	}	
	
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		printf("%d %s %lf\n", i, nodes.array[sortedNodes.array[i]].name, nodes.array[sortedNodes.array[i]].connectivity);
	}
	*/
	
	// return the sorted nodes
	return sortedNodes;	// successful return of sortNodesByConnectivity
}

void decluster(hypergraph *H, int level){
	// increase the next level nodes
	H->array[level - 1].blockNumber += H->array[level].blockNumber;
	
	// change the ranges
	H->array[level - 1].range[1] = H->array[level].range[1];
	
	return;	// successful return of decluster
}

void calculateMbPb(int horizontally, int vertically, binGrids *B, chip chip){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;		// counters to calculate Pb
	
	// initialize the center, Db, Mb and Pb for each bin
	for(j = 0; j < vertically; j++){
		for(w = 0; w < horizontally; w++){
			if(w == (horizontally - 1)){
				B->array[j][w].xCenter = (double) (B->xRange[w] + (B->xRange[w] + B->wb)) / 2;
			}else{
				B->array[j][w].xCenter = (double) (B->xRange[w + 1] + B->xRange[w]) / 2;
			}
			
			if(j == (vertically - 1)){		
				B->array[j][w].yCenter = (double) (B->yRange[j] + (B->yRange[j] + B->hb)) / 2;
			}else{
				B->array[j][w].yCenter = (double) (B->yRange[j + 1] + B->yRange[j]) / 2;
			}
			//printf("%lf %lf\t", B->array[j][w].xCenter, B->array[j][w].yCenter);
			
			B->array[j][w].Db = 0;
			B->array[j][w].Mb = 0;
			B->array[j][w].Pb = 0;
		}
		//printf("\n\n");
	}

	// calculate Pb
	// if circuit have preplaced nodes in the chip, we have to calculate the Pb for each bin
	if(chip.pbInChip == 1){
		// initialize y counter
		yCounter = 0;
			
		// start scanning for the preplaced areas
		for(j = 0; j < chip.numberOfRows; j++){
			for(w = 0; w < chip.array[j].height; w++){
				// initialize x counter
				xCounter = 0;
				
				for(k = 0; k < chip.array[j].width; k++){
					// check if slot in the chip is unavailable and increase the Pb in the corresponding bin
					if(chip.array[j].mixedArray[w][k] == notAvailable){
						// increase the preplaced area in the bin (j,w)
						B->array[yCounter][xCounter].Pb++;	
					}
					
					// move to the next bin in x direction
					if(k == round(B->xRange[xCounter])){
						xCounter++;	
					}
				}
			}
				
			// move to the next bin in y direction	
			if(w == round(B->yRange[yCounter])){
				yCounter++;	
			}
		}
	}
	
	// calculate Mb for each bin
	for(j = 0; j < vertically; j++){
		for(w = 0; w < horizontally; w++){
			B->array[j][w].Mb = tdensity * ((B->wb * B->hb) - B->array[j][w].Pb);
			//printf("%lf \t", B->array[j][w].Mb);
		}
		//printf("\n");
	}

	return;	// successful return of calculateMbPb	
}

void calculateDb(int horizontally, int vertically, binGrids *B, hypergraph H, int i, nodes nodes, connectivitySortedNodes sortedNodes){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;			// counters to calculate Pb
	double dx, dy;					// center (of node) to center (of bin) distance, in x and y directions
	double px, py, ax, bx, ay, by;	// variables to calculate Db
	
	for(k = 0; k < H.array[i].blockNumber; k++){
		// calculate ax
		ax = (double) 4 / ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 2 * B->wb) * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 4 * B->wb));
	
		// calculate bx
		bx = (double) 2 / (B->wb * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 4 * B->wb));
			
		// calculate ay
		ay = (double) 4 / ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 2 * B->wb) * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 4 * B->wb));
				
		// calculate by
		by = (double) 2 / (B->wb * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 4 * B->wb));
							
		// calculate Db for each bin
		for(j = 0; j < vertically; j++){
			for(w = 0; w < horizontally; w++){
				// calculate center to center distance in x direction
				dx = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xCenter - B->array[j][w].xCenter);
				
				if((dx >= 0) && 
				   (dx <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb))){
					
					// calculate px
					px = 1 - ax * pow(dx, 2);				
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb) <= dx) &&
						   (dx <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb)))){
					
					// calculate px
					px = bx * pow(dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb), 2);
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb))) <= dx){
					// calculate px
					px = 0;
				}
				
				// calculate center to center distance in y direction
				dy = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yCenter - B->array[j][w].yCenter);
				
				if((dy >= 0) && 
				   (dy <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb))){
	
					// calculate py
					py = 1 - ay * pow(dy, 2);				
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb) <= dy) &&
						   (dy <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb)))){

					// calculate py
					py = by * pow(dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb), 2);
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb))) <= dy){
					// calculate py
					py = 0;
				}
				
				// calculate Db for the bin (j,w) and node k
				B->array[j][w].Db += nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength
									 * nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength
									 * px * py; 
				
				/*
				// px or py cant be greater than 1
				if(px > 1 || py > 1){
					printf("px %lf\npy %lf\n\n", px, py);
					getch();
				}
				*/		
			}
		}
	}

	return;	// successful return of calculateDb
}

void DbGradientX(int horizontally, int vertically, binGrids *B, hypergraph H, int i, nodes nodes, connectivitySortedNodes sortedNodes, double *gradient){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;			// counters to calculate Pb
	double dx, dy;					// center (of node) to center (of bin) distance, in x and y directions
	double px, py, ax, bx, ay, by;	// variables to calculate Db
	
	// initialize gradient vector to zero
	for(j = 0; j < H.array[i].blockNumber; j++){
		gradient[j] = 0;
	}
	
	for(k = 0; k < H.array[i].blockNumber; k++){
		// calculate ax
		ax = (double) 4 / ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 2 * B->wb) * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 4 * B->wb));
	
		// calculate bx
		bx = (double) 2 / (B->wb * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 4 * B->wb));
			
		// calculate ay
		ay = (double) 4 / ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 2 * B->wb) * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 4 * B->wb));
				
		// calculate by
		by = (double) 2 / (B->wb * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 4 * B->wb));
					
		// calculate Db gradient for each bin
		for(j = 0; j < vertically; j++){
			for(w = 0; w < horizontally; w++){
				// calculate center to center distance in x direction
				dx = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xCenter - B->array[j][w].xCenter);
				
				if((dx >= 0) && 
				   (dx <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb))){
					
					// calculate px
					px = -2 * ax * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xCenter - B->array[j][w].xCenter);				
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb) <= dx) &&
						   (dx <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb)))){

					// calculate px
					if(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xCenter > B->array[j][w].xCenter){
						px = 2 * bx * (dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb));
					}else{
						px = -2 * bx * (dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb));
					}	
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb))) <= dx){
					// calculate px
					px = 0;
				}
				
				// calculate center to center distance in y direction
				dy = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yCenter - B->array[j][w].yCenter);
				
				if((dy >= 0) && 
				   (dy <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb))){
					
					// calculate py
					py = 1 - ay * pow(dy, 2);				
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb) <= dy) &&
						   (dy <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb)))){					
					
					// calculate py
					py = by * pow(dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb), 2);		
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb))) <= dy){
					// calculate py
					py = 0;
				}
				
				// calculate Db gradient for the bin (j,w) and node k
				gradient[k] += nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength
								* nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength
								* px * py;		
			}
		}
	}

	return;	// successful return of DbGradientX
}

void DbGradientY(int horizontally, int vertically, binGrids *B, hypergraph H, int i, nodes nodes, connectivitySortedNodes sortedNodes, double *gradient){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;			// counters to calculate Pb
	double dx, dy;					// center (of node) to center (of bin) distance, in x and y directions
	double px, py, ax, bx, ay, by;	// variables to calculate Db
	
	// initialize gradient vector to zero
	for(j = 0; j < H.array[i].blockNumber; j++){
		gradient[j] = 0;
	}
		
	for(k = 0; k < H.array[i].blockNumber; k++){
		// calculate ax
		ax = (double) 4 / ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 2 * B->wb) * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 4 * B->wb));
	
		// calculate bx
		bx = (double) 2 / (B->wb * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength + 4 * B->wb));
			
		// calculate ay
		ay = (double) 4 / ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 2 * B->wb) * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 4 * B->wb));
				
		// calculate by
		by = (double) 2 / (B->wb * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength + 4 * B->wb));
					
		// calculate Db gradient for each bin
		for(j = 0; j < vertically; j++){
			for(w = 0; w < horizontally; w++){
				// calculate center to center distance in x direction
				dx = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xCenter - B->array[j][w].xCenter);
				
				if((dx >= 0) && 
				   (dx <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb))){
					
					// calculate px
					px = 1 - ax * pow(dx, 2);				
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb) <= dx) &&
						   (dx <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb)))){
					
					// calculate px
					px = bx * pow(dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb), 2);		
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb))) <= dx){
					// calculate px
					px = 0;
				}
				
				// calculate center to center distance in y direction
				dy = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yCenter - B->array[j][w].yCenter);
				
				if((dy >= 0) && 
				   (dy <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb))){
					
					// calculate py
					py = -2 * ay * (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yCenter - B->array[j][w].yCenter);				
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb) <= dy) &&
						   (dy <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb)))){					
					
					// calculate py
					if(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yCenter > B->array[j][w].yCenter){
						py = 2 * by * (dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb));
					}else{
						py = -2 * by * (dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb));
					}		
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb))) <= dy){
					// calculate py
					py = 0;
				}
				
				// calculate Db gradient for the bin (j,w) and node k
				gradient[k] += nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength
								* nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength
								* px * py;		
			}
		}
	}

	return;	// successful return of DbGradientY
}

double W(nodes nodes, nets nets, chip chip){
	int i, j;	// counters for the loops
	double g;	// g parameter for the log-sum-exp model
	double totalWireLength = 0.0;	// total wirelength
	double logxk, log_xk, logyk, log_yk;	// temporary variables for the wirelength calculation
	
	// g have the value of the 1% of the chip width
	g = 0.01 * averageWidth;

	// calculate wire length for each net
	// with log-sum-exp model
	for(i = 0; i < nets.numberOfNets; i++){
		// initialize the sum of four paramaters for each net
		logxk = log_xk = logyk = log_yk = 0.0;

		for(j = 0; j < nets.array[i].netDegree; j++){
			// calculate the sum of four paramater for each net
			logxk  += exp((double) nodes.array[nets.array[i].netNodes[j]].x / g);
			log_xk += exp((double) -nodes.array[nets.array[i].netNodes[j]].x / g);
			logyk  += exp((double) nodes.array[nets.array[i].netNodes[j]].y / g);
			log_yk += exp((double) -nodes.array[nets.array[i].netNodes[j]].y / g);
			/*
			printf("e^(%.2f / %.2lf) = %.20lf\n"
				   "e^(%.2f / %.2lf) = %.20lf\n"
				   "e^(%.2f / %.2lf) = %.20lf\n"
				   "e^(%.2f / %.2lf) = %.20lf\n", nodes.array[nets.array[i].netNodes[j]].x,  g, logxk,
												  -nodes.array[nets.array[i].netNodes[j]].x, g, log_xk, 
											      nodes.array[nets.array[i].netNodes[j]].y,  g, logyk, 
												  -nodes.array[nets.array[i].netNodes[j]].y, g, log_yk);
			getch();
			*/
		}

		// calculate the log value of four paramaters for net j
		logxk  = log(logxk);
		log_xk = log(log_xk);
		logyk  = log(logyk);
		log_yk = log(log_yk);
		
		// add the net i wirelength to the total wirelength
		totalWireLength += (logxk + log_xk + logyk + log_yk);
	}
	
	// multiply the total wirelegth with g 	
	totalWireLength *= g;
	
	// return the total wirelength
	return totalWireLength;	// successful return of W 
}

double PDxW(nodes nodes, nets nets, chip chip, char *nodeName){
	int i, j;	// counters for the loops
	double g;	// g parameter for the log-sum-exp model
	double logxk1, log_xk1, logxk2, log_xk2, logyk, log_yk;	// temporary variables for the wirelength calculation
	
	const double delta = 1.0e-6;	// small step to find the derivative
	double sumSubDelta = 0.0, sumAddDelta = 0.0;	// x - delta and x + delta sum
	
	// g have the value of the 1% of the chip width
	g = 0.01 * averageWidth;

	// calculate wire length for each net
	// with log-sum-exp model
	for(i = 0; i < nets.numberOfNets; i++){
		// initialize the sum of four paramaters for each net
		logxk1 = log_xk1 = logxk2 = log_xk2 = logyk = log_yk = 0.0;
		
		for(j = 0; j < nets.array[i].netDegree; j++){
			// calculate the sum of four paramater for each net
			if(strcmp(nodeName, nodes.array[nets.array[i].netNodes[j]].name) == 0){	// derivative node
				// x - delta
				logxk1  += exp((double) (nodes.array[nets.array[i].netNodes[j]].x - delta) / g);
				log_xk1 += exp((double) (-nodes.array[nets.array[i].netNodes[j]].x - delta) / g);
				
				// x + delta
				logxk2  += exp((double) (nodes.array[nets.array[i].netNodes[j]].x + delta) / g);
				log_xk2 += exp((double) (-nodes.array[nets.array[i].netNodes[j]].x + delta) / g);					
			}else{	// normal node
				// x
				logxk1  += exp((double) nodes.array[nets.array[i].netNodes[j]].x / g);
				log_xk1 += exp((double) -nodes.array[nets.array[i].netNodes[j]].x / g);
				
				logxk2  = logxk1;
				log_xk2 = log_xk1;
			}
			
			// y is fixed since we derive the node for x
			logyk  += exp((double) nodes.array[nets.array[i].netNodes[j]].y / g);
			log_yk += exp((double) -nodes.array[nets.array[i].netNodes[j]].y / g);
		}
		
		// calculate the log value of four paramaters for net j
		logxk1  = log(logxk1);
		log_xk1 = log(log_xk1);
		logxk2  = log(logxk2);
		log_xk2 = log(log_xk2);
		
		logyk  = log(logyk);
		log_yk = log(log_yk);
		
		// add the net i wirelength to the total wirelength
		sumSubDelta += (logxk1 + log_xk1 + logyk + log_yk);
		
		// add the net i wirelength to the total wirelength
		sumAddDelta += (logxk2 + log_xk2 + logyk + log_yk);
	}
	
	// multiply the total wirelegth with g 
	sumSubDelta *= g;	
	sumAddDelta *= g;

	// return the partial derivative x for the node
	return (sumAddDelta - sumSubDelta) / (2 * delta);	// successful return of PDxW 
}

void wGradientX(nodes nodes, nets nets, chip chip, double *gradient){
	int i, j;	// counters for the loops
	double g;	// g parameter for the log-sum-exp model
	double totalWireLength = 0.0;	// total wirelength
	
	double commonPos, commonNeg;	// common constant for all net
	double *expDerivativesPos, *expDerivativesNeg;	// exponential derivative for each node k in the net (e^(xk / g))
	
	// initialize gradient vector to zero
	for(i = 0; i < nodes.numberOfNodes; i++){
		gradient[i] = 0;
	}
	
	// g have the value of the 1% of the chip width
	g = 0.01 * averageWidth;

	// calculate wire length for each net
	// with log-sum-exp model
	for(i = 0; i < nets.numberOfNets; i++){
		// initialize common constant for net i
		commonPos = commonNeg = 0.0;
		
		// creating the exponential vectors
		expDerivativesPos = malloc(nets.array[i].netDegree * sizeof(double));
		expDerivativesNeg = malloc(nets.array[i].netDegree * sizeof(double));
		
		// calculate exponential derivatives and common constant for each node in the net i
		for(j = 0; j < nets.array[i].netDegree; j++){
			// find the e^(xk / g) and e^(-xk / g)
			expDerivativesPos[j] = exp((double) nodes.array[nets.array[i].netNodes[j]].x / g);		
			expDerivativesNeg[j] = exp((double) -nodes.array[nets.array[i].netNodes[j]].x / g);	
			
			// calculate the denominator of common constant
			commonPos += expDerivativesPos[j];
			commonNeg += expDerivativesNeg[j];
		}
		
		// calculate the common constant
		commonPos = (double) 1 / commonPos;
		commonNeg = (double) 1 / commonNeg;
		
		// calculate gradient for each node in the net i
		for(j = 0; j < nets.array[i].netDegree; j++){
			// add to gradient the new gradient
			/*
			printf("commonPos %20.25lf\n"
				   "1/g       %20.25lf\n"
				   "expDPos   %20.25lf\n\n", commonPos, (double) 1 / g, expDerivativesPos[j]);
				   
			printf("commonNeg %20.25lf\n"
				   "1/g       %20.25lf\n"
				   "expDNeg   %20.25lf", commonNeg, (double) 1 / g, expDerivativesNeg[j]);getch();
			*/	   
			gradient[nets.array[i].netNodes[j]] += commonPos * expDerivativesPos[j];
			gradient[nets.array[i].netNodes[j]] += commonNeg * -expDerivativesNeg[j];
			
			//printf("\n\n%d node %s == %.50lf\n\n", nets.array[i].netNodes[j], nodes.array[nets.array[i].netNodes[j]].name, gradient[nets.array[i].netNodes[j]]);getch();
		}
		
		// free the exponential derivatives
		free(expDerivativesPos);
		free(expDerivativesNeg);
	}
	
	return;	// successful return of gradient 
}

void wGradientY(nodes nodes, nets nets, chip chip, double *gradient){
	int i, j;	// counters for the loops
	double g;	// g parameter for the log-sum-exp model
	double totalWireLength = 0.0;	// total wirelength
	
	double commonPos, commonNeg;	// common constant for all net
	double *expDerivativesPos, *expDerivativesNeg;	// exponential derivative for each node k in the net (e^(xk / g))

	// initialize gradient vector to zero
	for(i = 0; i < nodes.numberOfNodes; i++){
		gradient[i] = 0;
	}
	
	// g have the value of the 1% of the chip width
	g = 0.01 * averageWidth;

	// calculate wire length for each net
	// with log-sum-exp model
	for(i = 0; i < nets.numberOfNets; i++){
		// initialize common constant for net i
		commonPos = commonNeg = 0.0;
		
		// creating the exponential vectors
		expDerivativesPos = malloc(nets.array[i].netDegree * sizeof(double));
		expDerivativesNeg = malloc(nets.array[i].netDegree * sizeof(double));
		
		// calculate exponential derivatives and common constant for each node in the net i
		for(j = 0; j < nets.array[i].netDegree; j++){
			// find the e^(yk / g) and e^(-yk / g)
			expDerivativesPos[j] = exp((double) nodes.array[nets.array[i].netNodes[j]].y / g);		
			expDerivativesNeg[j] = exp((double) -nodes.array[nets.array[i].netNodes[j]].y / g);	
			
			// calculate the denominator of common constant
			commonPos += expDerivativesPos[j];
			commonNeg += expDerivativesNeg[j];
		}
		
		// calculate the common constant
		commonPos = (double) 1 / commonPos;
		commonNeg = (double) 1 / commonNeg;
		
		// calculate gradient for each node in the net i
		for(j = 0; j < nets.array[i].netDegree; j++){
			// add to gradient the new gradient
			/*
			printf("commonPos %20.25lf\n"
				   "1/g       %20.25lf\n"
				   "expDPos   %20.25lf\n\n", commonPos, (double) 1 / g, expDerivativesPos[j]);
				   
			printf("commonNeg %20.25lf\n"
				   "1/g       %20.25lf\n"
				   "expDNeg   %20.25lf", commonNeg, (double) 1 / g, expDerivativesNeg[j]);getch();
			*/	   
			gradient[nets.array[i].netNodes[j]] += commonPos * expDerivativesPos[j];
			gradient[nets.array[i].netNodes[j]] += commonNeg * -expDerivativesNeg[j];
			
			//printf("\n\n%d node %s == %.50lf\n\n", nets.array[i].netNodes[j], nodes.array[nets.array[i].netNodes[j]].name, gradient[nets.array[i].netNodes[j]]);getch();
		}
		
		// free the exponential derivatives
		free(expDerivativesPos);
		free(expDerivativesNeg);
	}
	
	return;	// successful return of gradient 
}

float ourPlacerGP(nodes nodes, nets nets, chip chip, int nmax, int totalMovableArea){
	int i, j, w;	// counters for the loops
	int level;		// counter for the hypergaph levels	
	hypergraph H;	// hypergraph
	binGrids B;		// bins
	connectivitySortedNodes sortedNodes;	// sorted nodes based on their connectivity
	int grid_num;	// dimensions of grids
	int numberOfBinsVertically, numberOfBinsHorizontally;	// number of bins vertically and horizontally
	double xCenter, yCenter;	// center of the chip
	double l;			// l constant	 
	double *wGXt, *wGYt;// temporary partial derivative vectors of function W
	double *wGX, *wGY;	// partial derivative vectors of function W
	double *dGX, *dGY;	// partial derivative vectors of function Db
	double sN, sD;		// numerator and denominator sum
	double sum;			// sum to calculate average chip width and overflow ratio
	double overflow_ratio;	// over flow ratio
	
	// calculate average chip width
	sum = 0.0;
	for(i = 0; i < chip.numberOfRows; i++){
		sum += chip.array[i].width;
	}
	averageWidth = sum / chip.numberOfRows;
	//printf("%lf", averageWidth);
	
	// mixed size circuit
	createH0(&H, nodes.numberOfNodes - nodes.numberOfTerminals);
	
	// get the sorted nodes based on their connectivity	
	sortedNodes = sortNodesByConnectivity(nodes);	

	// initialize level to zero
	level = 0;
	
	// create the hypergraph
	while(H.array[level].blockNumber > nmax){
		level++;	// increase the number of levels
		firstChoiceClustering(&H, level);	// create the new hypergraph level, by dividing the previous level
	}
	
	for(i = 0; i < H.numberOfLevels; i++){
		printf("H%d\tNumber Of Blocks: %ld\t Blocks ID: %ld - %ld\n", i, H.array[i].blockNumber, H.array[i].range[0], H.array[i].range[1]);
	}

	// set the nodes level
	for(i = 0; i < H.numberOfLevels; i++){
		for(j = H.array[i].range[0]; j <= H.array[i].range[1]; j++){
			nodes.array[sortedNodes.array[j]].level = i;
			//printf("%d %d\n", j, sortedNodes.array[j]);
		}
		
	}
	
	/*
	for(i = 0; i < nodes.numberOfNodes - nodes.numberOfTerminals; i++){
		printf("%d %lf\n", nodes.array[i].level, nodes.array[i].connectivity);
		if((i % 200) == 0){
			getch();
		}
	}
	*/
	
	// find the center of the chip
	xCenter = (double) averageWidth / 2;
	yCenter = (double) (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height) / 2;
	//printf("center: %lf %lf\n", xCenter, yCenter);

	// solve qp
	for(i = 0; i < H.array[level].blockNumber; i++){
		// set the nodes id in H level
		nodes.array[sortedNodes.array[H.array[level].range[0] + i]].id = i;
	}

	// call QP for the H level
	// if we can solve QP, save the solutions 
	// else we keep the center coordinates as initial solution
	//QP(&nodes, nets, H, sortedNodes, 30);
	
	// for the nodes which isnt in the last level
	// initialize nodes positions as the center of the chip
	for(i = 0; i < nodes.numberOfNodes - nodes.numberOfTerminals; i++){
		if(nodes.array[i].level != H.numberOfLevels - 1){
			nodes.array[i].x = nodes.array[i].xCenter = xCenter;
			nodes.array[i].y = nodes.array[i].yCenter = yCenter;
		}
	}
	
	// start counting the execution clocks of our algorithm
	clock_t start = clock();
	
	for(i = level; i >= 0; i--){
		//printf("%d\n", i);
		
		// calculate the dimensions of bins vertically and horizontally
		grid_num = round(sqrt(H.array[i].blockNumber));
		printf("\n\nDimensions of bins horizontally and vertically: %d\n", grid_num);
		
		// calculate the number of bins
		numberOfBinsVertically = chip.array[0].width / grid_num;
		numberOfBinsHorizontally = (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height) / grid_num;
		printf("Number of bins horizontally %d and vertically %d\n", numberOfBinsVertically, numberOfBinsHorizontally);
		
		// calculate the number of bins
		B.numberOfBins = numberOfBinsVertically * numberOfBinsHorizontally;
		printf("Total number of bins %d\n", B.numberOfBins);
			
		// creating the all rows for slot array
		B.array = malloc(numberOfBinsVertically * sizeof(bin *));
		
		for(j = 0; j < numberOfBinsVertically; j++){
			// creating the all columns in each row for slot array
			B.array[j] = malloc(numberOfBinsHorizontally * sizeof(bin));
		}
		
		// calculate the width and height of each bin
		B.wb = (double) chip.array[0].width / numberOfBinsHorizontally;
		B.hb = (double) (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height) / numberOfBinsVertically;
		printf("Width of each bin (Wb) = %lf\nHeight of each bin(Hb) = %lf\n", B.wb, B.hb);

		// set the sizes of the x and y ranges
		B.xRange = malloc(numberOfBinsHorizontally * sizeof(double));
		B.yRange = malloc(numberOfBinsVertically * sizeof(double));
	
		// set the x ranges
		for(j = 0; j < numberOfBinsHorizontally; j++){
			B.xRange[j] = B.wb * j;
			//printf("%lf\t", B.xRange[j]);
		}
		//printf("\n\n");
		
		// set the y ranges
		for(j = 0; j < numberOfBinsVertically; j++){
			B.yRange[j] = B.hb * j;
			//printf("%lf\t", B.yRange[j]);
		}
		//printf("\n\n");
		
		// calculate Mb and Pb
		calculateMbPb(numberOfBinsHorizontally, numberOfBinsVertically, &B, chip);
		
		// calculate Db
		calculateDb(numberOfBinsHorizontally, numberOfBinsVertically, &B, H, i, nodes, sortedNodes);
		
		// create temporary the partial derivative vectors
		wGXt = malloc(nodes.numberOfNodes * sizeof(double));
		wGYt = malloc(nodes.numberOfNodes * sizeof(double));
		
		// calculate the partial derivatives of function W
		wGradientX(nodes, nets, chip, wGXt);
		wGradientY(nodes, nets, chip, wGYt);
		
		// create the partial derivative vectors
		wGX = malloc(H.array[i].blockNumber * sizeof(double));
		wGY = malloc(H.array[i].blockNumber * sizeof(double));
		
		// keep only derivatives of Hlevel i
		for(j = 0; j < H.array[i].blockNumber; j++){
			wGX[i] = wGXt[sortedNodes.array[H.array[i].range[0] + j]];
			wGY[i] = wGYt[sortedNodes.array[H.array[i].range[0] + j]];
		}

		// we dont need the temporary partial derivatives
		// so we free them
		free(wGXt);
		free(wGYt);
		
		// create the partial derivative vectors
		dGX = malloc(H.array[i].blockNumber * sizeof(double));
		dGY = malloc(H.array[i].blockNumber * sizeof(double));
		
		// calculate the partial derivatives of function Db
		DbGradientX(numberOfBinsHorizontally, numberOfBinsVertically, &B, H, i, nodes, sortedNodes, dGX);
		DbGradientY(numberOfBinsHorizontally, numberOfBinsVertically, &B, H, i, nodes, sortedNodes, dGY);

		// initialize sums
		sN = sD = 0.0;
		
		// calculate numerator and denominator
		// to calculate the l constant
		for(j = 0; j < H.array[i].blockNumber; j++){
			// add to the numerator sum, the absolut value of the partial x and y derivative of function W 
			sN += fabs(wGX[j] + wGY[j]);
			
			// add to the denominator sum, the absolut value of the partial x and y derivative of function Db
			sD += fabs(dGX[j] + dGY[j]);
		}
		//printf("%lf / %lf\n", sN, sD);
		
		// calculate l
		l = (double) sN / sD;
		
		//printf("l = %lf\n", l);
		//getch();
		
		// we calculated l, and we dont need the partial derivatives anymore
		// so we free them
		free(wGX);
		free(wGY);
		free(dGX);
		free(dGY);
		
		// objective function minimization
		do{
		
		// increase l by two times
		l *= 2;

		// calculate overflow ratio
		sum = 0.0;
		for(j = 0; j < numberOfBinsVertically; j++){
			for(w = 0; w < numberOfBinsHorizontally; w++){
				sum += max((B.array[j][w].Db - B.array[j][w].Mb), 0);
			}
		}
		
		overflow_ratio = (double) sum / totalMovableArea;
		//printf("OR = %lf\n", overflow_ratio);	
		//getch();	
		break;
		
		// until overflow ratio its small enough
		}while(overflow_ratio > 0.05);
			
		//printf("\n%d - %d", H.array[i].range[0], H.array[i].range[1]);
		if(i != 0){	// if we are not in the last spread (level 0)
			// decluster nodes from current level (i) to the next level (i-1) 
			decluster(&H, i);
		}
		
		// free the x and y ranges of the bins
		free(B.xRange);
		free(B.yRange);
		
		// free the bins
		free(B.array);
		
		//getch();
	}
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	// return the execution time in seconds of our gp algorithm
	return seconds;	// successful return of ourPlacerGP 
}

float ourPlacerGPtest(nodes *nodes, nets nets, chip chip, int nmax){
	int i, j, n;	// counters for the loops
	int level;		// counter for the hypergaph levels	
	hypergraph H;	// hypergraph
	double sum;		// sum to calculate average chip width
	connectivitySortedNodes sortedNodes;	// sorted nodes based on their connectivity
	double *x, *y;	// temporary vectors for x and y coordinates
	cg_stats Stats;	// cg stats	
	int val;		// returned value of cg minimization
	int iterations;	// how many times we will minimize the W
	double WNew, WOld;	// new and old value of W
	double *xOld, *yOld;	// coordinates before minimization
	
	// start counting the execution clocks of our algorithm
	clock_t start = clock();
	
	// calculate average chip width
	sum = 0.0;
	for(i = 0; i < chip.numberOfRows; i++){
		sum += chip.array[i].width;
	}
	averageWidth = sum / chip.numberOfRows;
	
	// mixed size circuit
	createH0(&H, nodes->numberOfNodes - nodes->numberOfTerminals);
	
	// get the sorted nodes based on their connectivity	
	sortedNodes = sortNodesByConnectivity(*nodes);
	
	// initialize level to zero
	level = 0;
	
	// create the hypergraph
	while(H.array[level].blockNumber > nmax){
		level++;	// increase the number of levels
		firstChoiceClustering(&H, level);	// create the new hypergraph level, by dividing the previous level
	}
	
	for(i = 0; i < H.numberOfLevels; i++){
		printf("H%d\tNumber Of Blocks: %ld\t Blocks ID: %ld - %ld\n", i, H.array[i].blockNumber, H.array[i].range[0], H.array[i].range[1]);
	}
	
	// solve qp
	for(i = 0; i < H.array[level].blockNumber; i++){
		// set the nodes id in H level
		nodes->array[sortedNodes.array[H.array[level].range[0] + i]].id = i;
	}

	// test if we can minimize the (H level) nodes with linear method
	// if we can solve QP, save the solutions 
	// else we keep the random coordinates as initial solution
	QP(&(*nodes), nets, H, sortedNodes, 30);
	
	// check if node coordinate its out of chip
	for(i = 0; i < nodes->numberOfNodes - nodes->numberOfTerminals; i++){
		// check x
		if(nodes->array[i].x < 0.0){
			nodes->array[i].x = 0.0;
		}else if(nodes->array[i].x > chip.array[0].width){
			nodes->array[i].x = chip.array[0].width;
		}
		
		// check y
		if(nodes->array[i].y < 0.0){
			nodes->array[i].y = 0.0;
		}else if(nodes->array[i].y > chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height){
			nodes->array[i].y = chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height;
		}
	}
	
	// minimize for all nodes with non-linear method
	// set the number of nodes
	n = 2 * (nodes->numberOfNodes - nodes->numberOfTerminals);
	
	// create the result array
	x = malloc(n * sizeof(double));

	// set starting guess 
	for(i = 0; i < n; i++){
		if(i < (n / 2)){
			if((x[i] < 0.0) || (x[i] > chip.array[0].width)){
				x[i] = 0.0;
				//printf("x out of chip: %lf\n", x[i]);
			}else{
				// x starting guess
				x[i] = nodes->array[i].x;
				//x[i] = (double) chip.array[0].width / 2;
				//x[i] = 1.;
				
			}
		}else{
			if((x[i] < 0.0) || (x[i] > chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height)){
				//printf("y out of chip: %lf\n", x[i]);
				x[i] = 0.0;
			}else{
				// y starting guess
				x[i] = nodes->array[i - (n / 2)].y;
				//x[i] = (double) (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height) / 2;
				//x[i] = 1.;
			}
		}
	}
	
	// set minimization loops
	// depends on the number of nodes of each circuit
	if(H.numberOfLevels == 1){
		// coordinate fund already for all nodes with linear minimization
		iterations = 0;
	}else{
		switch(nodes->numberOfNodes){
			
			case 10000 ... 30000:
				iterations = 15;
				break;
					
			case 30001 ... 80000:
				iterations = 20;
				break;
				
			case 80001 ... 170000:
				iterations = 25;
				break;
				
			case 170001 ... 220000:
				iterations = 30;
				break;
				
			default:
				iterations = 25;
				break;
		}
	}
	
	for(i = 0; i < iterations; i++){
		// calculate W
		WOld = W(*nodes, nets, chip);
		
		// create the result array
		xOld = malloc((nodes->numberOfNodes - nodes->numberOfTerminals) * sizeof(double));
		yOld = malloc((nodes->numberOfNodes - nodes->numberOfTerminals) * sizeof(double));
		
		for(j = 0; j < nodes->numberOfNodes - nodes->numberOfTerminals; j++){
			xOld[j] = nodes->array[j].x;
			yOld[j] = nodes->array[j].y;
		}
		
		// minimize the W
		val = cg_descent(x, n, &Stats, NULL, 1.e-2, myvalue, mygrad, NULL, NULL, *(&nodes), nets, chip);
		printf("W minimized. Returned value: %d\n\n", val);
		
		// calculate W
		WNew = W(*nodes, nets, chip);
		
		// if after minimization we have worst results we break the process
		if(WOld <= WNew){
			for(j = 0; j < nodes->numberOfNodes - nodes->numberOfTerminals; j++){
				nodes->array[j].x = xOld[j];
				nodes->array[j].y = yOld[j];
			}
			
			break;
		}
		
		// free the two temporary arrays
		free(xOld);
		free(yOld);
	}
	
	// check if node coordinate its out of chip
	for(i = 0; i < nodes->numberOfNodes - nodes->numberOfTerminals; i++){
		// check x
		if(nodes->array[i].x < 0.0){
			nodes->array[i].x = 0.0;
		}else if(nodes->array[i].x > chip.array[0].width){
			nodes->array[i].x = chip.array[0].width;
		}
		
		// check y
		if(nodes->array[i].y < 0.0){
			nodes->array[i].y = 0.0;
		}else if(nodes->array[i].y > chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height){
			nodes->array[i].y = chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height;
		}
	}
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	// return the execution time in seconds of our gp algorithm
	return seconds;	// successful return of ourPlacerGP 
}

double myvalue(nodes *nodes, nets nets, chip chip, double *x, CG_INT n){
	CG_INT i;
	double f;

	for(i = 0; i < n; i++){
		if(i < (n / 2)){
			if((x[i] < 0.0) || (x[i] > chip.array[0].width)){
				//printf("x out of chip: %lf\n", x[i]);
			}else{
				nodes->array[i].x = x[i];
			}
		}else{
			if((x[i] < 0.0) || (x[i] > chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height)){
				//printf("y out of chip: %lf\n", x[i]);
			}else{
				nodes->array[i - (n / 2)].y = x[i];
			}
		}
	}
	
	// calculate the value of W
	f = W(*nodes, nets, chip);
	printf("f : %lf\n", f);

	// return the value of wirelength function
	return 	f;// successful return of myvalue 
}

void mygrad(nodes *nodes, nets nets, chip chip, double *g, double *x, CG_INT n){
	double t;
	CG_INT i;
	double *gX, *gY;
	
	gX = malloc(nodes->numberOfNodes * sizeof(double));
	gY = malloc(nodes->numberOfNodes * sizeof(double));
	
	for(i = 0; i < n; i++){
		if(i < (n / 2)){
			if((x[i] < 0.0) || (x[i] > chip.array[0].width)){
				//printf("x out of chip: %lf\n", x[i]);
			}else{
				nodes->array[i].x = x[i];
			}
		}else{
			if((x[i] < 0.0) || (x[i] > chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height)){
				//printf("y out of chip: %lf\n", x[i]);
			}else{
				nodes->array[i - (n / 2)].y = x[i];
			}
		}
	}

	// calculate the gradient for each x
	wGradientX(*nodes, nets, chip, gX);

	// calculate the gradient for each y
	wGradientY(*nodes, nets, chip, gY);
	
	// create the complete g vector
	for(i = 0; i < n; i++){
		if(i < (n / 2)){
			g[i] = gX[i];
		}else{
			g[i] = gY[i - (n / 2)];
		}
	}
	
	return;	// successful return of mygrad 
}
