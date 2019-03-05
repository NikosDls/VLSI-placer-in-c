#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "parser.h"
#include "quadraticPlacer.h"
#include "ourPlacer.h"
#define CG_DESCENT_IMPLEMENTATION
#include "cg_nl/cg_descent.h"

#include "quadraticPlacer.c"


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
	H->array[level].range[0] = H->array[level - 1].range[0] + nodes;
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

void swap(int *a, int *b){
	int temp = *a; 
	*a = *b; 
	*b = temp; 
}

int partition(nodes nodes, int arr[], int low, int high){
	int pivot = nodes.array[arr[high]].connectivity;	// pivot 
	int i = (low - 1);	// index of smaller element 
	int j;
	
	for (j = low; j <= (high - 1); j++){ 
		// if current element is smaller than or equal to pivot
		if(nodes.array[arr[j]].connectivity <= pivot){ 
			i++;    // increment index of smaller element 
			
			// swap i and j element
			swap(&arr[i], &arr[j]); 
		}
	}

	// swap i + 1 and high element
	swap(&arr[i + 1], &arr[high]); 
	
	return (i + 1); 
} 
  
void quickSort(nodes nodes, int arr[], int low, int high){ 
	if(low < high){ 
		// pi is partitioning index, arr[p] is now 
		// at right place
		int pi = partition(nodes, arr, low, high); 
		
		// separately sort elements before 
		// partition and after partition 
		quickSort(nodes, arr, low, pi - 1); 
		quickSort(nodes, arr, pi + 1, high); 
	}
}

connectivitySortedNodes sortNodesByConnectivity(nodes nodes){	// sort based on x coordinate
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
	
	// sort the nodes based on their connectivity
	quickSort(nodes, sortedNodes.array, 0, sortedNodes.numberOfNodes - 1);
	
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
		printf("%d %s %d\n", i, nodes.array[sortedNodes.array[i]].name, nodes.array[sortedNodes.array[i]].connectivity);
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

void calculateMbPb(int grid_num, binGrids *B, chip chip){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;		// counters to calculate Pb
	
	// initialize the center, Db, Mb and Pb for each bin
	for(j = 0; j < grid_num; j++){
		for(w = 0; w < grid_num; w++){
			if(w == (grid_num - 1)){
				B->array[j][w].xCenter = (double) (B->xRange[w] + B->wb) / 2;
			}else{
				B->array[j][w].xCenter = (double) B->xRange[w + 1] / 2;
			}
			
			if(j == (grid_num - 1)){
				B->array[j][w].yCenter = (double) (B->yRange[j] + B->hb) / 2;
			}else{
				B->array[j][w].yCenter = (double) B->yRange[j + 1] / 2;
			}
			
			//printf("%lf %lf\n", B.array[j][w].xCenter, B.array[j][w].yCenter);
			
			B->array[j][w].Db = 0;
			B->array[j][w].Mb = 0;
			B->array[j][w].Pb = 0;
		}
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
					if((k - 1) == B->xRange[xCounter]){
						xCounter++;	
					}
				}
			}	
				
			// move to the next bin in y direction	
			if((w - 1) == B->yRange[yCounter]){
				yCounter++;	
			}	
		}
	}
		
	// calculate Mb for each bin
	for(j = 0; j < grid_num; j++){
		for(w = 0; w < grid_num; w++){
			B->array[j][w].Mb = tdensity * ((B->wb * B->hb) - B->array[j][w].Pb);
		}
	}
	
	return;	// successful return of calculateMbPb	
}

void calculateDb(int grid_num, binGrids *B, hypergraph H, int i, nodes nodes, connectivitySortedNodes sortedNodes){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;			// counters to calculate Pb
	double dx, dy;					// center (from node) to center (from bin) distance, in x and y directions
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
		for(j = 0; j < grid_num; j++){
			for(w = 0; w < grid_num; w++){
				// calculate center to center distance in x direction
				dx = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xCenter - B->array[j][w].xCenter);
				
				if((dx >= 0) && 
				   (dx <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb))){
					
					// calculate px
					//px = 1 - ax * pow(dx, 2);	
					px = 1 - ax * dx * dx;			
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + B->wb) <= dx) &&
						   (dx <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb)))){
					
					// calculate px
					//px = bx * pow(dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb), 2);
					px = bx * (dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb)) * (dx - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) - (2 * B->wb));		
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength / 2) + 2 * B->wb))) <= dx){
					// calculate px
					px = 0;
				}
				
				// calculate center to center distance in y direction
				dy = fabs(nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yCenter - B->array[j][w].yCenter);
				
				if((dy >= 0) && 
				   (dy <= ((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb))){
	
					// calculate py
					//py = 1 - ay * pow(dy, 2);	
					py = 1 - ay * dy * dy;			
				}else if((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + B->wb) <= dy) &&
						   (dy <= (((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb)))){

					// calculate py
					//py = by * pow(dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb), 2);
					py = by * (dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb)) * (dy - (nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) - (2 * B->wb));			
				}else if(((((nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength / 2) + 2 * B->wb))) <= dy){
					// calculate py
					py = 0;
				}
				
				// calculate Db for the bin (j,w)
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

void DbGradientX(int grid_num, binGrids *B, hypergraph H, int i, nodes nodes, connectivitySortedNodes sortedNodes, double *gradient){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;			// counters to calculate Pb
	double dx, dy;					// center (from node) to center (from bin) distance, in x and y directions
	double px, py, ax, bx, ay, by;	// variables to calculate Db
	
	// initialize gradient vector to zero
	for(j = 0; j < grid_num * grid_num; j++){
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
		for(j = 0; j < grid_num; j++){
			for(w = 0; w < grid_num; w++){
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
				
				// calculate Db gradient for the bin (j,w)
				gradient[(grid_num * j) + w] += nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength
												* nodes.array[sortedNodes.array[H.array[i].range[0] + k]].yLength
												* px * py;		
			}
		}
	}

	return;	// successful return of DbGradientX
}

void DbGradientY(int grid_num, binGrids *B, hypergraph H, int i, nodes nodes, connectivitySortedNodes sortedNodes, double *gradient){
	int j, w, k;	// counter for the loops
	int xCounter, yCounter;			// counters to calculate Pb
	double dx, dy;					// center (from node) to center (from bin) distance, in x and y directions
	double px, py, ax, bx, ay, by;	// variables to calculate Db
	
	// initialize gradient vector to zero
	for(j = 0; j < grid_num * grid_num; j++){
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
		for(j = 0; j < grid_num; j++){
			for(w = 0; w < grid_num; w++){
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
				
				// calculate Db gradient for the bin (j,w)
				gradient[(grid_num * j) + w] += nodes.array[sortedNodes.array[H.array[i].range[0] + k]].xLength
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
	g = 0.01 * chip.array[0].width;

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
	g = 0.01 * chip.array[0].width;

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
	g = 0.01 * chip.array[0].width;

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
			gradient[nets.array[i].netNodes[j]] += commonPos * (double) (1 / g) * expDerivativesPos[j];
			gradient[nets.array[i].netNodes[j]] += commonNeg * (double) (1 / g) * expDerivativesNeg[j];
			
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
	g = 0.01 * chip.array[0].width;

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
			gradient[nets.array[i].netNodes[j]] += commonPos * (double) (1 / g) * expDerivativesPos[j];
			gradient[nets.array[i].netNodes[j]] += commonNeg * (double) (1 / g) * expDerivativesNeg[j];
			
			//printf("\n\n%d node %s == %.50lf\n\n", nets.array[i].netNodes[j], nodes.array[nets.array[i].netNodes[j]].name, gradient[nets.array[i].netNodes[j]]);getch();
		}
		
		// free the exponential derivatives
		free(expDerivativesPos);
		free(expDerivativesNeg);
	}
	
	return;	// successful return of gradient 
}

float ourPlacerGP(nodes nodes, nets nets, chip chip, int nmax){
	int i, j, w, k;	// counters for the loops
	int level;		// counter for the hypergaph levels	
	hypergraph H;	// hypergraph
	binGrids B;		// bins
	connectivitySortedNodes sortedNodes;	// sorted nodes based on their connectivity
	int grid_num;
	double xCenter, yCenter;	// center of the chip
	int xCounter, yCounter;		// counters to calculate Pb
	double dx, dy;				// center (from node) to center (from bin) distance, in x and y directions
	double px, py, a, b;		// variables to calculate Db	
	double l;			// l constant	 
	double *wGX, *wGY;	// partial derivative vectors of function W
	double *dGX, *dGY;	// partial derivative vectors of function Db
	double sN, sD;		// numerator and denominator sum
	
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
	
	// find the center of the chip
	xCenter = (double) chip.array[0].width / 2;
	yCenter = (double) (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height) / 2;
	//printf("center: %lf %lf\n", xCenter, yCenter);

	// initialize nodes center positions
	for(i = 0; i < nodes.numberOfNodes - nodes.numberOfTerminals; i++){
		nodes.array[i].x = nodes.array[i].xCenter = xCenter;
		nodes.array[i].y = nodes.array[i].yCenter = yCenter;
	}
	
	// solve qp
	for(i = 0; i < H.array[level].blockNumber; i++){
		// set the nodes id in H level
		nodes.array[sortedNodes.array[H.array[level].range[0] + i]].id = i;
	}

	// call QP for the H level
	// if we can solve QP, save the solutions 
	// else we keep the center coordinates as initial solution
	QP(&nodes, nets, H);
	
	// start counting the execution clocks of tetris algorithm
	clock_t start = clock();
	
	for(i = level; i >= 0; i--){
		//printf("%d\n", i);

		// calculate the number of bins vertically and horizontally
		grid_num = sqrt(H.array[i].blockNumber);
		printf("\n\nNumber of grids vertically and horizontally: %d\n", grid_num);
		
		// creating the all rows for slot array
		B.array = malloc(grid_num * sizeof(bin *));
		
		for(j = 0; j < grid_num; j++){
			// creating the all columns in each row for slot array
			B.array[j] = malloc(grid_num * sizeof(bin));
		}
		
		// calculate the width and height of each bin
		B.wb = (double) chip.array[0].width / grid_num;
		B.hb = (double) (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height) / grid_num;
		printf("Width of each bin (Wb) = %lf\nHeight of each bin(Hb) = %lf\n", B.wb, B.hb);
		
		// set the sizes of the x and y ranges
		B.xRange = malloc(grid_num * sizeof(double));
		B.yRange = malloc(grid_num * sizeof(double));
		
		// set the x and y ranges
		for(j = 0; j < grid_num; j++){
			B.xRange[j] = B.wb * j;
			B.yRange[j] = B.hb * j;
			//printf("%lf ", B.xRange[j]);
			//printf("%lf ", B.yRange[j]);
		}
		
		// calculate Mb and Pb
		calculateMbPb(grid_num, &B, chip);
		
		// calculate Db
		//calculateDb(grid_num, &B, H, i, nodes, sortedNodes);
		
		// create the partial derivative vectors
		wGX = malloc(nodes.numberOfNodes * sizeof(double));
		wGY = malloc(nodes.numberOfNodes * sizeof(double));
		dGX = malloc(grid_num * grid_num * sizeof(double));
		dGY = malloc(grid_num * grid_num * sizeof(double));
		
		// calculate the partial derivatives of function W
		wGradientX(nodes, nets, chip, wGX);
		wGradientY(nodes, nets, chip, wGY);
		
		// calculate the partial derivatives of function Db
		DbGradientX(grid_num, &B, H, i, nodes, sortedNodes, dGX);
		DbGradientY(grid_num, &B, H, i, nodes, sortedNodes, dGY);

		// initialize sums
		sN = sD = 0.0;
		
		// calculate numerator and denominator
		// to calculate the l constant
		for(j = 0; j < H.array[i].blockNumber; j++){
			// add to the numerator sum, the absolut value of the partial x and y derivative y of function W 
			sN += fabs(wGX[sortedNodes.array[H.array[i].range[0] + j]] + wGY[sortedNodes.array[H.array[i].range[0] + j]]);
		}
		
		for(j = 0; j <grid_num * grid_num; j++){
			// add to the denominator sum, the absolut value of the partial x and y derivative y of function Db
			sD += fabs(dGX[j] + dGY[j]);
		}
		
		//printf("%lf / %lf\n", sN, sD);
		//getch();
		
		// calculate l
		l = (double) sN / sD;
		
		//printf("l = %lf", l);
		//getch();
		
		// objective function minimization
		//do{
		
		// increase l by two times
		// l *= 2;
		
		// calculate overflow ratio
		// ..
		
		// until overflow ratio its 0
		//}while(1);
			
			
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

float GPtest(nodes nodes, nets nets, chip chip){
	double *x;
	CG_INT i, n;
	cg_stats Stats;
    
	// start counting the execution clocks of tetris algorithm
	clock_t start = clock();
	
	// minimize the W
	cg_descent(x, n, &Stats, NULL, 1.e-8, myvalue, mygrad, NULL, NULL);
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	// return the execution time in seconds of our gp algorithm
	return seconds;	// successful return of ourPlacerGP 
}

double myvalue(nodes nodes, nets nets, chip chip, double *x, CG_INT n){
	
	
	return W(nodes, nets, chip);
}

void mygrad(nodes nodes, nets nets, chip chip, double *g, double *x, CG_INT n){
	double t ;
	CG_INT i ;
	for(i = 0; i < n; i++){


	}
	
	return;
}
