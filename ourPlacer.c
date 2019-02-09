#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "parser.h"
#include "quadraticPlacer.h"
#include "ourPlacer.h"

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
	
	/*
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

void ourPlacerGP(nodes nodes, nets nets, chip chip, int nmax){
	int i, j, w, k;	// counters for the loops
	int level;		// counter for the hypergaph levels	
	hypergraph H;	// hypergraph
	binGrids B;		// bins
	connectivitySortedNodes sortedNodes;	// sorted nodes based on their connectivity
	int grid_num;
	int xCounter, yCounter;	// counters to calculate Pb
	
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
		
		// initialize Db, Mb and Pb for each bin
		for(j = 0; j < grid_num; j++){
			for(w = 0; w < grid_num; w++){
				B.array[j][w].Db = 0;
				B.array[j][w].Mb = 0;
				B.array[j][w].Pb = 0;
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
							B.array[yCounter][xCounter].Pb++;	
						}
						
						// move to the next bin in x direction
						if((k - 1) == B.xRange[xCounter]){
							xCounter++;	
						}
					}
				}	
					
				// move to the next bin in y direction	
				if((w - 1) == B.yRange[yCounter]){
					yCounter++;	
				}	
			}
		}
		
		// calculate Mb for each bin
		for(j = 0; j < grid_num; j++){
			for(w = 0; w < grid_num; w++){
				B.array[j][w].Mb = tdensity * ((B.wb * B.hb) - B.array[j][w].Pb);
			}
		}
		
		// calculate Db for each bin
		
		
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
		
		getch();
	}
	
	
	/*
	solveQP(&H, level);
	*/
	
	return;	// successful return of NTUplace3GP
}
