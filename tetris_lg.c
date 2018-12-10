#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "placer.h"

// max distance for the tetris algorithm
#define MAXDISTANCE 9223372036854775807

// sorted array with nodes (by x coordinate)
typedef struct sortedNodes{
	int numberOfNodes;	// total number of sorted nodes (without terminal pads)
	int *array;			// this number indicates the id of the node
}sortedNodes;

sortedNodes sortNodes(nodes nodes){	// sort based on x coordinate
	int i, j;	// counters for the loops
	int temp;	// temporary variable for the swaps
	sortedNodes sortedNodes;	// sorted nodes
	
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
	
	// sort the nodes for the tetris algorithm
	for(i = sortedNodes.numberOfNodes - 1; i >= 0; i--){
		for(j = 0; j < i; j++){
			if(nodes.array[sortedNodes.array[i]].x < nodes.array[sortedNodes.array[j]].x){	// values must swap
				// swap the values
            	temp = sortedNodes.array[i];
            	sortedNodes.array[i] = sortedNodes.array[j];
            	sortedNodes.array[j] = temp;
			}
		}
	}
	
	/*
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		printf("%d %d\n", i, nodes.array[sortedNodes.array[i]].x);
	}
	*/
	
	// return the sorted nodes
	return sortedNodes;	// successful return of sortNodes
}

float tetrisLG(nodes *nodes, chip chip){
	int i, j, w, k, z, x;	// counters for the loops
	int found;		// temporary variable to find the coordinates of the possible node position in the row
	int xCoordinate, yCoordinate;	// possible position coordinates for the node
	int height = 0;	// temporary save the height up to row we are
	
	long long best;		// best displacement cost
	long long cost;		// displacement cost of moving cell i to row j
	int bestRow;		// the row we fund the minimum cost
	int bestX, bestY;	// the coordinates of the best displacement in best row
	
	sortedNodes sortedNodes;	// sorted list of the nodes
	
	// initialize all rows slots to available
	for(i = 0; i < chip.numberOfRows; i++){
		// creating the all rows for slot array
		chip.array[i].array = malloc(chip.array[i].height * sizeof(slot *));

		for(j = 0; j < chip.array[i].height; j++){
			// creating the all columns in each row for slot array
			chip.array[i].array[j] = malloc(chip.array[i].width * sizeof(slot));
			
			for(w = 0; w < chip.array[i].width; w++){
				chip.array[i].array[j][w] = available;
			}
		}
	}

	// get the sorted node list
	sortedNodes = sortNodes(*nodes);
	
	// start counting the execution clocks of tetris algorithm
	clock_t start = clock();
	
	// run tetris algorithm for all nodes
	// terminals (pads) have fixed position
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		// initialize best cost to infinite (max long long)
		best = MAXDISTANCE;
		
		// find the distance of the cell displacement for each row and we hold the best one
		for(j = 0; j < chip.numberOfRows; j++){
			// check if node i fits in row j
			// "walk" through the slots column by column 
			for(w = 0; w < chip.array[j].width; w++){
				for(k = 0; k < chip.array[j].height; k++){
					// initialize the possible coordinates
					xCoordinate = -1;
					yCoordinate	= -1;
					
					// possible position found
					// check the boundaries too
					if((chip.array[j].array[k][w] == available) 
						&& ((k + nodes->array[sortedNodes.array[i]].yLength) <= chip.array[j].height)
						&& ((w + nodes->array[sortedNodes.array[i]].xLength) <= chip.array[j].width)){
						
						found = 0;	// initialize found to 0 which means position found
						
						// check if node i, fits in the possible position with w,k coordinates
						for(z = 0; z < nodes->array[sortedNodes.array[i]].yLength; z++){
							for(x = 0; x < nodes->array[sortedNodes.array[i]].xLength; x++){
								// check if all slots are available which node needs to fit
								if(chip.array[j].array[k + z][w + x] == notAvailable){
									// slot isnt available and break the loop
									found = 1;
									break;
								}								
							}
							// break and the second loop
							if(found == 1){
								break;
							}
						}
						
						// node fits in the position with w,k coordinates
						if(found == 0){
							// save the w,k coordinates for row j
							xCoordinate = w;
							yCoordinate	= k + chip.array[j].coordinate;
							//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
						}
						// else node cant fit in row j
					}
					
					// we already found position for the node in row j so we dont have to continue the lookup
					if(found == 0){
						break;
					}
				}
				
				// break this loop aswell
				if(found == 0){
					break;
				}
			}
			
			// at this point we have the possible coordinates for the node i in row j
			// if coordinates have the value -1, that means node i doesnt fit in the row j
			// calculate the displacement cost for node i in row j
			if((xCoordinate != -1) && (yCoordinate != -1)){
				cost = sqrt((pow(nodes->array[sortedNodes.array[i]].x - xCoordinate, 2)) + (pow(nodes->array[sortedNodes.array[i]].y - yCoordinate, 2)));
				//printf("%lld\n", cost);
				
				// save the best cost
				// and best row (plus the node coordinates in the row)
				if(cost <= best){
					best = cost;
					bestRow = j;
					bestX = xCoordinate;
					bestY = yCoordinate - chip.array[j].coordinate;
				}
			}
		}
		//printf("%d %d %d\n", bestRow, nodes->array[sortedNodes.array[i]].x, nodes->array[sortedNodes.array[i]].y);
		//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
		//printf("xCoordinate %d\tyCoordinate %d\n", bestX, bestY);

		// move cell to best row
		nodes->array[sortedNodes.array[i]].x = xCoordinate;
		nodes->array[sortedNodes.array[i]].y = yCoordinate;
		
		// mark cell as placed
		nodes->array[sortedNodes.array[i]].placed = 1;
		
		// mark all slots in the best row as notAvailable
		for(w = 0; w < nodes->array[sortedNodes.array[i]].xLength; w++){
			for(k = 0; k < nodes->array[sortedNodes.array[i]].yLength; k++){
				chip.array[bestRow].array[bestY + k][bestX + w] = notAvailable;
			}
		}
	}
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	//printf("\ntime: %lf seconds", seconds);
	
	// return the execution time in seconds of tetris algorithm
	return seconds;	// successful return of tetrisLG 
}
