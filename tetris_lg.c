#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "parser.h"
#include "tetris_lg.h"

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
	long long cost;		// displacement cost of moving node i to row j
	int bestRow;		// the row we fund the minimum cost
	int bestX, bestY;	// the coordinates of the best displacement in best row
	
	sortedNodes sortedNodes;	// sorted list of the nodes
	
	int cellHeight, chipHeight;		// temporary height to check if we have standard or mix size nodes and preplaced nodes in the chip
	int markFrom, markUntil;	
	
	int temp = 0;	// temporary variable to initialize first values, for the check
	
	/*
	// initialize that we dont have terminal node in the chip
	chip.pbInChip = 0;
	
	// initialize type of the circuit as standard-size nodes
	chip.standardCells = 1;
	
	// check if we have standard-size nodes or mix-size nodes
	// also we check if we have preplaced (terminal - fixed) nodes in the chip area
	for(i = 0; i < nodes->numberOfNodes; i++){
		// check for the circuit type (standard or mix) size
		if(nodes->array[i].terminal == 0){	// node i isnt terminal (its movable)
			if(temp == 0){
				// save the first node height to compare it with the others
				cellHeight = nodes->array[i].yLength;

				// initialize done
				temp = 1;
			}else{
				// if any node have different size from the first node, the circuit have mix-sized nodes
				if(cellHeight != nodes->array[i].yLength){
					chip.standardCells = 0;
				}
			}
		}else{	// check if the terminal node is in the chip area
			// calculate chip height
			chipHeight = chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height;
			
			if(((nodes->array[i].x >= 0) && (nodes->array[i].x < chip.array[0].width)) &&
				((nodes->array[i].y >= 0) && (nodes->array[i].y < chipHeight))
				){
					// terminal node is in the chip area
					chip.pbInChip = 1;
					
					// mark the terminal node as preplaced in the chip
					nodes->array[i].preplaced = 1;
					//printf("%s\t", nodes->array[i].name);
					//printf("%lf %lf\n", nodes->array[i].x, nodes->array[i].y);	
			}
		}
	}
	
	//printf("pb chip %d\tstandard cells: %d",chip.pbInChip, chip.standardCells);
	*/
	
	// get the sorted node list
	sortedNodes = sortNodes(*nodes);
	
	// start counting the execution clocks of tetris algorithm
	clock_t start = clock();
	if(nodes->array[0].connectivity==-1) chip.pbInChip = 1;
	// check which tetris algorithm we will execute 
	if(chip.standardCells == 1){	// normal tetris (standard-size nodes), without preplaced nodes
		if(chip.pbInChip == 0){
			printf("\nSTANDARD-SIZE CELLS WITHOUT PREPLACED CELLS\n"); 
			
			// initialize all rows indicators to 0
			for(i = 0; i < chip.numberOfRows; i++){
				chip.array[i].standardArray = 0;
			}
			
			// run tetris algorithm for all nodes
			// terminals (pins) have fixed position
			for(i = 0; i < sortedNodes.numberOfNodes; i++){
				// initialize best cost to infinite (max long long)
				best = MAXDISTANCE;
				
				// find the distance of the node displacement for each row and we hold the best one
				for(j = 0; j < chip.numberOfRows; j++){
					cost = sqrt((pow(nodes->array[sortedNodes.array[i]].x - chip.array[j].standardArray, 2)) + (pow(nodes->array[sortedNodes.array[i]].y - chip.array[j].coordinate, 2)));
					//printf("%lld\n", cost);	
					
					if(cost <= best &&
					  ((chip.array[j].standardArray + nodes->array[sortedNodes.array[i]].xLength) <= chip.array[0].width)){
						
						best = cost;
						bestRow = j;
					}
				}
				
				// move node to best row
				nodes->array[sortedNodes.array[i]].x = chip.array[bestRow].standardArray;
				nodes->array[sortedNodes.array[i]].y = chip.array[bestRow].coordinate;
				
				// mark node as placed
				nodes->array[sortedNodes.array[i]].placed = 1;
				
				// move the indicator in the best row
				chip.array[bestRow].standardArray += nodes->array[sortedNodes.array[i]].xLength;
			}
		}else{	// tetris-like (standard-size nodes), with preplaced nodes
			printf("\nSTANDARD-SIZE CELLS WITH PREPLACED CELLS\n"); 
			
			// initialize all rows indicators to 0
			for(i = 0; i < chip.numberOfRows; i++){
				chip.array[i].standardArray = 0;
			}
			
			/*
			// initialize all rows slots to available
			for(i = 0; i < chip.numberOfRows; i++){
				// creating the all rows for slot array
				chip.array[i].mixedArray = malloc(chip.array[i].height * sizeof(slot *));
		
				for(j = 0; j < chip.array[i].height; j++){
					// creating the all columns in each row for slot array
					chip.array[i].mixedArray[j] = malloc(chip.array[i].width * sizeof(slot));
					
					for(w = 0; w < chip.array[i].width; w++){
						chip.array[i].mixedArray[j][w] = available;
					}
				}
			}
			
			// if circuit have preplaced nodes in the chip, we have to mark those chip areas as notAvailable
			if(chip.pbInChip == 1){
				for(i = 0; i < nodes->numberOfNodes; i++){
					if(nodes->array[i].preplaced == 1){
						printf("\n%s size(%d,%d)\n", nodes->array[i].name, nodes->array[i].xLength, nodes->array[i].yLength);
						// set the y coordinate target area
						// x coordinate target area is the node x length
						markFrom = nodes->array[i].y;
						markUntil = markFrom + nodes->array[i].yLength;
						
						// temp height
						chipHeight = 0;
						
						// start scanning for the preplaced nodes
						for(j = 0; j < chip.numberOfRows; j++){
							for(w = 0; w < chip.array[j].height; w++){
								// preplaced node position found in the chip row j
								if((markFrom <= chipHeight) && (markUntil > chipHeight)){
									//printf("\n\nmark row %d ->\t", j);
									// starting mark the chip area of preplaced node as not available
									for(k = 0; k < nodes->array[i].xLength; k++){
										chip.array[j].mixedArray[chipHeight % chip.array[j].height][k] = notAvailable;
										//printf("(%d, %d) ", k, chipHeight % chip.array[j].height);
									}
								}
								
								// increase the height
								chipHeight++;
							}
						}
					}
				}
			}
			*/
			
			// run tetris-like algorithm for all nodes
			// terminals (pads) have fixed position
			for(i = 0; i < sortedNodes.numberOfNodes; i++){
				// initialize best cost to infinite (max long long)
				best = MAXDISTANCE;
			
				// find the distance of the node displacement for each row and we hold the best one
				for(j = 0; j < chip.numberOfRows; j++){		
					// initialize the possible coordinates
					xCoordinate = -1;
					yCoordinate	= -1;
					
					// check if node i fits in row j
					// "walk" through the slots column by column 
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){	// here is the difference between full left to right check tetris
						for(k = 0; k < chip.array[j].height; k++){
		
							found = 1;	// initialize found to 1 which means position not found yet
							
							// possible position found
							// check the boundaries too
							if((chip.array[j].mixedArray[k][w] == available) 
								&& ((k + nodes->array[sortedNodes.array[i]].yLength) <= chip.array[j].height)
								&& ((w + nodes->array[sortedNodes.array[i]].xLength) <= chip.array[j].width)){
								
								// check if node i, fits in the possible position with w,k coordinates
								for(z = 0; z < nodes->array[sortedNodes.array[i]].yLength; z++){
									for(x = 0; x < nodes->array[sortedNodes.array[i]].xLength; x++){
										// check if all slots are available which node needs to fit
										if(chip.array[j].mixedArray[k + z][w + x] == notAvailable){
											// slot isnt available and break the loop
											found = 2;
											break;
										}								
									}
									// break and the second loop
									if(found == 2){
										break;
									}
								}
								
								// node fits in the position with w,k coordinates
								if(found == 1){
									found = 0;
									// save the w,k coordinates for row j
									xCoordinate = w;
									yCoordinate	= k + chip.array[j].coordinate;
									//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
								}// node cant fit in row j
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
							//printf("bestRow %d\tbestX %d\tbestY %d\n", bestRow, bestX, bestY);
						}
					}
				}
				//printf("%d %d %d\n", bestRow, nodes->array[sortedNodes.array[i]].x, nodes->array[sortedNodes.array[i]].y);
				//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
				//printf("xCoordinate %d\tyCoordinate %d\n", bestX, bestY);
		
				// move node to best row
				nodes->array[sortedNodes.array[i]].x = bestX;
				nodes->array[sortedNodes.array[i]].y = bestY + chip.array[bestRow].coordinate;
		
				// mark node as placed
				nodes->array[sortedNodes.array[i]].placed = 1;
				
				// move the indicator in the best row
				chip.array[bestRow].standardArray += nodes->array[sortedNodes.array[i]].xLength;
				
				// mark all slots in the best row as notAvailable
				for(w = 0; w < nodes->array[sortedNodes.array[i]].xLength; w++){
					for(k = 0; k < nodes->array[sortedNodes.array[i]].yLength; k++){
						chip.array[bestRow].mixedArray[bestY + k][bestX + w] = notAvailable;
					}
				}
			}
		}
	}else{	// full left to right check tetris
		if((chip.standardCells == 0) && (chip.pbInChip == 1)){
			printf("\nMIXED-SIZE CELLS WITH PREPLACED CELLS\n"); 
		}else{	
			printf("\nMIXED-SIZE CELLS WITHOUT PREPLACED CELLS\n"); 
		}
		
		/*
		// initialize all rows slots to available
		for(i = 0; i < chip.numberOfRows; i++){
			// creating the all rows for slot array
			chip.array[i].mixedArray = malloc(chip.array[i].height * sizeof(slot *));
	
			for(j = 0; j < chip.array[i].height; j++){
				// creating the all columns in each row for slot array
				chip.array[i].mixedArray[j] = malloc(chip.array[i].width * sizeof(slot));
				
				for(w = 0; w < chip.array[i].width; w++){
					chip.array[i].mixedArray[j][w] = available;
				}
			}
		}
		
		// if circuit have preplaced nodes in the chip, we have to mark those chip areas as notAvailable
		if(chip.pbInChip == 1){
			for(i = 0; i < nodes->numberOfNodes; i++){
				if(nodes->array[i].preplaced == 1){
					printf("\n%s size(%d,%d)\n", nodes->array[i].name, nodes->array[i].xLength, nodes->array[i].yLength);
					// set the y coordinate target area
					// x coordinate target area is the node x length
					markFrom = nodes->array[i].y;
					markUntil = markFrom + nodes->array[i].yLength;
					
					// temp height
					chipHeight = 0;
						
					// start scanning for the preplaced nodes
					for(j = 0; j < chip.numberOfRows; j++){
						for(w = 0; w < chip.array[j].height; w++){
							// preplaced node position found in the chip row j
							if((markFrom <= chipHeight) && (markUntil > chipHeight)){
								//printf("\n\nmark row %d ->\t", j);
								// starting mark the chip area of preplaced node as not available
								for(k = 0; k < nodes->array[i].xLength; k++){
									chip.array[j].mixedArray[chipHeight % chip.array[j].height][k] = notAvailable;
									//printf("(%d, %d) ", k, chipHeight % chip.array[j].height);
								}
							}
								
							// increase the height
							chipHeight++;
						}
					}
				}
			}
		}
		*/
		
		// run tetris algorithm for all nodes
		// terminals (pads) have fixed position
		for(i = 0; i < sortedNodes.numberOfNodes; i++){
			// initialize best cost to infinite (max long long)
			best = MAXDISTANCE;
		
			// find the distance of the node displacement for each row and we hold the best one
			for(j = 0; j < chip.numberOfRows; j++){		
				// initialize the possible coordinates
				xCoordinate = -1;
				yCoordinate	= -1;
				
				// check if node i fits in row j
				// "walk" through the slots column by column 
				for(w = 0; w < chip.array[j].width; w++){
					for(k = 0; k < chip.array[j].height; k++){
	
						found = 1;	// initialize found to 1 which means position not found yet
						
						// possible position found
						// check the boundaries too
						if((chip.array[j].mixedArray[k][w] == available) 
							&& ((k + nodes->array[sortedNodes.array[i]].yLength) <= chip.array[j].height)
							&& ((w + nodes->array[sortedNodes.array[i]].xLength) <= chip.array[j].width)){
							
							// check if node i, fits in the possible position with w,k coordinates
							for(z = 0; z < nodes->array[sortedNodes.array[i]].yLength; z++){
								for(x = 0; x < nodes->array[sortedNodes.array[i]].xLength; x++){
									// check if all slots are available which node needs to fit
									if(chip.array[j].mixedArray[k + z][w + x] == notAvailable){
										// slot isnt available and break the loop
										found = 2;
										break;
									}								
								}
								// break and the second loop
								if(found == 2){
									break;
								}
							}
							
							// node fits in the position with w,k coordinates
							if(found == 1){
								found = 0;
								// save the w,k coordinates for row j
								xCoordinate = w;
								yCoordinate	= k + chip.array[j].coordinate;
								//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
							}// node cant fit in row j
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
						//printf("bestRow %d\tbestX %d\tbestY %d\n", bestRow, bestX, bestY);
					}
				}
			}
			//printf("%d %d %d\n", bestRow, nodes->array[sortedNodes.array[i]].x, nodes->array[sortedNodes.array[i]].y);
			//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
			//printf("xCoordinate %d\tyCoordinate %d\n", bestX, bestY);
	
			// move node to best row
			nodes->array[sortedNodes.array[i]].x = bestX;
			nodes->array[sortedNodes.array[i]].y = bestY + chip.array[bestRow].coordinate;
	
			// mark node as placed
			nodes->array[sortedNodes.array[i]].placed = 1;
			
			// mark all slots in the best row as not available
			for(w = 0; w < nodes->array[sortedNodes.array[i]].xLength; w++){
				for(k = 0; k < nodes->array[sortedNodes.array[i]].yLength; k++){
					chip.array[bestRow].mixedArray[bestY + k][bestX + w] = notAvailable;
				}
			}
		}
	}
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	//printf("\ntime: %lf seconds", seconds);
	
	/*
	for(i = 0; i < nodes->numberOfNodes; i++){
		printf("%7s %10d %10d\n", nodes->array[i].name, nodes->array[i].x, nodes->array[i].y);
	}
	*/
	
	// return the execution time in seconds of tetris lg algorithm
	return seconds;	// successful return of tetrisLG 
}
