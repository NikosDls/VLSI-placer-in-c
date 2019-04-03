#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "parser.h"
#include "tetris_lg.h"

void swapP(int *a, int *b){
	int temp = *a; 
	*a = *b; 
	*b = temp; 
}

// the presence of equal elements requires special care in quicksort
// so we use Hoare's partition method
int partitionP(nodes nodes, int arr[], int low, int high){ 
	double pivot = nodes.array[arr[high]].priority;	// pivot
	int i = low - 1, j = high + 1;
	
	while(1){ 
		// find leftmost element greater than 
		// or equal to pivot 
		do{ 
			i++; 
		}while(nodes.array[arr[i]].priority < pivot); 
   
		// find rightmost element smaller than 
		// or equal to pivot 
		do{ 
			j--; 
        }while(nodes.array[arr[j]].priority > pivot); 

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
int partitionP(nodes nodes, int arr[], int low, int high){
	double pivot = nodes.array[arr[high]].priority;	// pivot 
	int i = (low - 1);	// index of smaller element 
	int j;
	
	for (j = low; j <= (high - 1); j++){ 
		// if current element is smaller than or equal to pivot
		if(nodes.array[arr[j]].priority <= pivot){ 
			i++;    // increment index of smaller element 
			
			// swap i and j element
			swapP(&arr[i], &arr[j]); 
		}
	}

	// swap i + 1 and high element
	swapP(&arr[i + 1], &arr[high]); 
	
	return (i + 1);
}
*/
  
void quickSortP(nodes nodes, int arr[], int low, int high){ 
	if(low < high){ 
		// pi is partitioning index, arr[p] is now 
		// at right place
		int pi = partitionP(nodes, arr, low, high); 
		
		// separately sort elements before 
		// partition and after partition 
		quickSortP(nodes, arr, low, pi - 1); 
		quickSortP(nodes, arr, pi + 1, high); 
	}
}

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

	// calculate the priority for the sorted nodes
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		nodes.array[i].priority = k1 * nodes.array[i].x - k2 * nodes.array[i].xLength - k3 * nodes.array[i].yLength;
	}

	// sort the nodes for the tetris algorithm
	quickSortP(nodes, sortedNodes.array, 0, sortedNodes.numberOfNodes - 1);

	/*
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
	
	for(i = 0; i < sortedNodes.numberOfNodes; i++){
		printf("%d %s %lf\n", i, nodes.array[sortedNodes.array[i]].name, nodes.array[sortedNodes.array[i]].priority);
		if((i % 200) == 0){
			getch();
		}
	}
	*/

	// return the sorted nodes
	return sortedNodes;	// successful return of sortNodes
}

float tetrisLG(nodes *nodes, chip chip){
	int i, j, w, k, z, x, l;// counters for the loops
	int found, found1;		// temporary variable to find the coordinates of the possible node position in the row
	int flag;				// temporary flag for the mixed sized tetris
	int xCoordinate, yCoordinate;	// possible position coordinates for the node
	
	long long best;		// best displacement cost
	long long cost;		// displacement cost of moving node i to row j
	int bestRow;		// the row we fund the minimum cost
	int bestX, bestY;	// the coordinates of the best displacement in best row
	
	sortedNodes sortedNodes;	// sorted list of the nodes
	
	int c;			// temporary counter
	
	int minimumXLength;	// minimum x length of the nodes
	int chipHeight;		// temporary height to count height of the checked area
	int markFromY, markUntilY;		// until where we have to check (in y direction)
	int markFromX, markUntilX;		// until where we have to check (in x direction)
	
	// get the sorted node list
	sortedNodes = sortNodes(*nodes);
	
	// start counting the execution clocks of tetris algorithm
	clock_t start = clock();
	
	// initialize all rows indicators to 0
	for(i = 0; i < chip.numberOfRows; i++){
		chip.array[i].standardArray = 0;
	}
	
	// initialize the minimum as the first node x length
	minimumXLength = nodes->array[0].xLength;
		
	// find the minimum x size of the nodes
	for(i = 1; i < nodes->numberOfNodes - nodes->numberOfTerminals; i++){
		if(nodes->array[i].xLength < minimumXLength){
			minimumXLength = nodes->array[i].xLength; 
		}
	}
	
	// check which tetris algorithm we will execute 
	if(chip.standardCells == 1){	// normal tetris (standard-size nodes), without preplaced nodes
		if(chip.pbInChip == 0){
			printf("\nSTANDARD-SIZE CELLS WITHOUT PREPLACED CELLS\n"); 

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
					  ((chip.array[j].standardArray + nodes->array[sortedNodes.array[i]].xLength) <= chip.array[j].width)){
						
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

					found1 = 1;	// initialize found to 1 which means there is no preplaced node after the standard array 
					c = 0;		// initialize counter to 0
					
					// check if we have preplaced block after the standard array 
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){
						for(k = 0; k < chip.array[j].height; k++){
							if(chip.array[j].mixedArray[k][w] == notAvailable){
								// increase standard array
								chip.array[j].standardArray += (c + 1);
								
								// change the value of flag found1 and re-initialize counter
								found1 = 0;
								c = 0;
								
								// break the height loop							
								break;
							}
							found1 = 1;
						}
						// increase the counter
						c++;
						
						// there is no preplaced node after the standard array
						// so we stop the search
						if((found1 == 1) && (c >= minimumXLength)){
							break;
						}
					}
					
					// check if node i fits in row j
					// "walk" through the slots column by column 
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){	// here is the difference between full left to right check tetris
						if((w + nodes->array[sortedNodes.array[i]].xLength) > chip.array[j].width){
							break;
						}
						
						for(k = 0; k < chip.array[j].height; k++){
							if((k + nodes->array[sortedNodes.array[i]].yLength) > chip.array[j].height){
								break;
							}
							
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
	}else{	// mixed sized tetris
		if((chip.standardCells == 0) && (chip.pbInChip == 1)){
			printf("\nMIXED-SIZE CELLS WITH PREPLACED CELLS\n"); 
		}else{	
			printf("\nMIXED-SIZE CELLS WITHOUT PREPLACED CELLS\n"); 
		}
		
		// run tetris algorithm for all nodes
		// terminals (pads) have fixed position
		for(i = 0; i < sortedNodes.numberOfNodes; i++){
			// initialize best cost to infinite (max long long)
			best = MAXDISTANCE;
			
			// set the flag
			// check if node have size bigger than rows (node isnt "standard")
			// flag is 1, when we have mixed size cell
			if(nodes->array[sortedNodes.array[i]].yLength != chip.array[0].height){
				flag = 1;
			}else{
				flag = 0;
			}
			
			// deal with standard nodes
			if(flag == 0){
				// find the distance of the node displacement for each row and we hold the best one
				for(j = 0; j < chip.numberOfRows; j++){
					// initialize the possible coordinates
					xCoordinate = -1;
					yCoordinate	= -1;
					
					found1 = 1;	// initialize found to 1 which means there is no preplaced node after the standard array 
					c = 0;		// initialize counter to 0
					
					// check if we have preplaced block after the standard array 
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){
						for(k = 0; k < chip.array[j].height; k++){
							if(chip.array[j].mixedArray[k][w] == notAvailable){
								// increase standard array
								chip.array[j].standardArray += (c + 1);
								
								// change the value of flag found1 and re-initialize counter
								found1 = 0;
								c = 0;
								
								// break the height loop							
								break;
							}
							found1 = 1;
						}
						// increase the counter
						c++;
						
						// there is no preplaced node after the standard array
						// so we stop the search
						if((found1 == 1) && (c >= minimumXLength)){
							break;
						}
					}
					
					// check if node i fits in row j
					// "walk" through the slots column by column
					//for(w = 0; w < chip.array[j].width; w++){		// full left to right check tetris (FULL LTR)
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){
						if((w + nodes->array[sortedNodes.array[i]].xLength) > chip.array[j].width){
							break;
						}
						
						for(k = 0; k < chip.array[j].height; k++){
							if((k + nodes->array[sortedNodes.array[i]].yLength) > chip.array[j].height){
								break;
							}
							
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
			}else{	// deal with non standard nodes
				//printf("%s\n", nodes->array[sortedNodes.array[i]].name);
				//getch();

				// find the distance of the node displacement for each row and we hold the best one
				for(j = 0; j < chip.numberOfRows; j++){
					// initialize the possible coordinates
					xCoordinate = -1;
					yCoordinate	= -1;
					
					found1 = 1;	// initialize found to 1 which means there is no preplaced node after the standard array 
					c = 0;		// initialize counter to 0
					
					// check if we have preplaced block after the standard array 
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){
						for(k = 0; k < chip.array[j].height; k++){
							if(chip.array[j].mixedArray[k][w] == notAvailable){
								// increase standard array
								chip.array[j].standardArray += (c + 1);
								
								// change the value of flag found1 and re-initialize counter
								found1 = 0;
								c = 0;
								
								// break the height loop							
								break;
							}
							found1 = 1;
						}
						// increase the counter
						c++;
						
						// there is no preplaced node after the standard array
						// so we stop the search
						if((found1 == 1) && (c >= minimumXLength)){
							break;
						}
					}

					// check if node i fits in row j
					// "walk" through the slots column by column
					//for(w = 0; w < chip.array[j].width; w++){		// full left to right check tetris (FULL LTR)
					for(w = chip.array[j].standardArray; w < chip.array[j].width; w++){
						if((w + nodes->array[sortedNodes.array[i]].xLength) > chip.array[j].width){
							break;
						}
									
						// temp height
						chipHeight = 0;
						
						for(k = 0; k < chip.array[j].height; k++){
							// mixed sized nodes have different size from chip row height
							// so we skip the chech below
							/*
							if((k + nodes->array[sortedNodes.array[i]].yLength) > chip.array[j].height){
								break;
							}
							*/
							
							// set the y coordinate target area
							//markUntilY = chip.array[j].coordinate + k + nodes->array[i].yLength;
							
							found = 1;	// initialize found to 1 which means position not found yet
							
							// possible position found
							// check the boundaries too
							if((chip.array[j].mixedArray[k][w] == available)
								&& ((k + nodes->array[sortedNodes.array[i]].yLength) <= (chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height))
								&& ((w + nodes->array[sortedNodes.array[i]].xLength) <= chip.array[j].width)){
								
								// set the y coordinate target area
								markFromY = chip.array[j].coordinate + k;
								markUntilY = markFromY + nodes->array[sortedNodes.array[i]].yLength;
								
								// temp height
								chipHeight = chip.array[j].coordinate;
								
								// initialize counter to 0
								c = 0;		
								
								for(l = j; l < chip.numberOfRows; l++){	// here is the difference between standard and mix sized tetris
									// check if node i, fits in the possible position with w,k coordinates
									for(z = 0; z < chip.array[l].height; z++){
										for(x = 0; x < nodes->array[sortedNodes.array[i]].xLength; x++){
											if((markFromY <= chipHeight) && (markUntilY > chipHeight)){									
												// check if all slots are available which node needs to fit
												if(chip.array[l].mixedArray[z][w + x] == notAvailable){
													// deal with unnecessary loops
													while(1){
														c++;
														
														if(chip.array[l].mixedArray[z][w + x + c] == available){
															break;
														}
													}
													
													// move the w pointer (outer width loop)
													w += c;
													
													// slot isnt available and break the loop
													found = 2;
													break;
												}
											}
										}
									
										// break and the second loop
										if(found == 2){
											break;
										}
										
										// increase the height
										chipHeight++;
									}
									
									// break and the third loop
									if(found == 2){
										break;
									}
									
									// mark Y for node i done
									if(chipHeight > markUntilY){
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
							
							// increase the height
							chipHeight++;
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
			}

			//printf("%d %d %d\n", bestRow, nodes->array[sortedNodes.array[i]].x, nodes->array[sortedNodes.array[i]].y);
			//printf("xCoordinate %d\tyCoordinate %d\n", xCoordinate, yCoordinate);
			//printf("xCoordinate %d\tyCoordinate %d\n", bestX, bestY);

			// move node to best row
			nodes->array[sortedNodes.array[i]].x = bestX;
			nodes->array[sortedNodes.array[i]].y = bestY + chip.array[bestRow].coordinate;
	
			// mark node as placed
			nodes->array[sortedNodes.array[i]].placed = 1;
			
			if(flag == 0){
				// move the indicator in the best row
				chip.array[bestRow].standardArray += nodes->array[sortedNodes.array[i]].xLength;
				
				// mark all slots in the best row as not available
				for(w = 0; w < nodes->array[sortedNodes.array[i]].xLength; w++){
					for(k = 0; k < nodes->array[sortedNodes.array[i]].yLength; k++){
						chip.array[bestRow].mixedArray[bestY + k][bestX + w] = notAvailable;
					}
				}
			}else{
				// set the y coordinate target area
				markFromY = (int) roundf(nodes->array[sortedNodes.array[i]].x);
				markUntilY = markFromY + nodes->array[sortedNodes.array[i]].yLength;
				
				// set the x coordinate target area
				markFromX = (int) roundf(nodes->array[sortedNodes.array[i]].x);
				markUntilX = markFromX + nodes->array[sortedNodes.array[i]].xLength;
				
				// temp height
				chipHeight = chip.array[bestRow].coordinate;
				
				// mark all slots in the best row as not available
				for(l = bestRow; l < chip.numberOfRows; l++){
					for(w = 0; w < chip.array[l].height; w++){
						if((markFromY <= chipHeight) && (markUntilY > chipHeight)){
							for(k = markFromX; k < chip.array[l].width; k++){
								if(markUntilX > k){
									chip.array[l].mixedArray[w][k] = notAvailable;
								}
								
								// mark X for node i done
								if(k >= markUntilX){
									break;
								}
							}
						}
							
						// increase the height
						chipHeight++;
					}
					
					// mark Y for node i done
					if(chipHeight > markUntilY){
						break;
					}
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
		printf("%7s %10f %10f\n", nodes->array[i].name, nodes->array[i].x, nodes->array[i].y);
	}
	*/
	
	// return the execution time in seconds of tetris lg algorithm
	return seconds;	// successful return of tetrisLG 
}
