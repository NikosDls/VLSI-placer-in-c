#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "parser.h"
#include "random_gp.h"

float randomGP(nodes *nodes, chip chip){
	int i, j;	// counters for the loops
	int height;	// chip height
	int row;	// viariable sum to find the specific row
	int xRandom, yRandom;	// random coordinates for the node random global placement in the chip
	
	// seed for the "random" numbers
	srand((unsigned int) time(0));
	
	/*
	// calculate chip height
	height = 0;
	for(i = 0; i < chip.numberOfRows; i++){
		height += chip.array[i].height;	
	}
	//printf("%d\n", height);
	*/
	// calculate chip height
	height = chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height;
	//printf("%d", chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height);
	
	// start counting the execution clocks of tetris algorithm
	clock_t start = clock();
	
	// generate random cooridantes for all non terminal nodes (terminal nodes PADS, have fixed position in the chip)
	for(i = 0; i < nodes->numberOfNodes - nodes->numberOfTerminals; i++){
		// generate the random x and y coordinate
		yRandom = ((long long) (height - nodes->array[i].yLength) * rand()) / ((long long)RAND_MAX);
		
		// finding the specific row. Based on y coordinate
		row = 0;
		for(j = 0; j < chip.numberOfRows; j++){
			if(yRandom > (row + chip.array[j].height)){	// this isnt the the specific row j
				row += chip.array[j].height;
			}else{	// this is the the specific row j
				// we calculate the x coordinate, based on row j width
				xRandom = ((long long) (chip.array[j].width - nodes->array[i].xLength) * rand()) / ((long long)RAND_MAX);
				//printf("ROW j: %d\n", j);
				
				// we done so we break the loop
				break;
			}
		}
		//printf("%d\t%d\n", xRandom, yRandom);
				
		// set the coordinates for the node
		nodes->array[i].x = xRandom;
		nodes->array[i].y = yRandom;
	}
	
	// end of clocks counting
	clock_t end = clock();
	
	// calculate the time in seconds
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	// return the execution time in seconds of random gp algorithm
	return seconds;	// successful return of randomGP 
}
