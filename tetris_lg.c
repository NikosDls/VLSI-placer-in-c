#include <stdio.h>
#include <stdlib.h>

#include "placer.h"

// max distance for the tetris algorithm
#define MAXDISTANCE 9223372036854775807

// sorted array with nodes (by x coordinate)
typedef struct sortedNodes{
	int numberOfNodes;	// total number of sorted nodes (without terminal pads)
	int *array;			// this number indicates the id of the node
}sortedNodes;

void sortNodes(nodes nodes){	// sort based on x coordinate
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

	return;	
}
