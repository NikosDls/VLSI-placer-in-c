#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "parser.h"
#include "quadraticPlacer.h"

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
	long int nodes; // number of nodes for the next level
	
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

void NTUplace3GP(nodes nodes, nets nets, connectivitySortedNodes sortedNodes, hypergraph *H, int nmax){
	int i;			// counter for the loop
	int level;		// counter for the hypergaph levels
	
	// initialize level to zero
	level = 0;
	
	// create the hypergraph
	while(H->array[level].blockNumber > nmax){
		level++;	// increase the number of levels
		firstChoiceClustering(&(*H), level);	// create the new hypergraph level, by dividing the previous level
	}
	
	for(i = 0; i < H->numberOfLevels; i++){
		printf("H%d\tNumber Of Blocks: %ld\t Blocks ID: %ld - %ld\n", i, H->array[i].blockNumber, H->array[i].range[0], H->array[i].range[1]);
	}
	
	/*
	solveQP(&H, level);
	*/
	
	return;	// successful return of NTUplace3GP
}
