#ifndef TETRISLG_H
#define TETRISLG_H

#include "parser.h"

// max distance for the tetris algorithm
#define MAXDISTANCE 9223372036854775807

// constant for the priority calculation
int k1 = 1000;
int k2 = 1;
int k3 = 1;

// sorted array with nodes (by x coordinate)
typedef struct sortedNodes{
	int numberOfNodes;	// total number of sorted nodes (without terminal pads)
	int *array;			// this number indicates the id of the node
}sortedNodes;

// function prototypes
void swapC(int *, int *);
int partitionC(nodes, int [], int, int);
void quickSortC(nodes, int [], int, int);

sortedNodes sortNodes(nodes);
float tetrisLG(nodes *, chip);

#endif // TETRISLG_H
