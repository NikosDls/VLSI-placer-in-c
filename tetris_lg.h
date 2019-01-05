#ifndef TETRISLG_H
#define TETRISLG_H

#include "parser.h"

// max distance for the tetris algorithm
#define MAXDISTANCE 9223372036854775807

// sorted array with nodes (by x coordinate)
typedef struct sortedNodes{
	int numberOfNodes;	// total number of sorted nodes (without terminal pads)
	int *array;			// this number indicates the id of the node
}sortedNodes;

// function prototypes
float tetrisLG(nodes *, chip);

#endif // TETRISLG_H
