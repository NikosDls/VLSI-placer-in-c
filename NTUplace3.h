#ifndef NTUplace3
#define NTUplace3

#include "placer.h"

// hypergraph level 
typedef struct level{
	long int blockNumber;	// total number of blocks in the level
	long int range[2];		// first number indicates the first node and the second the last one in the level
}level;

// hypergraph
typedef struct hypergraph{
	int numberOfLevels;	// total number of levels
	level *array;		// array with all levels
}hypergraph;

// each bin
typedef struct bin{
	float wb;	// width of the bin
	float hb;	// height of the bin
	float Pb;	// area of preplaced blocks in bin
	float Db;	// area of movable blocks in bin
	float Mb;	// maximum allowable area of movable blocks in bin
}bin;

// sorted array with nodes (by connectivity)
typedef struct connectivitySortedNodes{
	int numberOfNodes;	// total number of sorted nodes
	int *array;			// this number indicates the id of the node
}connectivitySortedNodes;

// all bin grids
typedef struct binGrids{
	int numberOfBins;	// total number of bins
	bin *array;			// array with all bins
}binGrids;

// NTUplace3 global placer function prototypes
void createH0(hypergraph *, int);
void firstChoiceClustering(hypergraph *, int);
connectivitySortedNodes sortNodesByConnectivity(nodes);
void NTUplace3GP(hypergraph *, int);

// NTUplace3 legalization placer


#endif // NTUplace3
