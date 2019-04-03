#ifndef ourPlacer_H
#define ourPlacer_H

#include "parser.h"
#include "cg_nl/cg_descent.h"

// user-specified target density for each bin
#define tdensity 0.5

// find maximum between two numbers
#define max(x,y) ((x) >= (y)) ? (x) : (y)

// hypergraph level 
typedef struct level{
	int blockNumber;	// total number of nodes in the level
	int range[2];		// first number indicates the first node and the second the last one in the level
}level;

// hypergraph
typedef struct hypergraph{
	int numberOfLevels;	// total number of levels
	level *array;		// array with all levels
}hypergraph;

// each bin
typedef struct bin{
	double Pb;	// area of preplaced nodes in bin
	double Db;	// area of movable nodes in bin
	double Mb;	// maximum allowable area of movable nodes in bin
	double xCenter;	// x center coordinate of the bin
	double yCenter;	// y center coordinate of the bin
}bin;

// all bin grids
typedef struct binGrids{
	int numberOfBins;	// total number of bins
	bin **array;		// array with all bins
	double wb;			// width of the bins
	double hb;			// height of the bins
	double *xRange;		// range of the bins in x directions
	double *yRange;		// range of the bins in y directions
}binGrids;

// sorted array with nodes (by connectivity)
typedef struct connectivitySortedNodes{
	int numberOfNodes;	// total number of sorted nodes
	int *array;			// this number indicates the id of the node
}connectivitySortedNodes;

// our global placer function prototypes
void createH0(hypergraph *, int);
void firstChoiceClustering(hypergraph *, int);
connectivitySortedNodes sortNodesByConnectivity(nodes);
void decluster(hypergraph *, int);

void swapC(int *, int *);
int partitionC(nodes, int [], int, int);
void quickSortC(nodes, int [], int, int);

void calculateMbPb(int, int, binGrids *, chip);
void calculateDb(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes);
void DbGradientX(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, double *);
void DbGradientY(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, double *);

double W(nodes, nets, chip);
void wGradientX(nodes, nets, chip, double *);
void wGradientY(nodes, nets, chip, double *);

float ourPlacerGP(nodes, nets, chip, int, int);

double myvalue(nodes *, nets, chip, double *, CG_INT);
void mygrad(nodes *, nets, chip, double *, double *, CG_INT);
float ourPlacerGPtest(nodes *, nets, chip, int);

#endif // ourPlacer_H
