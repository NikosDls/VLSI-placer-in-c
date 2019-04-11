#ifndef ourCG_H
#define ourCG_H

#include "parser.h"
#include "ourPlacer.h"

// our cg function prototypes
void ourCG(double (*W)(nodes, nets, chip),
			void (*wGradientX)(nodes, nets, chip, double *),
			void (*wGradientY)(nodes, nets, chip, double *),
			void (*calculateDb)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes),
			void (*DbGradientX)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),
			void (*DbGradientY)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),
			double,
			double,
			int, int, binGrids,
			nodes *, nets, chip,
			hypergraph,
			int,
			connectivitySortedNodes
			);

double objF(double (*W)(nodes, nets, chip),
			double,
			int, int, binGrids,
			nodes, nets, chip,
			hypergraph,
			int,
			connectivitySortedNodes
			);
			
void dObjF(void (*wGradientX)(nodes, nets, chip, double *),
			void (*wGradientY)(nodes, nets, chip, double *),
			void (*DbGradientX)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),
			void (*DbGradientY)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),
			double,
			int, int, binGrids,
			nodes, nets, chip,
			hypergraph,
			int,
			connectivitySortedNodes,
			double *,
			double *
			);
			
void vectorAdd(double *, double *, double *, int);
void vectorSub(double *, double *, double *, int);
double vectorMul(double *, double *, int);
void scalarMul(double, double *, double *, int);

#endif // ourCGr_H
