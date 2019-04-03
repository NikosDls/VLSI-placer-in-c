#ifndef QP_H
#define QP_H

#include "parser.h"
#include "ourPlacer.h"

// function prototypes
float solveQP(nodes *, nets, int);
void QP(nodes *, nets, hypergraph, connectivitySortedNodes, int);

#endif // QP_H
