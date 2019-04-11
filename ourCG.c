#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parser.h"
#include "ourPlacer.h"
#include "ourCG.h"

void ourCG(double (*W)(nodes, nets, chip),						// W(x,y)
			void (*wGradientX)(nodes, nets, chip, double *),	// dW(x,y) / dx
			void (*wGradientY)(nodes, nets, chip, double *),	// dW(x,y) / dy
			void (*calculateDb)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes),					// Db(x,y)
			void (*DbGradientX)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),	// dDb(x,y) / dx
			void (*DbGradientY)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),	// dDb(x,y) / dy
			double s,										// step size
			double l,										// l constant
			int horizontally, int vertically, binGrids B,	// bins informations
			nodes *nodes, nets nets, chip chip,				// chip informations
			hypergraph H,									// hypergraph H
			int i,											// level of minimization
			connectivitySortedNodes sortedNodes				// sorted nodes
			){

	int n;		// number of nodes in Hi level
	int j;		// counter for the loop
	double *objGX, *objGY;	// seperated x and y gradients of objective function
	double *gk, *dk, *gk1, *dk1;	// gradient vectors gk and gk-1, direction vectors dk and dk-1
	double bk;	// Polak-Ribiere parameter
	double *gk_gk1, gkTmulgk_gk1;	// gk - gk-1 and gk^T * (gk - gk-1)
	double gk1mulgk1;				// gk-1 * gk-1 (||gk-1||^2)
	double *bkmuldk1;		// bk * dk-1
	double *bkmuldk1_gk;	// -gk + bk * dk-1
	double ak;	// step size
	double fxk, fxk1;	// values of fxk and fxk-1
	double sum;	// temporary sum
	
	// set the number of variables of objective function
	// which is two times the number of nodes in Hi level (cuz each node have x and y variables)
	n = 2 * H.array[i].blockNumber;
	
	// creating the gradient and direction vectors
	objGX = malloc((n / 2) * sizeof(double));
	objGY = malloc((n / 2) * sizeof(double));
	gk = malloc(n * sizeof(double));
	dk = malloc(n * sizeof(double));
	gk1 = malloc(n * sizeof(double));
	dk1 = malloc(n * sizeof(double));
	
	// create all other vectors which we want for the caclulations
	gk_gk1 = malloc(n * sizeof(double));
	bkmuldk1 = malloc(n * sizeof(double));
	bkmuldk1_gk = malloc(n * sizeof(double));
	
	// initialize gradient and direction vectors to zero
	// g0 and d0 = 0
	for(j = 0; j < n; j++){
		gk1[j] = dk1[j] = 0.0;
	}
	
	//printf("%lf", W(nodes, nets, chip));getch();
	
	// run our conjugate gradient with dynamic step size 
	do{
		printf("Minimizing %lf\n", W(*nodes, nets, chip));
		
		// compute the gradient directions
		dObjF(wGradientX, wGradientY, DbGradientX, DbGradientY, l, horizontally, vertically, B, *nodes, nets, chip, H, i, sortedNodes, objGX, objGY);
		
		// merge the two gradient vectors into one
		// the 0...n-1 elements will be the x gradient values
		// and n...2*n-1 elements will be the y gradient values
		for(j = 0; j < n; j++){
			if(j < (n / 2)){
				gk[j] = objGX[j];
			}else{
				gk[j] = objGY[j - (n / 2)];
			}
		}
		
		//printf("\n\ncompute the gradient directions DONE!!\n\n");
		//getch();
		
		// compute the Polak-Ribiere parameter
		// calculate gk - gk-1
		vectorSub(gk, gk1, gk_gk1, n);
		
		// calculate gk^T * (gk - gk-1)
		gkTmulgk_gk1 = vectorMul(gk, gk_gk1, n);
		
		// calculate gk-1 * gk-1 (||gk-1||^2)
		gk1mulgk1 = vectorMul(gk1, gk1, n);
		
		if(gk1mulgk1 == 0.0){
			bk = 0.0;
		}else{
			bk = (double) gkTmulgk_gk1 / gk1mulgk1;
		}
		
		//printf("compute the Polak-Ribiere parameter DONE!!\n\n");
		//getch();
		
		// compute the conjugate directions
		// calculate bk * dk-1
		scalarMul(bk, dk1, bkmuldk1, n);
	
		// calculate dk = -gk + bk * dk-1
		vectorSub(bkmuldk1, gk, dk, n);
		
		//printf("compute the conjugate directions DONE!!\n\n");
		//getch();
		
		// compute the step size
		// calculate the denominator for ak
		sum = 0.0;
		for(j = 0; j < n; j++){
			sum += pow(dk[j], 2);
		}
		//printf("%lf", sum);
		sum = sqrt(sum);

		ak = s / sum;
		//printf("%.50lf\n", ak);
		//printf("compute the step size DONE!!\n\n");
		
		// calculate value of objective function before updating the solution
		fxk1 = objF(W, l, horizontally, vertically, B, *nodes, nets, chip, H, i, sortedNodes);
		
		// update the solution 
		for(j = 0; j < (n / 2); j++){
			nodes->array[sortedNodes.array[H.array[i].range[0] + j]].xCenter += ak * dk[j];
			//nodes->array[sortedNodes.array[H.array[i].range[0] + j]].x = nodes->array[sortedNodes.array[H.array[i].range[0] + j]].xCenter - nodes->array[sortedNodes.array[H.array[i].range[0] + j]].xLength / 2;
		
			nodes->array[sortedNodes.array[H.array[i].range[0] + j]].yCenter += ak * dk[j + (n / 2)];
			//nodes->array[sortedNodes.array[H.array[i].range[0] + j]].y = nodes->array[sortedNodes.array[H.array[i].range[0] + j]].yCenter - nodes->array[sortedNodes.array[H.array[i].range[0] + j]].yLength / 2;	
		}
		
		// calculate Db(x,y) for the new coordinates
		calculateDb(horizontally, vertically, &B, H, i, *nodes, sortedNodes);
		
		// calculate the new value of objective
		fxk = objF(W, l, horizontally, vertically, B, *nodes, nets, chip, H, i, sortedNodes);
		
		//printf("\n\nupdate the solution  DONE!!\n\n");
		//getch();

		// set gk and dk as gk-1 and gk-1 for the new loop
		for(j = 0; j < n; j++){
			gk1[j] = gk[j];
			dk1[j] = dk[j];
		}
	// until the new value of f is bigger than the previous one
	}while(fxk < fxk1);
	
	//printf("%lf", W(nodes, nets, chip));getch();
	
	return;	// successful return of ourCG
}

double objF(double (*W)(nodes, nets, chip),					// W(x,y)
			double l,										// l constant
			int horizontally, int vertically, binGrids B,	// bins informations
			nodes nodes, nets nets, chip chip,				// chip informations
			hypergraph H,									// hypergraph H
			int i,											// level of minimization
			connectivitySortedNodes sortedNodes				// sorted nodes
			){
	
	int j, w;	// counters for the loops
	double f;	// value of objective function
	
	// initialize the value of f to 0
	f = 0.0;
	
	// calculate the sum of (Db(x,y) - Mb)^2 for all bins
	for(j = 0; j < vertically; j++){
		for(w = 0; w < horizontally; w++){
			f += pow((B.array[j][w].Db - B.array[j][w].Mb), 2);
		}
	}
	
	// calculate l * sum((Db(x,y) - Mb)^2)
	f *= l;
	
	// calculate W(x,y) + l * sum((Db(x,y) - Mb)^2)
	f += W(nodes, nets, chip);
	
	// return the value of objective function
	return f;	// successful return of objF
}

void dObjF(void (*wGradientX)(nodes, nets, chip, double *),		// dW(x,y) / dx
			void (*wGradientY)(nodes, nets, chip, double *),	// dW(x,y) / dy
			void (*DbGradientX)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),	// dDb(x,y) / dx
			void (*DbGradientY)(int, int, binGrids *, hypergraph, int, nodes, connectivitySortedNodes, int, double **),	// dDb(x,y) / dy
			double l,										// l constant
			int horizontally, int vertically, binGrids B,	// bins informations
			nodes nodes, nets nets, chip chip,				// chip informations
			hypergraph H,									// hypergraph H
			int i,											// level of minimization
			connectivitySortedNodes sortedNodes,			// sorted nodes
			double *gradientX,								// x gradient of objective function
			double *gradientY								// y gradient of objective function
			){
	
	int n;				// number of nodes in Hi level
	int j, w, k;		// counter for the loop
	double *wGXt, *wGYt;// temporary partial derivative vectors of function W
	double *wGX, *wGY;	// partial derivative vectors of function W
	double **dGX, **dGY;// partial derivative vectors of function Db
	double sumX, sumY;	// to calculate sum(Db)
	
	// set the number of gradients of objective function
	// we seperate the x and y gradients in two vectors 
	// and each vector have the size equals to the number of nodes in Hi level
	n = H.array[i].blockNumber;
	
	// creating temporary the partial derivative vectors
	wGXt = malloc(nodes.numberOfNodes * sizeof(double));
	wGYt = malloc(nodes.numberOfNodes * sizeof(double));
	
	// calculate the partial derivatives of function W
	wGradientX(nodes, nets, chip, wGXt);
	wGradientY(nodes, nets, chip, wGYt);
	
	// creating the partial derivative vectors
	wGX = malloc(n * sizeof(double));
	wGY = malloc(n * sizeof(double));
	
	// keep only derivatives of Hlevel i
	for(j = 0; j < n; j++){
		wGX[j] = wGXt[sortedNodes.array[H.array[i].range[0] + j]];
		wGY[j] = wGYt[sortedNodes.array[H.array[i].range[0] + j]];
	}
	
	// we dont need the temporary partial derivatives
	// so we free them
	free(wGXt);
	free(wGYt);
	
	// creating the partial derivative vectors
	// creating the partial derivative vectors
	dGX = malloc(vertically * sizeof(double *));
	dGY = malloc(vertically * sizeof(double *));
	
	for(j = 0; j < vertically; j++){
		dGX[j] = malloc(horizontally * sizeof(double));
		dGY[j] = malloc(horizontally * sizeof(double));
	}
	
	// calculate the partial derivatives of function Db
	// and calculate denominator to calculate the l constant
	for(k = 0; k < n; k++){
		// x partial derivative for node k and all bins
		DbGradientX(horizontally, vertically, &B, H, i, nodes, sortedNodes, k, dGX);
		
		// y partial derivative for node k and all bins
		DbGradientY(horizontally, vertically, &B, H, i, nodes, sortedNodes, k, dGY);
		
		// initialize sums
		sumX = sumY = 0.0;
		
		// calculate sum((Db(x,y) / dx|y) * (Db(x,y) - Mb))
		for(j = 0; j < vertically; j++){
			for(w = 0; w < horizontally; w++){
				sumX += dGX[j][w] * (B.array[j][w].Db - B.array[j][w].Mb);
				sumY += dGY[j][w] * (B.array[j][w].Db - B.array[j][w].Mb);
			}
		}
		
		// calculate 2 * l * sum((Db(x,y) / dx|y) * (Db(x,y) - Mb))
		sumX *= 2.0 * l;
		sumY *= 2.0 * l;
		
		// calculate (W(x,y) + l * sum((Db(x,y) - Mb)^2)) / dx
		gradientX[k] = wGX[k] + sumX;
		
		// calculate  (W(x,y) + l * sum((Db(x,y) - Mb)^2)) / dy
		gradientY[k] = wGY[k] + sumY;
	}
	
	return;	// successful return of dObjF
}

void vectorAdd(double *a, double *b, double *c, int n){
	int i;	// counter for the loop
	
	// add the vectors
	for(i = 0; i < n; i++){
		c[i] = a[i] + b[i];
	}
	
	return;	// successful return of vectorAdd
}

void vectorSub(double *a, double *b, double *c, int n){
	int i;	// counter for the loop
	
	// subtract the vectors
	for(i = 0; i < n; i++){
		c[i] = a[i] - b[i];
	}
	
	return;	// successful return of vectorSub
}

double vectorMul(double *a, double *b, int n){
	int i;	// counter for the loop
	double sum = 0.0;	// result of vectors multiplication
	
	// multiplicate the vectors
	for(i = 0; i < n; i++){
		sum += a[i] * b[i];
	}
	
	// return the sum of two vectors multiplication
	return sum;	// successful return of vectorMul
}

void scalarMul(double a, double *b, double *c, int n){
	int i;	// counter for the loop
	
	// multiplicate the vector with constant a
	for(i = 0; i < n; i++){
		c[i] = a * b[i];
	}
	
	return;	// successful return of scalarMul
}
