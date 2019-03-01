#include <stdio.h>	

#include "parser.h"
#include "hpwl.h"

float HPWL(nodes nodes, nets nets){
	int i, j;	// counters for the loops
	int xLeft, xRight, yDown, yUp;	// four points to find perimeter for the total wire length
	float totalWireLength = 0.0;	// total wirelength
		
	// calculate wire length for each net
	// with Half Perimeter Wire Length model
	for(i = 0; i < nets.numberOfNets; i++){
		// initialize four points
		// and taking the values of the first node of the i net
		xLeft = nodes.array[nets.array[i].netNodes[0]].x;
		xRight = nodes.array[nets.array[i].netNodes[0]].x + nodes.array[nets.array[i].netNodes[0]].xLength;
		yDown = nodes.array[nets.array[i].netNodes[0]].y;
		yUp = nodes.array[nets.array[i].netNodes[0]].y + nodes.array[nets.array[i].netNodes[0]].yLength;
		
		for(j = 1; j < nets.array[i].netDegree; j++){
			// check if we must change any point
			// compare with j node of the i net
			if(nodes.array[nets.array[i].netNodes[j]].x < xLeft){
				xLeft =  nodes.array[nets.array[i].netNodes[j]].x;
			}
			
			if((nodes.array[nets.array[i].netNodes[j]].x + nodes.array[nets.array[i].netNodes[j]].xLength) > xRight){
				xRight = nodes.array[nets.array[i].netNodes[j]].x + nodes.array[nets.array[i].netNodes[j]].xLength;
			}
			
			if(nodes.array[nets.array[i].netNodes[j]].y < yDown){
				yDown = nodes.array[nets.array[i].netNodes[j]].y;
			}
			
			if((nodes.array[nets.array[i].netNodes[j]].y + nodes.array[nets.array[i].netNodes[j]].yLength) > yUp){
				yUp = nodes.array[nets.array[i].netNodes[j]].y + nodes.array[nets.array[i].netNodes[j]].yLength;
			}
		}
		
		// calculate the new total wire length
		// add the net i wirelength to the total wirelength
		totalWireLength += ((xRight - xLeft) + (yUp - yDown));
	}
	
	// return the total wirelength
	return totalWireLength;	// successful return of HPWL 
}
