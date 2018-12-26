#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "placer.h"
#include "NTUplace3.h"

#include "parser.c"
#include "random_gp.c"
#include "tetris_lg.c"
#include "NTUplace3.c"

int main(){
	// parse variables
	char filesFolder[64];		// folder name
	char filesNames[5][64];		// seperated all file names 
	nodes nodes;	// nodes
	nets nets;		// nets
	chip chip;		// chip
	// END OF parse variables
	
	int i;			// counter for the loop
	float seconds;	// to count algorithms runtime
	
	FILE *fp;		// file pointer
	char fname[64];	// output file name
	
	// NTUplace3 global placer 
	hypergraph H;	// hypergraph
	connectivitySortedNodes sortedNodes;	// sorted nodes based on their connectivity
	// END OF NTUplace3 global placer
	
	// read the auxiliary file (.aux)
	readAux(filesFolder, filesNames);
	
	// check the circuit name and then
	// set useless lines depends on what circuit we work
	if(strstr(filesFolder, "adaptec") != NULL){
		UL = ADAPTEC;
	}else if(strstr(filesFolder, "bigblue") != NULL){
		UL = BIGBLUE;
	}else if(strstr(filesFolder, "ibm") != NULL){
		UL = IBM;
	}else if(strstr(filesFolder, "newblue") != NULL){
		UL = NEWBLUE;
	}else{
		UL = SUPERBLUE;
	}

// parse
	// read the nodes file (.nodes)
	readNodes(filesNames[0], &nodes);

	// read the nets file (.nets)
	readNets(filesNames[1], &nets, &nodes);
	
	// read the chip file (.scl)
	readChip(filesNames[4], &chip);
	
	// read the pads (pins) coordinates file (.pl)
	readPads(filesNames[3], &nodes);

	// NTUplace3 global placement		
	createH0(&H, nodes.numberOfNodes);	// mixed size circuit
	sortedNodes = sortNodesByConnectivity(nodes);	// get the sorted nodes based on their connectivity
	NTUplace3GP(&H, 6000);	// GP Algorithm
	
	
/*
// global placement
	// random global placement
	randomGP(&nodes, chip);

// legalization
	// tetris algorithm
	seconds = tetrisLG(&nodes, chip);
	
	
// create output file
// read the user input
	printf("\nGive the file name to print the results: ");
	gets(fname);
	
	// clear the buffer
	fflush(stdin);
	
	// open the output file
	fp = fopen(fname, "w");
	
	// write first lines with comments and infos
	fprintf(fp, "# ISPD VLSI PLACEMENT BENCHMARK: %s\n", filesFolder);
	fprintf(fp, "# PLACMENT WITH RANDOM GP AND TETRIS LG ALGORITHM\n");
	fprintf(fp, "# TETRIS RUNTIME           : %lf seconds\n", seconds);
	//fprintf(fp, "# TOTAL WIRE LENGTH (HPWL) : %lld\n\n", totalWireLength);
	
	// write results to file
	for(i = 0; i < nodes.numberOfNodes; i++){
		if(nodes.array[i].terminal == 0){	// node isnt terminal	
			fprintf(fp, "%7s %10d %10d\n", nodes.array[i].name, nodes.array[i].x, nodes.array[i].y);
		}else{	// node is terminal
			fprintf(fp, "%7s %10d %10d / FIXED\n", nodes.array[i].name, nodes.array[i].x, nodes.array[i].y);
		}
	}
	
	// all done and we close the file
	fclose(fp);	
*/	
return 0;	// successful return of main
}
