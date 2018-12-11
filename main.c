#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "placer.h"
#include "parser.c"
#include "random_gp.c"
#include "tetris_lg.c"

int main(){
	char filesFolder[64];		// folder name
	char filesNames[5][64];		// seperated all file names 
	nodes nodes;	// nodes
	nets nets;		// nets
	chip chip;		// chip
	
	float seconds;	// algorithm runtime
	
	FILE *fp;		// file pointer
	char fname[64];	// output file name
	
	int i;			// counter for the loop
	
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
	
return 0;	// successful return of main
}
