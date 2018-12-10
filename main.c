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
	tetrisLG(&nodes, chip);
return 0;	// successful return of main
}
