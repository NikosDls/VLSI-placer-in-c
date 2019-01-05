#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

// header files
#include "parser.h"
#include "random_gp.h"
#include "tetris_lg.h"
#include "hpwl.h"
#include "NTUplace3.h"
#include "writeResults.h"

// executable files (.c) can be excluded if
// we compile the programm like gcc main.c (all .c files) -o parser
#include "parser.c"
#include "random_gp.c"
#include "tetris_lg.c"
#include "hpwl.c"

#include "NTUplace3.c"
#include "writeResults.c"

int main(){
	short int choice;	// user choice about placement algorithm
	
	// parse variables
	char filesFolder[64];		// folder name
	char filesNames[5][64];		// seperated all file names 
	nodes nodes;	// nodes
	nets nets;		// nets
	chip chip;		// chip
	// END OF parse variables
	
	int i;			// counter for the loop
	float GPseconds, LGseconds;	// to count gp and lg algorithm runtime
	float hpwl;		// half perimeter wirelength 
	
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
// end of parse

	// print the menu
	printf("\n1. Random global placement and tetris legalization"
		   "\n2. Quadratic global placement and tetris legalization"
	   	   "\n3. NTUplace3 placement\n\n");
	
	// read the users choice
	do{
		printf("Choice: ");
		scanf("%d", &choice);
	}while(choice < 1 || choice > 3);
	
	// do what the user wants
	switch(choice){
		case 1:	// random GP and tetris LG
			// random global placement
			GPseconds = randomGP(&nodes, chip);
			
			// tetris legalization
			LGseconds = tetrisLG(&nodes, chip);
			
			// compute the wirelegth
			hpwl = HPWL(nodes, nets);
			
			// write the results to file
			writeResults(nodes, choice, filesFolder, GPseconds, LGseconds, hpwl);
			
			break;
		
		case 2:	// QP and tetris LG
			
			break;
			
		case 3:	// NTUplace3
			
			// global placement		
			createH0(&H, nodes.numberOfNodes);	// mixed size circuit
			sortedNodes = sortNodesByConnectivity(nodes);	// get the sorted nodes based on their connectivity
			NTUplace3GP(&H, 6000);	// GP Algorithm
			
			break;
	}

	return 0;	// successful return of main
}
