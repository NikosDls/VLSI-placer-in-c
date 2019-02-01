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
// we compile the programm like gcc main.c (all .c files) -o placer
#include "parser.c"
#include "random_gp.c"
#include "tetris_lg.c"
#include "hpwl.c"

#include "NTUplace3.c"
#include "writeResults.c"

int main(){
	//nodes e;
	// parse variables
	char filesFolder[64];		// folder name
	char filesNames[5][64];		// seperated all file names 
	nodes nodes;	// nodes
	nets nets;		// nets
	chip chip;		// chip
	// END OF parse variables
	
	int cellHeight, chipHeight;		// temporary height to check if we have standard or mix size cells and preplaced cells in the chip
	int temp = 0;	// temporary variable to initialize first values, for the check
	
 	int choice;		// user choice about placement algorithm
	
	int i, j, w, k;		// counter for the loop
	int markFrom, markUntil;	// variables to set where the preplaced nodes is in the chip as not available
	
	// placer variables
	float GPseconds, LGseconds;	// to count gp and lg algorithm runtime
	float hpwl;		// half perimeter wirelength 
	
	hypergraph H;	// hypergraph
	connectivitySortedNodes sortedNodes;	// sorted nodes based on their connectivity
	// END OF placer variables
	
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

// set the chip area
	// initialize that we dont have terminal node in the chip
	chip.pbInChip = 0;
	
	// initialize type of the circuit as standard-size cells
	chip.standardCells = 1;
	
	// check if we have standard-size cells or mix-size cells
	// also we check if we have preplaced (terminal - fixed) cells in the chip area
	for(i = 0; i < nodes.numberOfNodes; i++){
		// check for the circuit type (standard or mix) size
		if(nodes.array[i].terminal == 0){	// node i isnt terminal (its movable)
			if(temp == 0){
				// save the first node height to compare it with the others
				cellHeight = nodes.array[i].yLength;

				// initialize done
				temp = 1;
			}else{
				// if any node have different size from the first node, the circuit have mix-sized nodes
				if(cellHeight != nodes.array[i].yLength){
					chip.standardCells = 0;
				}
			}
		}else{	// check if the terminal node is in the chip area
			// calculate chip height
			chipHeight = chip.array[chip.numberOfRows - 1].coordinate + chip.array[chip.numberOfRows - 1].height;
			
			if(((nodes.array[i].x >= 0) && (nodes.array[i].x < chip.array[0].width)) &&
				((nodes.array[i].y >= 0) && (nodes.array[i].y < chipHeight))
				){
					// terminal node is in the chip area
					chip.pbInChip = 1;
					
					// mark the terminal node as preplaced in the chip
					nodes.array[i].preplaced = 1;
					//printf("%s\t", nodes->array[i].name);
					//printf("%lf %lf\n", nodes->array[i].x, nodes->array[i].y);	
			}
		}
	}

	// initialize all rows slots to available
	for(i = 0; i < chip.numberOfRows; i++){
		// creating the all rows for slot array
		chip.array[i].mixedArray = malloc(chip.array[i].height * sizeof(slot *));
	
		for(j = 0; j < chip.array[i].height; j++){
			// creating the all columns in each row for slot array
			chip.array[i].mixedArray[j] = malloc(chip.array[i].width * sizeof(slot));
				
			for(w = 0; w < chip.array[i].width; w++){
					chip.array[i].mixedArray[j][w] = available;
			}
		}
	}
				
	// if circuit have preplaced nodes in the chip, we have to mark those chip areas as notAvailable
	if(chip.pbInChip == 1){
		for(i = 0; i < nodes.numberOfNodes; i++){
			if(nodes.array[i].preplaced == 1){
				//printf("\n%s size(%d,%d)\n", nodes.array[i].name, nodes.array[i].xLength, nodes.array[i].yLength);
				// set the y coordinate target area
				// x coordinate target area is the node x length
				markFrom = nodes.array[i].y;
				markUntil = markFrom + nodes.array[i].yLength;
					
				// temp height
				chipHeight = 0;
					
				// start scanning for the preplaced nodes
				for(j = 0; j < chip.numberOfRows; j++){
					for(w = 0; w < chip.array[j].height; w++){
						// preplaced node position found in the chip row j
						if((markFrom <= chipHeight) && (markUntil > chipHeight)){
							//printf("\n\nmark row %d ->\t", j);
							// starting mark the chip area of preplaced node as not available
							for(k = 0; k < nodes.array[i].xLength; k++){
								chip.array[j].mixedArray[chipHeight % chip.array[j].height][k] = notAvailable;
								//printf("(%d, %d) ", k, chipHeight % chip.array[j].height);
							}
						}
							
						// increase the height
						chipHeight++;
					}
				}
			}
		}
	}
// END OF set the chip area

	// print the menu
	printf("\n1. Random global placement and tetris legalization"
		   "\n2. Quadratic global placement and tetris legalization"
	   	   "\n3. NTUplace3 placement\n\n");

	// read the users choice
	do{
		fflush(stdin);
		printf("Choice: ");
		scanf("%d", &choice);
	}while(choice < 1 || choice > 3);

	//printf("w %d h %d", chip.array[0].width, chip.numberOfRows);getch();
	// do what the user wants
	switch(choice){
		case 1:	// random GP and tetris-like LG
			// random global placement
			GPseconds = randomGP(&nodes, chip);
		
			//e = nodes;
		
			// tetris-like legalization
			LGseconds = tetrisLG(&nodes, chip);

			// compute the wirelegth
			hpwl = HPWL(nodes, nets);
			
			// write the results to file
			writeResults(nodes, choice, filesFolder, GPseconds, LGseconds, hpwl);
			
			/*
			e.array[0].connectivity= -1;						
			// tetris-like legalization
			LGseconds = tetrisLG(&e, chip);
			// write the results to file
			writeResults(e, choice, filesFolder, GPseconds, LGseconds, hpwl);
			*/
			break;
		
		case 2:	// QP and tetris LG
			solveQP(nodes, nets);	

			//createAandB(nodes, nets);
			break;
			
		case 3:	// NTUplace3
			
			// global placement		
			createH0(&H, nodes.numberOfNodes);	// mixed size circuit
			sortedNodes = sortNodesByConnectivity(nodes);	// get the sorted nodes based on their connectivity
			//NTUplace3GP(&H, 6000);	// GP Algorithm
			break;
	}

	return 0;	// successful return of main
}
