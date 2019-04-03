#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

// header files
#include "parser.h"
#include "random_gp.h"
#include "quadraticPlacer.h"
#include "tetris_lg.h"
#include "hpwl.h"
#include "ourPlacer.h"
#include "writeResults.h"

// executable files (.c) can be excluded if
// we compile the programm like gcc main.c (all .c files) -o placer
#include "parser.c"
#include "random_gp.c"
#include "quadraticPlacer.c"
#include "tetris_lg.c"
#include "hpwl.c"
#include "ourPlacer.c"
#include "writeResults.c"
#include "cg\cg.c"

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
	
	int i, j, w, k;	// counters for the loops
	int markFromY, markUntilY;			// variables to set the chip area as not available, where the preplaced nodes is
	int markFromX, markUntilX;			// variables to set the chip area as not available, where the preplaced nodes is
	int unavailableArea = 0, total = 0;	// area of terminal nodes in the chip and total chip area
	int totalMovableArea;				// total movable area of the chip
	
	// placer variables
	int s;	// s is 1 when we do simple minimiation otherwise 0
	float GPseconds, LGseconds;	// to count gp and lg algorithms runtime
	float hpwl;		// half perimeter wirelength 
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
			
			for(j = 0; j < chip.numberOfRows; j++){
				if(((nodes.array[i].x >= 0) && (nodes.array[i].x < chip.array[j].width)) &&
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
				
				// increase the chip total area
				total++;
			}
		}
	}

	// if circuit have preplaced nodes in the chip, we have to mark those chip areas as notAvailable
	if(chip.pbInChip == 1){
		for(i = 0; i < nodes.numberOfNodes; i++){
			if(nodes.array[i].preplaced == 1){
				//printf("\n%s size(%d,%d)\n", nodes.array[i].name, nodes.array[i].xLength, nodes.array[i].yLength);
				// set the y coordinate target area
				markFromY = (int) roundf(nodes.array[i].y);
				markUntilY = markFromY + nodes.array[i].yLength;
				
				// set the x coordinate target area
				markFromX = (int) roundf(nodes.array[i].x);
				markUntilX = markFromX + nodes.array[i].xLength;
				
				// temp height
				chipHeight = 0;

				// start scanning for the preplaced nodes
				for(j = 0; j < chip.numberOfRows; j++){
					for(w = 0; w < chip.array[j].height; w++){
						// preplaced node position found in the chip row j
						if((markFromY <= chipHeight) && (markUntilY > chipHeight)){
							//printf("\n\nmark row %d ->\t", j);
							// starting mark the chip area of preplaced node as not available
							/*
							if(strcmp(nodes.array[i].name, "o1096433") == 0){
								printf("w %d\n", w);
								printf("%d %d\t", markFromX, markUntilX);
								printf("%d %d\n", markFromY, markUntilY);
								printf("%d\t", j);
								printf("row width %d - height %d\tcell %s %d\n", chip.array[j].width, chip.array[j].height, nodes.array[i].name, (int)nodes.array[i].x);
							getch();
							}
							*/
							
							// for some reason the loop below doesnt work
							// so we handle it with a simple loop and break statement
							//for(k = markFromX; k < markUntilX; k++){
							for(k = markFromX; k < chip.array[j].width; k++){
								if(markUntilX > k){
									chip.array[j].mixedArray[w][k] = notAvailable;
									/*
									if(strcmp(nodes.array[i].name, "o1096433") == 0){
										printf("%d - %d\n", w, k);
									}
									*/
									//printf("(%d, %d) ", k, chipHeight % chip.array[j].height);
								}
								
								// mark X for node i done
								if(k >= markUntilX){
									break;
								}
							}
							
							/*
							if(strcmp(nodes.array[i].name, "o1096433") == 0){
								printf("%d\t", w);
								printf("%d\n\n", chipHeight);
								getch();
							}
							*/
						}
							
						// increase the height
						chipHeight++;
					}
					
					// mark Y for node i done
					if(chipHeight > markUntilY){
						break;
					}
				}
			}
		}

		// calculate the available and unavailable area
		for(i = 0; i < chip.numberOfRows; i++){
			for(j = 0; j < chip.array[i].height; j++){
				for(w = 0; w < chip.array[i].width; w++){	
					if(chip.array[i].mixedArray[j][w] == notAvailable){
						unavailableArea++;
					}
				}
			}
		}
	}
	
	// calculate total movable area
	totalMovableArea = total - unavailableArea;
	
	// print the chip informations
	printf("\n\n----------------------------------------");
	printf("\nCHIP TOTAL AREA         : %d", total);
	printf("\nCHIP NOT AVAILABLE AREA : %d", unavailableArea);
	printf("\nCHIP AVAILABLE AREA     : %d", totalMovableArea);
	if(unavailableArea != 0){
		printf("\n\nSO THE %.2lf%% OF THE TOTAL CHIP AREA IS UNAVAILABLE", ((double) unavailableArea / total) * 100);
	}
	printf("\n----------------------------------------\n\n");
// END OF set the chip area

	// print the menu
	printf("\n1. Random global placement and tetris legalization"
		   "\n2. Quadratic global placement and tetris legalization"
	   	   "\n3. Our placement\n\n");

	// read the users choice
	do{
		fflush(stdin);
		printf("Choice: ");
		scanf("%d", &choice);
	}while(choice < 1 || choice > 3);

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
			
			/*
			// test the difference between HPWL and log-sum-exp model
			double testW = W(nodes, nets, chip);
		
			printf("HPWL: %f\n"
				   "W:    %lf\n"
				   "Difference: %lf", hpwl, testW, fabs(hpwl - testW));
				   	   
			return 1;
			*/
			/*
			// test gradient calculation
			double *gX, *gY;
			
			gX = malloc(nodes.numberOfNodes * sizeof(double));
			gY = malloc(nodes.numberOfNodes * sizeof(double));
			
			wGradientX(nodes, nets, chip, gX);
			wGradientY(nodes, nets, chip, gY);
			
			double s = 0.0;
			for(i = 0; i < nodes.numberOfNodes; i++){
				s += fabs(gX[i] + gY[i]);
			}
			printf("sum of the x and y gradient of function w: %lf\n", s);
			
			
			for(i = 0; i < nodes.numberOfNodes; i++){
				printf("X gradient for node %10s = %.40lf\n", nodes.array[i].name, gX[i]);
				
				// pause every 100 loops
				if((i %100) == 0){
					getch();
				}	
			}
			
			return 1;
			*/
			/*
			double test;
			for(i = 0; i < nodes.numberOfNodes; i++){
				test = PDxW(nodes, nets, chip, nodes.array[i].name);
				printf("PDxW for node %s = %lf\n", nodes.array[i].name, test);
			}
			
			return 1;
			*/
			
			// write the results to file
			writeResults(nodes, choice, filesFolder, GPseconds, LGseconds, hpwl);
			
			/*
			// test the execution times of the different LG algorithms
			e.array[0].connectivity = -1;						
			// tetris-like legalization
			LGseconds = tetrisLG(&e, chip);
			// write the results to file
			writeResults(e, choice, filesFolder, GPseconds, LGseconds, hpwl);
			*/
			break;
		
		case 2:	// QP and tetris LG
				// testing QP (for circuits with more than 15000 nodes, much memory required)
			// quadratic global placement
			GPseconds = solveQP(&nodes, nets, 30);

			// tetris-like legalization
			LGseconds = tetrisLG(&nodes, chip);

			// compute the wirelegth
			hpwl = HPWL(nodes, nets);
	
			// write the results to file
			writeResults(nodes, choice, filesFolder, GPseconds, LGseconds, hpwl);
	
			break;
			
		case 3:	// our placer
			s = 0;	// do minimization with density constraints
			//s = 1;	// do simple minimization
			
			if(s == 0){
				// our global placement with density constraints
				ourPlacerGP(nodes, nets, chip, 6000, totalMovableArea);
				
				return 1;
			
				// compute the wirelegth after our global non legal placement
				float after = HPWL(nodes, nets);
				
				// write the results to file
				writeResults(nodes, choice, filesFolder, GPseconds, 0.0, after);
				
				// tetris-like legalization
				LGseconds = tetrisLG(&nodes, chip);
				
				// compute the wirelegth
				hpwl = HPWL(nodes, nets);
				
				// write the final results to file
				writeResults(nodes, choice, filesFolder, GPseconds, LGseconds, hpwl);
			}else{
				// random initial positions for the nodes
				GPseconds = randomGP(&nodes, chip);
				
				// compute the wirelegth for the initial random non legal placement
				float before = HPWL(nodes, nets);
				
				// write the results to file
				writeResults(nodes, choice, filesFolder, GPseconds, 0.0, before);
				
				// our simple global placement
				GPseconds += ourPlacerGPtest(&nodes, nets, chip, 20000);
				
				// compute the wirelegth after our global non legal placement
				float after = HPWL(nodes, nets);
				
				// write the results to file
				writeResults(nodes, choice, filesFolder, GPseconds, 0.0, after);
				
				// tetris-like legalization
				LGseconds = tetrisLG(&nodes, chip);
				
				// compute the wirelegth
				hpwl = HPWL(nodes, nets);
				
				// write the final results to file
				writeResults(nodes, choice, filesFolder, GPseconds, LGseconds, hpwl);
			}
			
			
			break;
	}

	return 0;	// successful return of main
}
