#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "placer.h"

const char space[2] = " ";	// delimit token
const char tab[3] = "\t";	// delimit token

// reading the auxiliary file (.aux)
void readAux(char filesFolder[64], char filesNames[5][64]){
	char auxFile[64];	// auxiliary file name
	char allFiles[64];	// all file names in one string
	char temp[64];		// temporary variable to create file paths
	
	int flag = 0;	// temporary flag variable
	int i;			// counter for the loop
	
	FILE *fp;		// file pointer to read the files
	
	int counter;	// counter for the tokenization 
	char *token;	// variable to keep the tokens
	
	// reading the files if they are exists
	do{
		// if flag is 1 that means user gave wrong input file name 
		if(flag == 1){
			printf("FILE DOESNT EXIST ----> %s\n\n", auxFile);
		}
		
		// folder must be in the same directory with the executable 
		printf("Give the file folder name: ");
		
		// read the folder name
		gets(filesFolder);
		
		// creating the path for the auxiliary file
		// auxiliary file is "filesFolder"\"filesFolder".aux
		strcpy(auxFile, filesFolder);
		strcat(auxFile, "\\");
		strcat(auxFile, filesFolder);
		strcat(auxFile, ".aux");
		
		// reading auxiliary file message
		printf("READING CHECK FILE (.aux): %s\n", auxFile);
		
		// pause 2 seconds
		sleep(2);
		
		// open auxFile
		fp = fopen(auxFile, "r");
		
		// flag indicates that first loop is done
		flag = 1;
	}while(fp == NULL);
	
	// read auxFile content
	fgets(allFiles, 128, fp);
	
	// reading is done and we close the file
	fclose(fp);
	
	// setting counter to 0
	counter = 0;
	
	// we seperate the auxFile content
	// get the first token (RowBasedPlacement) 
	token = strtok(allFiles, space);
   
   	// walk through other tokens
   	// the second token will be the (:), so we have to pass it
	while(token != NULL && counter < 6){
		//printf("%s\n", token);
		// taking the next token
		token = strtok(NULL, space);
		
		// increase the counter by 1
		counter++;
		
		// if counter is 1, that means the token is : and we simply continue to the next loop
		if(counter == 1){
			continue;
		}
		
		// we store each file name
		strcpy(filesNames[counter - 2], token);
		//printf("%s\n", filesNames[counter - 2]);
	}

	// create and check all file paths
	for(i = 0; i < 5; i++){
		// we create all file paths
		strcpy(temp, filesNames[i]);
		strcpy(filesNames[i], filesFolder);	
		strcat(filesNames[i], "\\");
		strcat(filesNames[i], temp);
	
		// check if file exists and we have the access to read it
		if(access(filesNames[i], F_OK) != -1){
	    	printf("FILE EXISTS ----> %s\n", filesNames[i]);
		}else{
			printf("FILE DOESNT EXIST ----> %s\n", filesNames[i]);
		}	
	}
	
return;	// successful return of readAux
}

// reading nodes file (.nodes)
void readNodes(char *fileName, nodes *nodes){
	int i; 			// counter for the loops
	
	FILE *fp;		// file pointer
	char temp[128];	// temporary string to read the file, line by line
	
	int counter = 0;	// counter for the tokenization 
	char *token;		// variable to keep the tokens

	// open nodes file
	fp = fopen(fileName, "r");
	
	// read first lines with comments and infos
	for(i = 0; i < UL; i++){
		fgets(temp, 128, fp);
		//printf("%s", temp);
	}
	
	// read the number of nodes as string
	fgets(temp, 128, fp);
	//printf("%s", temp);
	
	// setting counter to 0
	counter = 0;
	
	// we seperate the number of nodes line
	// get the first token (NumNodes) 
	token = strtok(temp, space);
	
	// walk through other tokens
	// the second token will be the (:), so we have to pass it
	while(token != NULL && counter < 2){
		//printf("%s\n", token);
		// taking the next token
		token = strtok(NULL, space);
		
		// increase the counter by 1
		counter++;
		
		// if counter is 1, that means the token is : and we simply continue to the next loop
		if(counter == 1){
			continue;
		}
		
		// convert string to integer and set the number of nodes
		nodes->numberOfNodes = atoi(token);
		//printf("%d\n", nodes->numberOfNodes);
	}

	// creating the array of nodes
	nodes->array = malloc(nodes->numberOfNodes * sizeof(node));
	
	// read the number of terminals as string
	fgets(temp, 128, fp);
	//printf("%s", temp);
	
	// resetting the counter to 0
	counter = 0;
	
	// we seperate the number of terminals line
	// get the first token (NumTerminals) 
	token = strtok(temp, space);
	
	// walk through other tokens
	// the second token will be the (:), so we have to pass it
	while(token != NULL && counter < 2){
		//printf("%s\n", token);
		// taking the next token
		token = strtok(NULL, space);
		
		// increase the counter by 1
		counter++;
		
		// if counter is 1, that means the token is : and we simply continue to the next loop
		if(counter == 1){
			continue;
		}
		
		// convert string to number and set number of terminals
		nodes->numberOfTerminals = atoi(token);
		//printf("%d\n", nodes->numberOfTerminals);
	}
	
	// read all nodes sizes line by line
	for(i = 0; i < nodes->numberOfNodes; i++){
		// read the next whole line as string
		fgets(temp, 128, fp);
		//printf("%s\n", temp);
		
		// resetting the counter to 0
		counter = 0;
		
		// we seperate each line
		// get the first token (name of the node) 
		token = strtok(temp, space);
		
		// initialize the size of the node name
		nodes->array[i].name = malloc(strlen(token) + 1);
		
		// save the node name
		strcpy(nodes->array[i].name, token);
		
		// walk through other tokens
		while(token != NULL && counter < 4){
			//printf("%s\n", token);
			// taking the next token
			token = strtok(NULL, space);
			
			// increase the counter by 1
			counter++;
			
			switch(counter){
				case 1:	// x length of the node
					// convert string to integer and set the node x length
					nodes->array[i].xLength = atoi(token);
					break;
					
				case 2:	// y length of the node
					// convert string to integer and set the node y length
					nodes->array[i].yLength = atoi(token);
					
					// mark node as non-terminal
					nodes->array[i].terminal = 0;
					break;
		/*			
				case 3:	// this node is terminal
					// mark node as terminal
					nodes->array[i].terminal = 1;
					break;
		*/
			}

			// check if node is terminal and set it
			(i <= nodes->numberOfNodes - nodes->numberOfTerminals - 1) ? (nodes->array[i].terminal = 0) : (nodes->array[i].terminal = 1);
			
			// mark each node as non-placed
			nodes->array[i].placed = 0;
		}
		// set node connectivity to zero
		nodes->array[i].connectivity = 0;
	}
	
	// reading is done and we close the file
	fclose(fp);
	
	/*
	//	print to check if all nodes readed successfully
	for(i = 0; i < nodes->numberOfNodes; i++){
		printf("%10s%10d%10d%15s\n", nodes->array[i].name, nodes->array[i].xLength, nodes->array[i].yLength, (nodes->array[i].terminal == 0) ? "NON_TERMINAL" : "TERMINAL");	
	}
	*/
return;	// successful return of readNodes
}

// reading chip file (.scl)
void readChip(char *fileName, chip *chip){
	int i, j;		// counters for the loops
	int loops;		// number of the loops to read lines without any usefull information or the empty lines
	
	FILE *fp;		// file pointer
	char temp[128];	// temporary string to read the file, line by line
	
	int counter = 0;	// counter for the tokenization 
	char *token;		// variable to keep the tokens
	
	
	// open chip file
	fp = fopen(fileName, "r");
	
	// ibm circuits have one more extra line in this file before the number of rows
	// so we check that
	loops = (strstr(fileName, "ibm") != NULL) ? (UL + 1) : UL;
	
	// read first lines with comments and infos
	for(i = 0; i < loops; i++){
		fgets(temp, 128, fp);
		//printf("%s", temp);
	}
	
	// read the number of nodes as string
	fgets(temp, 128, fp);
	//printf("%s", temp);
	
	// setting counter to 0
	counter = 0;
	
	// we seperate the number of chip rows line
	// get the first token (NumRows) 
	token = strtok(temp, space);
	
	// walk through other tokens
	// the second token will be the (:), so we have to pass it
	while(token != NULL && counter < 2){
		//printf("%s\n", token);
		// taking the next token
		token = strtok(NULL, space);

		// increase the counter by 1
		counter++;
		
		// if counter is 1, that means the token is : and we simply continue to the next loop
		if(counter == 1){
			continue;
		}
		
		// convert string to number and set number of chip rows
		chip->numberOfRows = atoi(token);
		//printf("%d\n", chip->numberOfRows);
	}
	
	// creating the array of chip rows
	chip->array = malloc(chip->numberOfRows * sizeof(row));
	
	// read the next empty line
	fgets(temp, 128, fp);
	
	for(i = 0; i < chip->numberOfRows; i++){
		// read the first information of the row (CoreRow Horizontal)
		// and then we continue with the useful informations
		fgets(temp, 128, fp);
		
		// read all informations for the row i
		// each row have 7 lines of informations
		for(j = 0; j < 7; j++){
			// read the next whole line as string
			fgets(temp, 128, fp);
			
			// setting counter to 0
			counter = 0;
			
			// we seperate the number of chip rows line
			// get the first token (each time the name of element that will read) 
			token = strtok(temp, space);
			
			// walk through other tokens
			// the second token will be the (:), so we have to pass it
			// last line contains two informations, so we have to pass more than 3 tokens
			while(token != NULL && counter < ((j == 6) ? 5 : 2)){
				//printf("%s\n", token);
				// taking the next token
				token = strtok(NULL, space);
		
				// increase the counter by 1
				counter++;
				
				// if counter is 1 or 4, that means the token is :
				// if the counter is 3, that means we are in the last line and the token is Numsites
				// and we simply continue to the next loop
				if(counter == 1 || counter == 3 || counter == 4){
					continue;
				}
				
				switch(j){
					case 0:	// coordinate
						// convert string to number and set coordinate of the row i
						chip->array[i].coordinate = atoi(token);
						//printf("%d\n", chip->array[i].coordinate);
						
						break;
					
					case 1:	// height
						// convert string to number and set height of the row i
						chip->array[i].height = atoi(token);
						//printf("%d\n", chip->array[i].height);
						
						break;
					
					case 2:	// sitewidth
						// convert string to number and set sitewidth of the row i
						chip->array[i].sitewidth = atoi(token);
						//printf("%d\n", chip->array[i].sitewidth);
					
						break;
					
					case 3:	// sitespacing
						// convert string to number and set sitespacing of the row i
						chip->array[i].sitespacing = atoi(token);
						//printf("%d\n", chip->array[i].sitespacing);
					
						break;
					
					case 4:	// siteorient
						// set siteorient of the row i
						strcpy(chip->array[i].siteorient, token);
						//printf("%s\n", chip->array[i].siteorient);
					
						break;
					
					case 5:	// sitesymmetry
						// set sitesymmetry of the row i
						strcpy(chip->array[i].sitesymmetry, token);
						//printf("%s\n", chip->array[i].sitesymmetry);
					
						break;
					
					case 6:	// subrowOrigin or numsites
						// if counter is two then we are in the first element of the last line
						if(counter == 2){
							// convert string to number and set sub-row origin of the row i
							chip->array[i].subrowOrigin = atoi(token);
							//printf("%d\n", chip->array[i].subrowOrigin);
								
						}else{	// counter here is five and we are in the last element of the last line
							// convert string to number and set numsites of the row i
							chip->array[i].width = atoi(token);
							//printf("%d\n", chip->array[i].width);
						}
						
						break;
				}
			}
		}
		
		// read the last information of the row (End)
		fgets(temp, 128, fp);
	}
	
	
return;	// successful return of readChip
}

// reading nets file (.nets)
void readNets(char *fileName, nets *nets, nodes *nodes){
	int i, j;		// counters for the loops
	
	FILE *fp;		// file pointer
	char temp[128];	// temporary string to read the file, line by line
	
	int counter = 0;	// counter for the tokenization 
	char *token;		// variable to keep the tokens
	
	int id;	// temporary node id
	int counterI, counterO;	// counters to save the number of inputs(I) and outputs(O) in the net
	
	// open nets file
	fp = fopen(fileName, "r");
	
	// read first lines with comments and infos
	for(i = 0; i < UL; i++){
		fgets(temp, 128, fp);
		//printf("%s", temp);
	}
	
	// read the number of nets as string
	fgets(temp, 128, fp);
	//printf("%s", temp);
	
	// setting counter to 0
	counter = 0;
	
	// we seperate the number of nets line
	// get the first token (NumNets) 
	token = strtok(temp, space);

	// walk through other tokens
	// the second token will be the (:), so we have to pass it
	while(token != NULL && counter < 2){
		//printf("%s\n", token);
		// taking the next token
		token = strtok(NULL, space);

		// increase the counter by 1
		counter++;
		
		// if counter is 1, that means the token is : and we simply continue to the next loop
		if(counter == 1){
			continue;
		}
		
		// convert string to number and set number of nets
		nets->numberOfNets = atoi(token);
		//printf("%d\n", nets->numberOfNets);
	}
	
	// creating the array of nets
	nets->array = malloc(nets->numberOfNets * sizeof(net));
	
	// read and pass the next line (NumPins)
	fgets(temp, 128, fp);
	
	for(i = 0; i < nets->numberOfNets; i++){
		// read the next whole line as string
		fgets(temp, 128, fp);
		
		// resetting counter to 0
		counter = 0;
		
		// we seperate the net degree line
		// get the first token (NetDegree) 
		token = strtok(temp, space);
	
		// walk through other tokens
		// the second token will be the (:), so we have to pass it
		while(token != NULL && counter < 2){
			//printf("%s\n", token);
			// taking the next token
			token = strtok(NULL, space);
	
			// increase the counter by 1
			counter++;
			
			// if counter is 1, that means the token is : and we simply continue to the next loop
			if(counter == 1){
				continue;
			}

			// convert string to number and set net i degree
			nets->array[i].netDegree = atoi(token);
			//printf("%d\n", nets->array[i].netDegree);
		}
			
		// creating the array of ids for the net
		nets->array[i].netNodes = malloc(nets->array[i].netDegree * sizeof(int));
		
		// creating the array of attributes for the net
		nets->array[i].netNodesAttributes = malloc(nets->array[i].netDegree * sizeof(char));
		
		// initialize the two counters
		counterI = counterO = 0;
		
		for(j = 0; j < nets->array[i].netDegree; j++){
			// read the next whole line as string
			fgets(temp, 128, fp);
			//printf("%s\n", temp);
			
			// we seperate each line
			// get the first token (name of the node in the net) 
			token = strtok(temp, space);
			
			// get the node j id of the net i (at this point we dont know if its terminal or not)
			id = atoi(&token[1]);

			// check if node its terminal or not
			// then we add to the net i the id (index position in nodes array) of the node j
			if(strcmp(nodes->array[id].name, token) == 0){	// node isnt terminal
				nets->array[i].netNodes[j] = id;
			}else if(strcmp(nodes->array[(nodes->numberOfNodes - nodes->numberOfTerminals - 1) + id].name, token) == 0){	// node is terminal
				nets->array[i].netNodes[j] = (nodes->numberOfNodes - nodes->numberOfTerminals - 1) + id;
			}
			
			// taking the next token (indicates if node in INPUT or OUTPOUT) 
			token = strtok(NULL, space);
			
			// increase number of inputs of outputs
			if(token[0] == 'I'){
				nets->array[i].netNodesAttributes[j] = 'I';
				counterI++;
			}else if(token[0] == 'O'){
				nets->array[i].netNodesAttributes[j] = 'O';
				counterO++;
			}else{	// what B means. BOTH ???
				nets->array[i].netNodesAttributes[j] = 'B';
				counterI++;
				counterO++;
			}
		}
		
		// now we find the connectivity for each node in net i
		for(j = 0; j < nets->array[i].netDegree; j++){
			if(nets->array[i].netNodesAttributes[j] == 'I'){
				nodes->array[nets->array[i].netNodes[j]].connectivity += counterO;
			}else if(nets->array[i].netNodesAttributes[j] == 'O'){
				nodes->array[nets->array[i].netNodes[j]].connectivity += counterI;
			}else{
				nodes->array[nets->array[i].netNodes[j]].connectivity += nets->array[i].netDegree - 1;
			}
			//printf("%s %c %d\n", nodes->array[nets->array[i].netNodes[j]].name ,nets->array[i].netNodesAttributes[j], nodes->array[nets->array[i].netNodes[j]].connectivity);
		}
		//printf("\n");getch();
	}
	
return;	// successful return of readNets
}

// reading coordinates file (.pl)
void readPads(char *fileName, nodes *nodes){
	int i; 			// counter for the loops

	FILE *fp;		// file pointer
	char temp[128];	// temporary string to read the file, line by line
	
	int counter = 0;	// counter for the tokenization 
	char *token;		// variable to keep the tokens

	// open coordinates file
	fp = fopen(fileName, "r");
	
	// ibm circuits have one more extra line in this file before the pad coordinates
	// and also ibm circuits have in the start all pads coordinates, instead after the nodes (like all other circuits) 
	// so we check that
	if((strstr(fileName, "ibm") != NULL)){
		// read first lines with comments and infos
		for(i = 0; i < UL + 1; i++){
			fgets(temp, 128, fp);
			//printf("%s", temp);
		}
	}else{	// other circuits
		// read first lines with comments and infos
		for(i = 0; i < UL + 1; i++){
			fgets(temp, 128, fp);
			//printf("%s", temp);
		}
		
		// skip all nodes until we will read pad coordinate
		for(i = 0; i < nodes->numberOfNodes - nodes->numberOfTerminals - 1; i++){
			fgets(temp, 128, fp);
			//printf("%s", temp);
		}
	}
	
	// read all pads coordinates line by line
	for(i = nodes->numberOfNodes - nodes->numberOfTerminals; i < nodes->numberOfNodes; i++){
		// read the next whole line as string
		fgets(temp, 128, fp);
		//printf("%s", temp);
	
		// resetting counter to 0
		counter = 0;
		
		// we seperate each line
		// get the first token (name of the pad) 
		token = strtok(temp, space);
	
		// walk through other tokens
		// the second token will be the x coordinate and third the y coordinate
		while(token != NULL && counter < 2){
			//printf("%s\n", token);
			// taking the next token
			token = strtok(NULL, space);
	
			// increase the counter by 1
			counter++;
			
			if(counter == 1){	// x coordinate
				nodes->array[i].x = atoi(token);
				//printf("x %d\t", nodes->array[i].x);
			}else if(counter == 2){	// y coordinate
				nodes->array[i].y = atoi(token);
				//printf("y %d\n", nodes->array[i].y);
			}
		}
	}

	// reading is done and we close the file
	fclose(fp);
	
return;	// successful return of readPads
}
