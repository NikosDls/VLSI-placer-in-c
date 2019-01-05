#include <stdio.h>

#include "parser.h"
#include "writeResults.h"

void writeResults(nodes nodes, int choice, char *filesFolder, float GPseconds, float LGseconds, float totalWireLength){
	int i;	// counter for the loop
	FILE *fp;		// file pointer
	char fname[64];	// output file name	
	
	// clear the buffer
	fflush(stdin);
	
	// read the user input
	printf("\nGive the file name to print the results: ");
	gets(fname);
		
	// open the output file
	fp = fopen(fname, "w");
	
	// write first lines with comments and infos
	fprintf(fp, "# ISPD VLSI PLACEMENT BENCHMARK: %s\n", filesFolder);
	
	switch(choice){
		case 1:
			fprintf(fp, "# PLACMENT WITH RANDOM GP AND TETRIS LG ALGORITHM\n");
			fprintf(fp, "# RANDOM GP RUNTIME       : %lf seconds\n", GPseconds);
			fprintf(fp, "# TETRIS LG RUNTIME       : %lf seconds\n", LGseconds);
			fprintf(fp, "# TOTAL WIRELENGTH (HPWL) : %lf\n\n", totalWireLength);
			break;
		
		case 2:
			fprintf(fp, "# PLACMENT WITH QUADRATIC GP AND TETRIS LG ALGORITHM\n");
			fprintf(fp, "# QUADRATIC GP RUNTIME    : %lf seconds\n", GPseconds);
			fprintf(fp, "# TETRIS LG RUNTIME       : %lf seconds\n", LGseconds);
			fprintf(fp, "# TOTAL WIRELENGTH (HPWL) : %lf\n\n", totalWireLength);
			break;
			
		case 3:
			fprintf(fp, "# PLACMENT WITH NTUplace3 ALGORITHM\n");
			fprintf(fp, "# NTUplace3 GP RUNTIME           : %lf seconds\n", GPseconds);
			fprintf(fp, "# TETRIS LIKE LG RUNTIME         : %lf seconds\n", LGseconds);
			fprintf(fp, "# TOTAL WIRELENGTH (log-sum-exp) : %lf\n\n", totalWireLength);
			break;
	}
	
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
	
	return;	// successful return of writeResults
}
