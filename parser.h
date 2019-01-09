#ifndef PLACER_H
#define PLACER_H

// lines without any usefull information or empty lines
// for each different circuit the lines varies
#define ADAPTEC 4
#define BIGBLUE 4
#define IBM 5
#define NEWBLUE 4
#define SUPERBLUE 4

// "useless lines"
short int UL;

// delimit token
extern const char space[2];	

// chip slot
typedef enum {available, notAvailable} slot;

// a simple node
typedef struct node{
	char *name;	// node name
	
	// node coordinates in the chip (bottom-left)
	float x;
	float y;
	
	// node coordinates in the chip (center)
	float xCenter;
	float yCenter;
	
	// node length. Horizontally (x) and vertically (y) 
	int xLength;
	int yLength;
	
	int connectivity;	// the number of nodes which are connected with the node
	
	int placed;		// 1 if the node its placed on the chip otherwise 0 
	int terminal;	// 1 if the node its terminal otherwise 0
}node;

// all nodes
typedef struct nodes{
	int numberOfNodes;		// total number of nodes
	int numberOfTerminals;	// number of terminal nodes
	node *array;			// array with all nodes
}nodes;

// a simple net
typedef struct net{
	int netDegree;		// the number of the nodes in the net
	int *netNodes;		// array with ids from the nodes in the net
	char *netNodesAttributes;	// save all nodes attributes (Input, Output, B??)
}net;

// all nets
typedef struct nets{
	int numberOfNets;	// total number of nets
	net *array;			// array with all nets
}nets;

// chip row
typedef struct row{
	int coordinate; 		// row coordinate
	int height;				// row height
	int sitewidth;			// row site width
	int sitespacing;		// row site spacing
	char siteorient[3];		// row site orient
	char sitesymmetry[3];	// row site symmetry
	int subrowOrigin;		// sub-row origin
	int width;				// row width (numsites)
	slot **array;			// row slots (size: height x width)
}row;

// chip
typedef struct chip{
	int numberOfRows;	// total number of rows
	row *array;			// array with all chip rows
}chip;

// function prototypes
void readAux(char [64], char [5][64]);
void readNodes(char *, nodes *);
void readChip(char *, chip *);
void readNets(char *, nets *, nodes *);
void readPads(char *, nodes *);

#endif // PLACER_H
