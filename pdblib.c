/*
 * 
 * 
 * 
 */
#include "pdblib.h"
#include "stringlib.h"


void read_pdb(FILE *pdbFile, double **coord, int nAtoms)
/*
    read_pdb_selectatoms reads the coordinates for a single step from a pdb file.
    Additionally, only the positions of select atoms (having a one in nonHatoms array) are stored.
*/
{
	int   atom;
	int   k;
	int   junk;
	char buffer[1024];
	char posChar[10];
	atom = 0;

	while (fgets(buffer,1024,pdbFile)!=NULL) {

		if (strncmp(string_firstword(buffer),"ATOM",4)==0) {
      				
			strncpy(posChar,buffer+30,8);
			posChar[8]='\0';
      			coord[atom][0] = atof(posChar);
			strncpy(posChar,buffer+38,8);
      			coord[atom][1] = atof(posChar);
			strncpy(posChar,buffer+46,8);
      			coord[atom][2] = atof(posChar);
			atom++;
   		}
	}	

}

