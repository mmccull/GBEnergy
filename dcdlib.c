/*
 * a subroutine to read the dcd header from an input trajectory and write it to an output trajectory.
 * the subtourinte requires the input file and output file pointers (files must be opened prior to call)
 * The number of atoms and number of steps from the input trajectory are returned.
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int dim=3;

extern double gridmin[3];
extern double gridmax[3];

void copy_dcd_header(FILE *in, FILE *out, int *nAtoms, int *nSteps)
/*
    read_dcd_header reads the header of a dcd file returning the number of atoms and number of steps
*/
{
	int magicnumber;
	char hdr[4];
	int junk;
	int i;
	int ntitle;
	char title[80];
	float junkf;

	// magic number
	fread(&magicnumber,sizeof(int),1,in);
	fwrite(&magicnumber,sizeof(int),1,out);
	// header
	fread(hdr,4*sizeof(char),1,in);
	fwrite(hdr,4*sizeof(char),1,out);
	// number of steps
	fread(nSteps,sizeof(int),1,in);
	fwrite(nSteps,sizeof(int),1,out);
	// a bunch of crap integers I dont care about
	for(i=0;i<8;i++) {
		fread(&junk,sizeof(int),1,in);
		fwrite(&junk,sizeof(int),1,out);
	}
	// a junk float I dont care about
	fread(&junkf,sizeof(float),1,in);
	fwrite(&junkf,sizeof(float),1,out);
	// a bunch of crap integers I dont care about
	for(i=0;i<12;i++) {
		fread(&junk,sizeof(int),1,in);
		fwrite(&junk,sizeof(int),1,out);
	}
	// the number of title strings
	fread(&ntitle,sizeof(int),1,in);
	fwrite(&ntitle,sizeof(int),1,out);
	// the title
	for(i=0;i<ntitle;i++) {
		fread(title,80*sizeof(char),1,in);
		fwrite(title,80*sizeof(char),1,out);
	}
	// a bunch of crap integers I dont care about
	for(i=0;i<2;i++) {
		fread(&junk,sizeof(int),1,in);
		fwrite(&junk,sizeof(int),1,out);
	}
	// number of atoms
	fread(nAtoms,sizeof(int),1,in);
	printf("Number of atoms: %d\n",*nAtoms);
	fwrite(nAtoms,sizeof(int),1,out);
	// magic number
	fread(&magicnumber,sizeof(int),1,in);
	fwrite(&magicnumber,sizeof(int),1,out);
	//  the end of the header

}

void read_dcd_header(FILE *in, int *nAtoms, int *nSteps)
/*
    read_dcd_header reads the header of a dcd file returning the number of atoms and number of steps
*/
{
	int magicnumber;
	char hdr[4];
	int junk;
	int i;
	int ntitle;
	char title[80];
	float junkf;

	// magic number
	fread(&magicnumber,sizeof(int),1,in);
	// header
	fread(hdr,4*sizeof(char),1,in);
	// number of steps
	fread(nSteps,sizeof(int),1,in);
	// a bunch of crap integers I dont care about
	for(i=0;i<8;i++) {
		fread(&junk,sizeof(int),1,in);
	}
	// a junk float I dont care about
	fread(&junkf,sizeof(float),1,in);
	// a bunch of crap integers I dont care about
	for(i=0;i<12;i++) {
		fread(&junk,sizeof(int),1,in);
	}
	// the number of title strings
	fread(&ntitle,sizeof(int),1,in);
	// the title
	for(i=0;i<ntitle;i++) {
		fread(title,80*sizeof(char),1,in);
	}
	// a bunch of crap integers I dont care about
	for(i=0;i<2;i++) {
		fread(&junk,sizeof(int),1,in);
	}
	// number of atoms
	fread(nAtoms,sizeof(int),1,in);
	// magic number
	fread(&magicnumber,sizeof(int),1,in);
	//  the end of the header

}



/* 
 * a subroutine to read the coordinates from a dcd file
 * this subroutine takes the file pointer, file position and number of atoms as input.
 * the filled out axis and coordinate arrays are returned.
 */

void read_dcd_step(FILE *d, long int *fileposition, double *axis, double **coord, int nAtoms, int box)
/*
    read_dcd_coord reads the coordinates and box information from a dcd file starting from a certain fileposition.
    The resulting fileposition is returned (thus fileposition must be sent as a pointer).
*/
{
	int   atom;
	int   k;
	int   junk;
	float temp;
	int step;


	if (box==1) {
		// skip the first input of this data section
		fseek(d,*fileposition+1*sizeof(int),0);
		// read box size informtaion
		fread(&axis[0],sizeof(double),1,d);
		fseek(d,ftell(d)+1*sizeof(double),0); // skip one double
		fread(&axis[1],sizeof(double),1,d); 
		fseek(d,ftell(d)+2*sizeof(double),0); // skip two doubles
		fread(&axis[2],sizeof(double),1,d);
       	 	// skip the last line of this data section
		fseek(d,ftell(d)+sizeof(int),0);
	} 
	// read coordinates
	for(k=0;k<dim;k++) {
	   	fread(&junk,sizeof(int),1,d);
   		for (atom=0;atom<nAtoms;atom++) {
      			fread(&temp,sizeof(float),1,d);
      			coord[atom][k] = (double) temp;
   		}
   		fread(&junk,sizeof(int),1,d);
	}	

	*fileposition = ftell(d);

}


/*
 * a subroutine to simply place the COM of a set of atoms at the origin
 */

void center_coord(double **coord, double *mass, int nAtoms)
/*
 * Moves the COM of the system to the origin
*/
{
   int atom;
   int k;
   double center[dim];
   double totalmass;

   totalmass=0;

   for(k=0;k<dim;k++) {
      center[k]=0;
      for(atom=0;atom<nAtoms;atom++) {
	      center[k]+=mass[atom]*coord[atom][k];
	      if (k==0) 
	         totalmass += mass[atom];
      }
   }
   for(k=0;k<dim;k++) {
      center[k]/=totalmass;
      for(atom=0;atom<nAtoms;atom++) {
	      coord[atom][k] -= center[k];
      }
   }


}


/* 
 * a subroutine to wrap the coordinates in a double array using the box dimensions from another double array
 * The number of atoms, masses, and the residue atom limits must also be given in order wrap using the COM of each residue.
 */
void wrap_coord(double *axis, double **coord, int nAtoms, int **resAtomLimits, int nres, double *mass)
/*
 * Makes sure the COM of each residue is inside the defined box
*/
{
   int atom;
   int k;
   int res;
   double totalmass;
   double center[dim];


   for(res=0;res<nres;res++) {
	   for (k=0;k<dim;k++) 
		   center[k]=0;
	   totalmass=0;
	   for (atom=resAtomLimits[res][0];atom<=resAtomLimits[res][1];atom++) {
		   totalmass+=mass[atom];
		   for (k=0;k<dim;k++) {
			   center[k] += mass[atom]*coord[atom][k];
		   }
	   }
	   for (k=0;k<dim;k++) {
		   center[k]/=totalmass;
	           if (center[k] < (-axis[k]/2.0) ) {
			   for (atom=resAtomLimits[res][0];atom<=resAtomLimits[res][1];atom++) {
				   coord[atom][k] += axis[k];
			   }
   	           } else if (center[k] > (axis[k]/2.0) ) {
			   for (atom=resAtomLimits[res][0];atom<=resAtomLimits[res][1];atom++) {
				   coord[atom][k] -= axis[k];
			   }
	           }
	   }
   }

}



/* a subroutine to write the positions from the input coord array to the file outtraj.  
 * outtraj must be opened prior to calling this routine.  The routine also requires the axis array and number of atoms as input
 */
void write_dcd_coord(FILE *outtraj, double *axis, double **coord,int nAtoms)
/*
 * write the atomic positions in array coord to a binary output file (outtraj)
 */
{
	int atom;
	int k;
	double junkd;
	float temp;
	int junk;

	junkd=0.0;
	junk = 48;
	// write the box size data
	fwrite(&junk,sizeof(int),1,outtraj);
	fwrite(&axis[0],sizeof(double),1,outtraj);
	fwrite(&junkd,sizeof(double),1,outtraj);
	fwrite(&axis[1],sizeof(double),1,outtraj);
	fwrite(&junkd,sizeof(double),1,outtraj);
	fwrite(&junkd,sizeof(double),1,outtraj);
	fwrite(&axis[2],sizeof(double),1,outtraj);
	fwrite(&junk,sizeof(int),1,outtraj);

        // write coordinates
	junk=nAtoms*4;
	for(k=0;k<dim;k++) {
           fwrite(&junk,sizeof(int),1,outtraj);
           for (atom=0;atom<nAtoms;atom++) {
	      temp = (float) coord[atom][k];
              fwrite(&temp,sizeof(float),1,outtraj);
	   }
           fwrite(&junk,sizeof(int),1,outtraj);
	}	

}



void read_dcd_COM_coord(FILE *d, long int *fileposition, double **axis, double **coord, double **avgCoord, int nAtoms, int nSteps, double *masses, int *resNumber, int nRes)
/*
    read_dcd_coord reads the coordinates and box information from a dcd file starting from a certain fileposition.
    The resulting fileposition is returned (thus fileposition must be sent as a pointer).
*/
{
	int   atom;
	int res;
	int   k;
	int   junk;
	float temp;
	int step;
	double resMass[nRes];
	double totalmass;
	double frameCOM[3];
	FILE* xyzFile;

	// Zero the average position
	for (res=0;res<nRes;res++) {
		resMass[res]=0;
		for (k=0;k<dim;k++) {
			avgCoord[res][k] = 0;
		}
	}

	totalmass=0;

//	xyzFile = fopen("original_cas.xyz","w");

	for (step=0;step<nSteps;step++) {
//		fprintf(xyzFile,"%d\n",nRes);
//		fprintf(xyzFile,"Step %d\n",step+1);

		if(step%1000==0) printf("Reading step %d from trajectory file\n",step);

		// skip the first input of this data section
		fseek(d,*fileposition+1*sizeof(int),0);
		// read box size informtaion
		fread(&axis[0][step],sizeof(double),1,d);
		fseek(d,ftell(d)+1*sizeof(double),0); // skip one double
		fread(&axis[1][step],sizeof(double),1,d); 
		fseek(d,ftell(d)+2*sizeof(double),0); // skip two doubles
		fread(&axis[2][step],sizeof(double),1,d);
       	 	// skip the last line of this data section
		fseek(d,ftell(d)+sizeof(int),0);

		frameCOM[0]=frameCOM[1]=frameCOM[2]=0;
        	// read coordinates
		for(k=0;k<dim;k++) {
        	   	fread(&junk,sizeof(int),1,d);
           		for (atom=0;atom<nAtoms;atom++) {
              			fread(&temp,sizeof(float),1,d);
	      			coord[resNumber[atom]-1+step*nRes][k] += (masses[atom] * ((double) temp));
				if (step ==0 && k==0) {
					resMass[resNumber[atom]-1] += masses[atom];
					totalmass +=masses[atom];
				}
				frameCOM[k] += masses[atom] * ((double) temp);
//				frameCOM[k] += ((double) temp);
	   		}
           		fread(&junk,sizeof(int),1,d);
		}
		// divide by residue masses
		for (k=0;k<dim;k++) {
			frameCOM[k]/=totalmass;
//			frameCOM[k]/= ((double) nAtoms);
			for (res=0;res<nRes;res++) {
				coord[res+step*nRes][k] /= resMass[res];
				coord[res+step*nRes][k] -= frameCOM[k];
				avgCoord[res][k] += coord[res+step*nRes][k];
			}
		}
/*
		for (res=0;res<nRes;res++) {
			fprintf(xyzFile,"C %8.3f%8.3f%8.3f\n",coord[res+step*nRes][0],coord[res+step*nRes][1],coord[res+step*nRes][2]);
		}
*/
		*fileposition = ftell(d);
	}

//	fclose(xyzFile);

	// Divide by number of steps to get average
	for (res=0;res<nRes;res++) {
		for (k=0;k<dim;k++) {
			avgCoord[res][k] /= ((double) nSteps);
		}
	}

}



void read_dcd_step_selectatoms(FILE *d, long int *fileposition, double *axis, double **coord, int nAtoms, int *nonHatoms, int box)
/*
    read_dcd_step_selectatoms reads the coordinates and box information for a single step from a dcd file.
    Additionally, only the positions of select atoms (having a one in nonHatoms array) are stored.
*/
{
	int   atom;
	int   k;
	int   junk;
	float temp;
	int step;
	int nonHatomCount;

	if (box==1) {
		// skip the first input of this data section
		fseek(d,*fileposition+1*sizeof(int),0);
		// read box size informtaion
		fread(&axis[0],sizeof(double),1,d);
		fseek(d,ftell(d)+1*sizeof(double),0); // skip one double
		fread(&axis[1],sizeof(double),1,d); 
		fseek(d,ftell(d)+2*sizeof(double),0); // skip two doubles
		fread(&axis[2],sizeof(double),1,d);
		// skip the last line of this data section
		fseek(d,ftell(d)+sizeof(int),0);
	}

	// read coordinates
	for(k=0;k<dim;k++) {
		nonHatomCount=0;
	   	fread(&junk,sizeof(int),1,d);
   		for (atom=0;atom<nAtoms;atom++) {
//			if (k==0) {printf("nonHatoms[%3d]=%3d\n",atom,nonHatoms[atom]);}
      			fread(&temp,sizeof(float),1,d);
			if (nonHatoms[atom]>=0) {
      				coord[nonHatomCount][k] = (double) temp;
				nonHatomCount++;
			}
   		}
   		fread(&junk,sizeof(int),1,d);
	}	

	*fileposition = ftell(d);

}

