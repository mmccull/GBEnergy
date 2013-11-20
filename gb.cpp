
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include "dcdlib.h"
#include "psflib.h"
#include "pdblib.h"
#include "stringlib.h"

using namespace std;

void shift_pos_positive(double **, int, int, double *); 
void dump_pos(double **, int, int, char *);
void matmul(double **, int, int, double **, int, int, double **);
double det_3_3(double **);
void read_cfg_data(double *,double *, double *, double *, double *,char *,char *, char *, char *,char *,char *, int *); 
void read_par_file(char **, int, char *, double *, double *,int *, int);
void make_neighbor_cells(double **, int, double, int *, int *, int *,int **);
void compute_gb_energy(double **, int, double, double, double, double, double, double *, double*, double *,int *, int *, int**, int *, FILE *, FILE *, int **);
double compute_H(double, double, double, double, double, double);
double compute_dH(double, double, double, double, double, double);
int exclude_check(int, int, int **);

int main (char* argv[]) {
	
	FILE *psfFile;
	FILE *dcdFile;
	FILE *pdbFile;
	FILE *gbFile;
	FILE *gbForceFile;
	char psfFileName[1024];
	char dcdFileName[1024];
	char pdbFileName[1024];
	char parFileName[1024];
	char gbOutFileName[1024];
	char gbForceOutFileName[1024];
	char buffer[1024];
	int nPoints;

	int i,j,k;

	double **pos;
	double *axis;
	double *charge;
	double rmsd;
	char **uniqueAtomTypes;
	int *uniqueAtomTypeNum;
	int **excludeList;
	double rCut;
	double aCut;
	double kappa;
	double epsSolute,epsSolvent;
	int box;

	int nAtoms;
	int nAtomTypes;
	int nSteps;
	int step;
	long int dcdFilePosition;

	double *atomicRadius;
	double *atomicScaling;

	int *cellAssign;
	int *cellSize;
	int **cellList;
	int nCells[3];
	int totalCells;
	double max[3];

	time_t currentTime;
	time_t previousTime;

	pdbFileName[0] = dcdFileName[0] = '\0';

	// read cfg file
	read_cfg_data(&rCut,&aCut,&kappa,&epsSolute,&epsSolvent,psfFileName,parFileName,dcdFileName,pdbFileName,gbOutFileName,gbForceOutFileName,&box);

	// read in number of positions	
	psfFile = fopen(psfFileName,"r");
	read_psf_header(psfFile,&nAtoms);
	printf("Number of atoms in psf file= %d\n",nAtoms);
	//now can allocate arrays
	uniqueAtomTypeNum = new int [nAtoms];
	charge = new double [nAtoms];
	atomicRadius = new double [nAtoms];
	atomicScaling = new double [nAtoms];
	excludeList = new int *[nAtoms];
	for (i=0;i<nAtoms;i++) {
		excludeList[i] = new int [53]; // 4+4*3+4*3*3+1
		excludeList[i][0] = 0;
	}
	//read atom type data from psf file
	uniqueAtomTypes = read_psf_atom_data(psfFile, nAtoms, charge, &nAtomTypes, uniqueAtomTypeNum, atomicRadius, atomicScaling);
//	printf("Number of unique atom types = %d\n",nAtomTypes);
	read_psf_bond_data(psfFile,excludeList);

	//allocate GB arrays
	pos = new double* [nAtoms];
	cellAssign = new int [nAtoms];
	for(i=0;i<nAtoms;i++) {
		pos[i] = new double [3];
	}

	//close psf file
	fclose(psfFile);

	//open output files
	gbFile = fopen(gbOutFileName,"w");
	gbForceFile = fopen(gbForceOutFileName,"w");

	// read coordinates 
	if (pdbFileName[0] != '\0') {
		
		pdbFile = fopen(pdbFileName, "r");
		read_pdb(pdbFile,pos,nAtoms);
		fclose(pdbFile);

		// move the positions into positive quadrant
		shift_pos_positive(pos,nAtoms,3,max);
	
		// allocate neighbor cell arrays
		totalCells=0;
		for(k=0;k<3;k++) {
			nCells[k] = (int) (max[k]/16.0)+1;
		}
		totalCells = nCells[0]*nCells[1]*nCells[2];
		cellSize = new int [totalCells];
		cellList = new int* [totalCells];
		for (i=0;i<totalCells;i++) {
			cellList[i] = new int [nAtoms];
			cellSize[i] = 0;
		}
	
		// make neighbor cells
		make_neighbor_cells(pos,nAtoms,16.0,nCells,cellAssign,cellSize,cellList);
	
		// compute GB
		compute_gb_energy(pos,nAtoms,rCut,aCut,kappa,epsSolute,epsSolvent,charge,atomicRadius,atomicScaling,nCells,cellSize,cellList,cellAssign,gbFile,gbForceFile,excludeList);
	
		delete [] cellSize;
		delete [] cellList;

	} else if (dcdFileName[0] != '\0')  {
		// read nsteps from dcd file
		dcdFile = fopen(dcdFileName,"r");
		read_dcd_header(dcdFile,&nAtoms,&nSteps);	
		printf("Number of atoms in dcd file= %d\n",nAtoms);
		
		// allocate arrays
		axis = new double [3];
		
		previousTime = clock();
		for (step=0;step<nSteps;step++) {
		
			//read in positions
			dcdFilePosition = ftell(dcdFile);
			read_dcd_step(dcdFile,&dcdFilePosition,axis,pos,nAtoms,box);
		
			// move the positions into positive quadrant
			shift_pos_positive(pos,nAtoms,3,max);
		
			// allocate neighbor cell arrays
			totalCells=0;
			for(k=0;k<3;k++) {
				nCells[k] = (int) (max[k]/rCut)+1;
			}
			totalCells = nCells[0]*nCells[1]*nCells[2];
//			printf("Total Number of Cells:%d (%d,%d,%d)\n",totalCells, nCells[0],nCells[1],nCells[2]);
			cellSize = new int [totalCells];
			cellList = new int* [totalCells];
			for (i=0;i<totalCells;i++) {
				cellList[i] = new int [nAtoms];
				cellSize[i] = 0;
			}
	
			// make neighbor cells
			make_neighbor_cells(pos,nAtoms,rCut,nCells,cellAssign,cellSize,cellList);
		
			// compute GB energy
			compute_gb_energy(pos,nAtoms,rCut,aCut,kappa,epsSolute,epsSolvent,charge,atomicRadius,atomicScaling,nCells,cellSize,cellList,cellAssign,gbFile,gbForceFile,excludeList);

			// print timing periodically	
			if(step%1==0) {
				printf("Reading step %d from trajectory file\n",step);
				// save time	
				currentTime = clock();
				printf("Step %4d to %4d took %f\n",step+1-1000, step+1,( (double) (currentTime-previousTime)) / CLOCKS_PER_SEC);
				previousTime = currentTime;
			}

			delete [] cellSize;
			delete [] cellList;
		
		}
		fclose(dcdFile);
	}
	fclose(gbFile);
	fclose(gbForceFile);


}

/*
 *                SUBROUTINES
 */


void read_cfg_data(double *rCut,double *aCut, double *kappa, double *epsSolute, double *epsSolvent,char *psfFileName,char *parFileName,char *dcdFileName, char *pdbFileName, char *gbOutFileName,char *gbForceOutFileName, int *box) {

	double C = 0.000395394; // Angstroms^2*mole/(K*L) combined constants in Debye screening length term

	char buffer[1024];
	char tempBuffer[1024];
	char check[15];
	char *firstWord;
	double T;
	double conc;

	*rCut = 0;
	*aCut = 0;

	while (fgets(buffer,1024,stdin) != NULL) {

		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
//		printf ("First word = %s\n",firstWord);
		if (strncmp(firstWord,"psfFile",7)==0) {
			strcpy(psfFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"rCut",4)==0) {
			*rCut = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"aCut",4)==0) {
			*aCut = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"temperature",11)==0) {
			T = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"ionConcentration",16)==0) {
			conc = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"epsSolvent",10)==0) {
			*epsSolvent = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"epsSolute",9)==0) {
			*epsSolute = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"box",3)==0) {
			*box = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"dcdFile",7)==0) {
			strcpy(dcdFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"pdbFile",7)==0) {
			strcpy(pdbFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"parFile",7)==0) {
			strcpy(parFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"gbFile",6)==0) {
			strcpy(gbOutFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"gbForceFile",11)==0) {
			strcpy(gbForceOutFileName,string_secondword(buffer));
		}
	

	}	

	if (*rCut==0) {
		printf("Cannot have a 0.0 cutoff value.  Using default value of 14.0 Angstroms\n");
		*rCut = 14.0;
	}
	if (*aCut==0) {
		printf("Cannot have a 0.0 GB cutoff value.  Using default value of 14.0 Angstroms\n");
		*aCut = 14.0;
	}

	*kappa = 1.0/sqrt(*epsSolvent*C*T/conc);

	printf("psf file: %s\n",psfFileName);
	printf("dcd file: %s\n",dcdFileName);
	printf("pdb file: %s\n",pdbFileName);
	printf("GB out: %s\n",gbOutFileName);
	printf("GB force out: %s\n",gbForceOutFileName);
	printf("real cutoff: %f\n",*rCut);
	printf("GB cutoff: %f\n",*aCut);
	printf("Solute dielectric: %f\n",*epsSolute);
	printf("Solvent dielectric: %f\n",*epsSolvent);
	printf("kappa: %f\n",*kappa);


}


void dump_pos(double **coord, int nAtoms, int nSteps, char *fileName) {

	int i, j, k;
	int atom;
	int step;
	FILE *xyzOut;

	xyzOut = fopen(fileName,"w");

	for (step=0;step<nSteps;step++) {

		for (atom=0;atom<nAtoms;atom++) {

			fprintf(xyzOut,"%12.6f%12.6f%12.6f\n",coord[atom+step*nAtoms][0],coord[atom+step*nAtoms][1],coord[atom+step*nAtoms][2]);

		}

	}

	fclose(xyzOut);

}

double det_3_3(double **mat) {

	double det;

	det = mat[0][0]*mat[1][1]*mat[2][2] - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0];
	det+= mat[0][2]*mat[1][0]*mat[2][1] - mat[0][2]*mat[1][1]*mat[2][0];

	return det;
}

void matmul(double **mat1, int n, int m, double **mat2, int p, int q, double **mat3) {

	int i,j,k;


	if (m!=p) {
		printf ("Error in matmul: dimensions do not match (m=%d,p=%d)\n",m,p);
		exit(1);
	}

	for(i=0;i<n;i++) {
		for(j=0;j<q;j++) {
			mat3[i][j]=0;
			for(k=0;k<m;k++) {
				mat3[i][j]+=mat1[i][k]*mat2[k][j];
			}
		}
	}

}


	// shift reference frame to be in quadrant 1 for hilbert curve
void shift_pos_positive(double **pos,int nRows, int nCols, double *max) {

	int i, j;
	double min[nCols];

	for (i=0;i<nRows;i++) {
		for (j=0;j<nCols;j++) {
			if (i==0) {
				min[j] = pos[i][j];
				max[j] = pos[i][j];
			} else if (pos[i][j] < min[j]) {
				min[j] = pos[i][j];
			} else if (pos[i][j] > max[j]) {
				max[j] = pos[i][j];
			}
		}
	}
	for (i=0;i<nRows;i++) {
		for (j=0;j<nCols;j++) {
			pos[i][j] -= min[j];
		}
	}

	for(j=0;j<nCols;j++) {
		max[j] -= min[j];
	}
}



// make neighbor cells
//
void make_neighbor_cells(double **pos, int numNonHatoms, double rCut, int *nCells, int *cellAssign, int *cellSize, int **cellList) {

	int i, j, k;
	int cell;
	int cellXYZ[3];
	int atom;

	for (atom=0;atom<numNonHatoms;atom++) {

		for (k=0;k<3;k++) {

			cellXYZ[k] = (int) (pos[atom][k]/rCut);

		}

		cell = cellXYZ[0]*nCells[1]*nCells[2]+cellXYZ[1]*nCells[2]+cellXYZ[2];

//		printf("%8.3f%8.3f%8.3f cell:%4d%4d%4d\n",pos[atom][0],pos[atom][1],pos[atom][2],cellXYZ[0],cellXYZ[1],cellXYZ[2]);

		cellAssign[atom] = cell;
		cellList[cell][cellSize[cell]]=atom;
		cellSize[cell]++;


	}

}

// compute GB energy
void compute_gb_energy(double **pos,int nAtoms, double rCut, double aCut, double kappa, double epsSolute, double epsSolvent, double *charge, double *atomicRadius, double *atomicScaling, int *nCells,int *cellSize, int **cellList, int *cellAssign, FILE *gbFile, FILE *gbForceFile,int **excludeList)
{

	double pi=3.1415926535;
	double offset=0.09; // Angstroms
//	double ke = 332.063711; //kcal/mol*Angstroms/e^2 coulombs constant
	double ke = 332.0636; //kcal/mol*Angstroms/e^2 coulombs constant
	double delta = 1.0; // parameter for born radius eqaution
	double beta = 0.8; // parameter for born radius eqaution
	double gamma = 4.85; // parameter for born radius eqaution
	double rCut2 = rCut*rCut;
	double aCut2 = aCut*aCut;

	int atom1;
	int atom2;
	int atom3;
	int i,j,k;
	int cell1;
	int cellXYZ[3];
	int cellMin[3],cellMax[3];
	int cellX1, cellY1, cellZ1;
	double dist1, dist1_2;
	double dx, dy, dz;
	double dist2, dist2_2;
	int neighborList[nAtoms][nAtoms+1];
	double psi[nAtoms];
	double bornRadius[nAtoms];
	double gbE;
	double gbEij;
	double selfE;
	double coulE;
	double Dij;
	double fij;
	double hij;
	double scale;
	double temp;
	double force[nAtoms][3];
	double dEijdrij;
	double fdEijdrij;
	double dDijdrij;
	double dfijdrij;
	double dalphaidrij[nAtoms];
	double dscaledrij;
	double qiqj;
	double aiaj;
	double dEda[nAtoms];
	double tmp_dEda;
	double tanh2;
	double daidr[nAtoms];
	double forceMag,fx,fy,fz;
	double dhij,dhji;

	coulE=0;

	/* GB Phase 1: Compute Born Radii */
	for (atom1=0;atom1<nAtoms;atom1++) {

		/* zero forces */
		force[atom1][0]=force[atom1][1]=force[atom1][2]=0.0;
		dEda[atom1]=0;
	
		/* zero first term in neighbor list.  This will be used as a counting term */
		neighborList[atom1][0]=0;

		/* Compute the neighbor cell that atom1 is in and what the surrounding cells are (no PBC)*/
		for (k=0;k<3;k++) {
			cellXYZ[k] = (int) (pos[atom1][k]/rCut);
			if (cellXYZ[k]==0) {
				cellMin[k]=0;
				if (nCells[k]==1) {
					cellMax[k] = 0;
				} else {
					cellMax[k] = 1;
				}
			} else if (cellXYZ[k] == (nCells[k]-1)) {
				cellMin[k] = nCells[k]-2;
				cellMax[k] = nCells[k]-1;
			} else {
				cellMin[k] = cellXYZ[k]-1;
				cellMax[k] = cellXYZ[k]+1;
			}
		}

		psi[atom1] = 0;
//		for (atom2=0;atom2<nAtoms;atom2++) {	
		for (cellX1=cellMin[0];cellX1<=cellMax[0];cellX1++) {
			for (cellY1=cellMin[1];cellY1<=cellMax[1];cellY1++) {
				for (cellZ1=cellMin[2];cellZ1<=cellMax[2];cellZ1++) {
					cell1 = cellX1*nCells[1]*nCells[2]+cellY1*nCells[2]+cellZ1;

					for (i=0;i<cellSize[cell1];i++) {
			
						atom2 = cellList[cell1][i];
						if (atom2 != atom1) {
							/* compute the distance between the atoms */
							dist1_2=0;
							for (k=0;k<3;k++) {
								dist1_2 += (pos[atom1][k]-pos[atom2][k])*(pos[atom1][k]-pos[atom2][k]);
							}
							dist1 = sqrt(dist1_2);
							if (dist1<rCut) {
								dist1 = sqrt(dist1_2);
		
								/* accumulate the GB screening function (also called psi) */
								psi[atom1] += compute_H(dist1,atomicRadius[atom1],atomicRadius[atom2],offset,atomicScaling[atom2],aCut);
	
								/* populate neighbor list while we are here */
								neighborList[atom1][neighborList[atom1][0]+1] = atom2;
								neighborList[atom1][0]++; // count number of atoms in atom1's neighbor list
							}
							if (dist1<rCut && atom2>atom1 && exclude_check(atom1,atom2,excludeList)==0) {
								/* MM debug */
//								printf("Atom %4d and atom %4d are a distance of %10.5f apart and not excluded\n",atom1, atom2, dist1);
								/* NAMD shift/smoothing function */
								scale = (1-dist1_2/rCut2);
								scale *= scale;
	
								/* accumulate coulomb energy */
//								coulE += ke*charge[atom1]*charge[atom2]/dist1*scale;
								temp = ke*charge[atom1]*charge[atom2]/dist1*scale;
								coulE += temp;
//								printf("%4d %4d %10.5f %10.5f %10.5f %10.5f\n",atom1, atom2, charge[atom1], charge[atom2], dist1, temp);

							} 
						}
					}
				}
			}
		}
		/* Finish screening of atom1 */
		psi[atom1] *= (atomicRadius[atom1]-offset);

		/* Compute born radius for atom1 using screening term */
		bornRadius[atom1] = 1.0/(1.0/(atomicRadius[atom1]-offset)-1.0/atomicRadius[atom1]*tanh(delta*psi[atom1]-beta*psi[atom1]*psi[atom1]+gamma*psi[atom1]*psi[atom1]*psi[atom1]));

//		printf("Atom: %5d Radius: %20.10f Psi: %20.10f\n",atom1,bornRadius[atom1],screening);

	}

	/* GB Phase 2: Compute dEijdrij, dai/drij and GB energy */
	gbE = 0;
	for (atom1=0;atom1<nAtoms;atom1++) {
		for(i=1;i<=neighborList[atom1][0];i++) {
			atom2 = neighborList[atom1][i];
			/* avoid double counting */
			if (atom2>atom1) {
				/* compute the distance between the atoms */
				dist1_2=0;
				dx = pos[atom1][0]-pos[atom2][0];
				dy = pos[atom1][1]-pos[atom2][1];
				dz = pos[atom1][2]-pos[atom2][2];
				dist1_2 = dx*dx+dy*dy+dz*dz;
				dist1 = sqrt(dist1_2);
				qiqj = ke*charge[atom1]*charge[atom2];
				aiaj = bornRadius[atom1]*bornRadius[atom2];

				/* compute GB energy term for this pair */

				fij = sqrt(dist1_2 + aiaj*exp(-dist1_2*0.25/aiaj));

//				hij = 1/dist1 - 1/fij;

				Dij = 1.0/epsSolute - exp(-kappa*fij)/epsSolvent;

				gbEij = -qiqj*Dij/fij;
				/* NAMD smoothing function */
				scale = dist1_2/rCut2-1.0;
				scale*=scale;
				/* accumulate energy */
				gbE += scale*gbEij;

				/* accumulate forces for dEijdrij terms */
				dfijdrij = dist1/fij*(1-0.25*exp(-dist1_2*0.25/aiaj));
				dDijdrij = kappa/epsSolvent*exp(-kappa*fij)*dfijdrij;
				dEijdrij = -qiqj/fij*(dDijdrij-Dij/fij*dfijdrij);
				dscaledrij = 4.0*dist1*(dist1_2-rCut2)/(rCut2*rCut2);
				fdEijdrij = -dEijdrij*scale-gbEij*dscaledrij;
				printf("dEijdrij:%20.10f between atom %d and %d. Eij: %20.10f dEijdrij: %20.10f scaled %20.10f\n",fdEijdrij,atom1,atom2,gbEij,dEijdrij,scale);

				/* add to atomic force array */
				force[atom1][0] +=  fdEijdrij*dx/dist1;
				force[atom1][1] +=  fdEijdrij*dy/dist1;
				force[atom1][2] +=  fdEijdrij*dz/dist1;
				force[atom2][0] += -fdEijdrij*dx/dist1;
				force[atom2][1] += -fdEijdrij*dy/dist1;
				force[atom2][2] += -fdEijdrij*dz/dist1;

				/* Accumulate component of dEda terms */
				tmp_dEda = 0.5*qiqj/fij/fij*(kappa/epsSolvent*exp(-kappa*fij)-Dij/fij)*(aiaj+0.25*dist1_2)*exp(-0.25*dist1_2/aiaj);
				dEda[atom1] += tmp_dEda/bornRadius[atom1]*scale;
				dEda[atom2] += tmp_dEda/bornRadius[atom2]*scale;

			}

		}
		/* Self Energy term */
		fij = bornRadius[atom1];
		Dij = 1.0/epsSolute - exp(-kappa*fij)/epsSolvent;
		selfE = -ke * Dij * charge[atom1]*charge[atom1]/(2*fij);
		gbE += selfE;
		/* Accumulate daidr terms */
		tanh2 = tanh(psi[atom1]*(delta+psi[atom1]*(-beta+gamma*psi[atom1])));
		tanh2*=tanh2;
		daidr[atom1] = bornRadius[atom1]*bornRadius[atom1]*(atomicRadius[atom1]-offset)/atomicRadius[atom1]*(1-tanh2)*(delta-2*beta*psi[atom1]+3*psi[atom1]*psi[atom1]);

	}

	printf("Generalized Born Energy: %20.10f Coulomb Energy: %20.10f Total: %20.10f\n",gbE,coulE,coulE+gbE);

	/* GB Phase 3: Compute dETGBddai */
	for (atom1=0;atom1<nAtoms;atom1++) {
		for(i=1;i<=neighborList[atom1][0];i++) {
			atom2 = neighborList[atom1][i];
			/* avoid double counting */
			if (atom2>atom1) {
				/* compute the distance between the atoms */
				dist1_2=0;
				dx = pos[atom1][0]-pos[atom2][0];
				dy = pos[atom1][1]-pos[atom2][1];
				dz = pos[atom1][2]-pos[atom2][2];
				dist1_2 = dx*dx+dy*dy+dz*dz;
				dist1 = sqrt(dist1_2);

				/* compute remaining components of daidr terms */
				dhij = compute_dH(dist1,atomicRadius[atom1],atomicRadius[atom2],offset,atomicScaling[atom2],aCut);
				dhji = compute_dH(dist1,atomicRadius[atom2],atomicRadius[atom1],offset,atomicScaling[atom1],aCut);

				forceMag = dEda[atom1]*daidr[atom1]*dhij+dEda[atom2]*daidr[atom2]*dhji;
				fx = dx/dist1*forceMag;
				fy = dy/dist1*forceMag;
				fz = dz/dist1*forceMag;

				/* add to atomic force array */
				force[atom1][0] +=  fx;
				force[atom1][1] +=  fy;
				force[atom1][2] +=  fz;
				force[atom2][0] += -fx;
				force[atom2][1] += -fy;
				force[atom2][2] += -fz;
			}
		}
	}


}


int exclude_check(int atom1,int atom2,int **excludeList) {

	int i;
	int flag=0;

	for (i=1;i<=excludeList[atom1][0];i++) {
		if (atom2==excludeList[atom1][i]) {
			flag = 1;
			break;
		}
	}
	return flag;

}

double compute_H(double rij, double rhoi, double rhoj, double offset, double Sj, double rCut) {

	double scaledRhoj = (rhoj-offset)*Sj;
	double sRhoj2 = scaledRhoj*scaledRhoj;
	double rhoi0 = rhoi-offset;
	double a=0.333333333333333;
	double b=0.4;
	double c=0.428571428571428;
	double d=0.444444444444444;
	double e=0.454545454545454;
	double rij2 = rij*rij;
	double rc2 = rCut*rCut;
	double k;

	if (rij > (scaledRhoj + rCut)) { //h0
		return 0.0;
	} else if (rij <= (scaledRhoj+rCut) && rij > (rCut - scaledRhoj)) { //h1
		return 1.0/(8.0*rij)*(1.0+2.0*rij/(rij-scaledRhoj)+1.0/(rc2)*(rij2-4.0*rCut*rij-sRhoj2)+2.0*log((rij-scaledRhoj)/rCut));
	} else if (rij <= (rCut-scaledRhoj) && rij > (4*scaledRhoj)) { //h2
		k = sRhoj2/rij2;
		return k*scaledRhoj/rij2 * (a+k* (b+k* (c+k* (d+k*e))));
	} else if (rij <= 4*scaledRhoj && rij > (rhoi0+scaledRhoj)) { //h3
		return 0.5*(scaledRhoj/(rij2-sRhoj2)+1.0/(2.0*rij)*log((rij-scaledRhoj)/(rij+scaledRhoj)));
	} else if ( rij <= (rhoi0+scaledRhoj) && rij > fabs(rhoi0-scaledRhoj)) { //h4
		return 0.25*(1.0/rhoi0*(2.0-1.0/(2.0*rij*rhoi0)*(rij2+rhoi0*rhoi0-sRhoj2))-1.0/(rij+scaledRhoj)+1.0/rij*log(rhoi0/(rij+scaledRhoj)));
	} else if (rhoi0 < scaledRhoj) { //h5
		return 0.5*(scaledRhoj/(rij2-sRhoj2)+2.0/rhoi0+1.0/(2.0*rij)*log((scaledRhoj-rij)/(rij+scaledRhoj)));
	} else { //h6
		return 0.0;
	}



}

double compute_dH(double rij, double rhoi, double rhoj, double offset, double Sj, double rCut) {

	double sRhoj = (rhoj-offset)*Sj;
	double sRhoj2 = sRhoj*sRhoj;
	double rhoi0 = rhoi-offset;
	double rhoi02 = rhoi0*rhoi0;
	double da = 1.333333333333333; // 4* 1/3
	double db = 2.4; // 6* 2/5
	double dc = 3.428571428571428; // 8* 3/7
	double dd = 4.444444444444444; // 10*4/9
	double de = 5.454545454545454; // 12*5/11
	double rij2 = rij*rij;
	double rc2 = rCut*rCut;
	double k;
	double k2;

	if (rij > (sRhoj + rCut)) { //h0
		return 0.0;
	} else if (rij <= (sRhoj+rCut) && rij > (rCut - sRhoj)) { //h1
		return -(rCut+sRhoj-rij)*(sRhoj2-rij2)/(8.0f*rc2*rij2*(sRhoj-rij)*(sRhoj-rij))-0.25/rij2*log((rij-sRhoj)/rCut);
	} else if (rij <= (rCut-sRhoj) && rij > (4*sRhoj)) { //h2
		k = sRhoj2/rij2; // rhojs^2/rij^2
		k2 = k*sRhoj/(rij2*rij); // rhojs^3/rij^5
		return -k2*(da+k*(db+k*(dc+k*(dd+k*de))));
	} else if (rij <= 4*sRhoj && rij > (rhoi0+sRhoj)) { //h3
		return -0.5*( sRhoj*(rij2+sRhoj2)/(rij*(rij2-sRhoj2)*(rij2-sRhoj2)) + 0.5/rij2*log((rij-sRhoj)/(rij+sRhoj)) );
	} else if ( rij <= (rhoi0+sRhoj) && rij > fabs(rhoi0-sRhoj)) { //h4
		return -0.25*( 0.5/rhoi02 - (rij2*(rhoi02-sRhoj2)-2.0*rij*sRhoj2*sRhoj+sRhoj2*(rhoi02-sRhoj2))/(2.0*rij2*rhoi02*(rij+sRhoj)*(rij+sRhoj)) + 1.0/rij2*log(rhoi0/(rij+sRhoj)) );
	} else if (rhoi0 < sRhoj) { //h5
		return -0.5*( sRhoj*(rij2+sRhoj2)/(rij*(rij-sRhoj2)*(rij-sRhoj2)) + 0.5/rij2*log((sRhoj-rij)/(rij+sRhoj)) );
	} else { //h6
		return 0.0;
	}



}

