#include "psflib.h"

void read_psf_header( FILE *psfFile, int *nAtoms)
{
	char buffer[1024];
	char key[8];
	char temp[10];
	
	
	
	/*
		Read the number of atoms from the header of the psf file
	*/
	while( fgets( buffer, 1024, psfFile ) != NULL )
	{
		strncpy(key,buffer+9,6);
//		key[6]='\0';
		if( strncmp(key,"!NATOM",6)==0 )
		{
			strncpy(temp,buffer,9);
			temp[9]='\0';
			*nAtoms = atoi( temp );
			break;
		}
	}
}

char** read_psf_atom_data(FILE  *psfFile, int nAtoms, double *charge, int *nAtomTypes, int *uniqueAtomTypeNum, double *atomicRadius, double *atomicScaling) {

	char buffer[1024];
	int i, j;
	char temp[6];
	int flag;
	int type;
	char atomType[6];
	char atomName[6];
	char chargeTemp[13];
	char massTemp[9];
	char uniqueAtomType[nAtoms][5];
	char **smallUniqueType;
	int uniqueTypeCount;
	double mass;


	uniqueTypeCount = 0;
	for (i=0;i<nAtoms;i++) {
		fgets(buffer,1024,psfFile);
		strncpy(temp,buffer+24,4);
		temp[4]='\0';
		strncpy(atomType,buffer+29,4);
		atomType[4]='\0';
		strncpy(chargeTemp,buffer+33,11);
		chargeTemp[11]='\0';
		charge[i] = atof(chargeTemp);
		strncpy(massTemp,buffer+50,8);
		massTemp[8]='\0';
		mass = atof(massTemp);
		if (mass < 2.50) { // H
			atomicRadius[i] = 1.20;
			atomicScaling[i] = 0.85;
		} else if (mass < 12.50) { // C
			atomicRadius[i] = 1.70;
			atomicScaling[i] = 0.72;
		} else if (mass < 14.50) { // N
			atomicRadius[i] = 1.55;
			atomicScaling[i] = 0.79;
		} else if (mass < 16.50) { // 0
			atomicRadius[i] = 1.50;
			atomicScaling[i] = 0.85;
		} else if (mass < 19.50) { // F
			atomicRadius[i] = 1.50;
			atomicScaling[i] = 0.88;
		} else if (mass < 23.64) { // Na
			atomicRadius[i] = 2.27;
			atomicScaling[i] = 0.80;
		} else if (mass < 31.50) { // P
			atomicRadius[i] = 1.85;
			atomicScaling[i] = 0.86;
		} else if (mass < 32.50) { // S
			atomicRadius[i] = 1.80;
			atomicScaling[i] = 0.96;
		} else if (mass < 37.28) { // Cl
			atomicRadius[i] = 1.70;
			atomicScaling[i] = 0.80;
		} else { // all others
			atomicRadius[i] = 1.50;
			atomicScaling[i] = 0.80;
		}
		/* convert mass to radius and radius scaling */
		if (uniqueTypeCount==0) {
			strncpy(uniqueAtomType[uniqueTypeCount],atomType,5);
			uniqueAtomTypeNum[i] = uniqueTypeCount;
			uniqueTypeCount++;
		} else {
			flag = 0;
			for (type=0;type<uniqueTypeCount;type++) {
				if (strncmp(uniqueAtomType[type],atomType,4)==0) {
					flag = 1;
					uniqueAtomTypeNum[i] = type;
					break;
				} 
			}
			if (flag == 0) {
				strncpy(uniqueAtomType[uniqueTypeCount],atomType,5);
				uniqueAtomTypeNum[i] = uniqueTypeCount;
				uniqueTypeCount++;
			}
		}

	}

	*nAtomTypes = uniqueTypeCount;

	smallUniqueType = new char* [*nAtomTypes];

	for (type=0;type<*nAtomTypes;type++) {
		smallUniqueType[type] = new char [5];
		strncpy(smallUniqueType[type],uniqueAtomType[type],5);
	}


	return smallUniqueType;

	delete [] smallUniqueType;

}



void read_psf_bond_data( FILE *psfFile, int **excludeList) 
{
	int bondsPerLine = 4;
	int anglesPerLine = 3;
	int dihsPerLine = 2;
	int dih;
	int dihLineCount;
	int nDihs;
	int angle;
	int angleLineCount;
	int nAngles;
	char buffer[1024];
	char key[8];
	char temp[10];
	int nBonds;
	int nLines;
	int line;
	int bond;
	int bondLineCount;
	int atom1;
	int atom2;
	int atom3;
	int atom4;
	int stringPos;
	int i, j;
	
	

	/*
		Read psf line by line
	*/
	while( fgets( buffer, 1024, psfFile ) != NULL )
	{
		strncpy(key,buffer+9,6);
		if( strncmp(key,"!NBOND",6)==0 )
		{
			strncpy(temp,buffer,9);
			temp[9]='\0';
			nBonds = atoi( temp );
			printf("number of bonds: %d\n",nBonds);
			break;
		}
	}
	nLines = (int) ( ((double) nBonds) / ((double) bondsPerLine) +0.75);

	bond = 0;
	for (line=0;line<nLines;line++) {
		
		fgets( buffer, 1024, psfFile );
		stringPos = 0;
		for(bondLineCount=0;bondLineCount<bondsPerLine;bondLineCount++) {
			if (bond>=nBonds) break;
			bond++;
			strncpy(temp,buffer+stringPos,8);
			temp[8] = '\0';
			atom1 = atoi(temp)-1;
			strncpy(temp,buffer+stringPos+8,8);
			atom2 = atoi(temp)-1;
			excludeList[atom1][excludeList[atom1][0]+1] = atom2;
			excludeList[atom1][0]++;
			excludeList[atom2][excludeList[atom2][0]+1] = atom1;
			excludeList[atom2][0]++;
			stringPos+=16;
		}	

	}


	/* Now read angles */

	while( fgets( buffer, 1024, psfFile ) != NULL )
	{
		strncpy(key,buffer+9,7);
		if( strncmp(key,"!NTHETA",7)==0 )
		{
			strncpy(temp,buffer,9);
			temp[9]='\0';
			nAngles = atoi( temp );
			printf("number of angles: %d\n",nAngles);
			break;
		}
	}
	nLines = (int) ( ((double) nAngles) / ((double) anglesPerLine) +0.66);

	angle = 0;
	for (line=0;line<nLines;line++) {
		
		fgets( buffer, 1024, psfFile );
		stringPos = 0;
		for(angleLineCount=0;angleLineCount<anglesPerLine;angleLineCount++) {
			if (angle>=nAngles) break;
			angle++;
			strncpy(temp,buffer+stringPos,8);
			temp[8] = '\0';
			atom1 = atoi(temp)-1;
			strncpy(temp,buffer+stringPos+8,8);
			atom2 = atoi(temp)-1;
			strncpy(temp,buffer+stringPos+16,8);
			atom3 = atoi(temp)-1;
			excludeList[atom1][excludeList[atom1][0]+1] = atom3;
			excludeList[atom1][0]++;
			excludeList[atom3][excludeList[atom3][0]+1] = atom1;
			excludeList[atom3][0]++;
			stringPos+=24;
		}	

	}

	/* Now read dihedrals (only if we want to exclude 1-4 nonbonded interactions */
/*
	while( fgets( buffer, 1024, psfFile ) != NULL )
	{
		strncpy(key,buffer+9,5);
		if( strncmp(key,"!NPHI",5)==0 )
		{
			strncpy(temp,buffer,9);
			temp[9]='\0';
			nDihs = atoi( temp );
			printf("number of dihedrals: %d\n",nDihs);
			break;
		}
	}
	nLines = (int) ( ((double) nDihs) / ((double) dihsPerLine) +0.50);

	dih = 0;
	for (line=0;line<nLines;line++) {
		
		fgets( buffer, 1024, psfFile );
		stringPos = 0;
		for(dihLineCount=0;dihLineCount<dihsPerLine;dihLineCount++) {
			if (dih>=nDihs) break;
			dih++;
			strncpy(temp,buffer+stringPos,8);
			temp[8] = '\0';
			atom1 = atoi(temp)-1;
			strncpy(temp,buffer+stringPos+8,8);
			atom2 = atoi(temp)-1;
			strncpy(temp,buffer+stringPos+16,8);
			atom3 = atoi(temp)-1;
			strncpy(temp,buffer+stringPos+24,8);
			atom4 = atoi(temp)-1;
			excludeList[atom1][excludeList[atom1][0]+1] = atom4;
			excludeList[atom1][0]++;
			excludeList[atom4][excludeList[atom4][0]+1] = atom1;
			excludeList[atom4][0]++;
			stringPos+=32;
		}	
	}
*/	
}

