#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

void read_psf_header(FILE *, int *);

char** read_psf_atom_data(FILE  *, int, double *, int *, int *, double *, double *);

void read_psf_bond_data(FILE *, int **);

