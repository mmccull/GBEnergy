
/*
 *  A set of subroutines used for read, writing and manipulating coordinates from a DCD file
 * 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void copy_dcd_header(FILE *, FILE *, int *, int *);

void read_dcd_header(FILE *, int *, int *);

void read_dcd_step(FILE *, long int *, double *, double **, int, int );

void read_dcd_COM_coord(FILE *, long int *, double **, double **, double **, int, int, double *, int *, int);

void read_dcd_step_selectatoms(FILE *, long int *, double *, double **, int, int *, int);

void center_coord(double **, double *, int );

void wrap_coord(double *, double **, int , int **, int , double *);

void write_dcd_coord(FILE *, double *, double **,int );
