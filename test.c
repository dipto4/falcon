#include <stdio.h>
#include <stdlib.h>

void check_2d_array(double** a, size_t N) {
    //double* a_x = a[0];
    //double* a_y = a[1];
    //double* a_z = a[2];
    int i = 0;

    for(i = 0; i<10; i++) {
        printf("%f ",a[0][i]);
        printf("%f ",a[1][i]);
        printf("%f ",a[2][i]);
    }
}



void main() {
    double **array;
    array = (double **) malloc(3*sizeof(double*));

    int row = 0;
    for(row = 0; row<3;row++) {
        array[row] = (double *) malloc(10*sizeof(double));
    }

    int col = 0;

    for (row = 0; row<3;row++) {
        for(col = 0; col<10;col++) {
            array[row][col] = row*10.0e0 + col;
        }
    } 
     int i = 0;

    for(i = 0; i<10; i++) {
        printf("%f ",array[0][i]);
        printf("%f ",array[1][i]);
        printf("%f ",array[2][i]);
    }
    printf("\n");
    check_2d_array(array,10);
}
