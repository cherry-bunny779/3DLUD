#include <stdio.h>
#include <stdlib.h>

// Function prototypes
void display(float **A, int size);
void init(float **A, int size);
void diag_func(float **A, int i, int size, int b_size, int &num_div, int &num_mac, int &num_subt);
void row_func(float **A, int i, int j, int size, int b_size, int &num_mac, int &num_subt);
void col_func(float **A, int i, int j, int size, int b_size, int &num_div, int &num_mac, int &num_subt);
void inner_func(float **A, int i, int j, int k, int size, int b_size, int &num_mac, int &num_subt);

// count ops in all helper functions now

void row_func(float **A, int i, int j, int size, int b_size, int &num_mac, int &num_subt) {
    for (int ii = i*b_size; ii < (i*b_size)+(b_size-1); ii++) {
        for(int jj = ii+1; jj < b_size; jj++) {
            for(int kk = j*b_size; kk < (j*b_size)+b_size; kk++) {
                A[jj][kk] -= A[jj][ii] * A[ii][kk];
                num_mac++;
                num_subt++;
            }
        }
    }
}

void col_func(float **A, int i, int j, int size, int b_size, int &num_div, int &num_mac, int &num_subt) {
    for(int ii = i*b_size; ii < (i*b_size)+b_size; ii++) {
        for(int jj = j*b_size; jj < (j*b_size)+b_size; jj++) {
            A[jj][ii] /= A[ii][ii];
            num_div++;

            for(int kk = ii + 1; kk < (i*b_size)+b_size; kk++) {
                A[jj][kk] -= A[jj][ii] * A[ii][kk];
                num_mac++;
                num_subt++;
            }
        }
    }
}

void inner_func(float **A, int i, int j, int k, int size, int b_size, int &num_mac, int &num_subt) {
    for(int ii = i*b_size; ii < (i*b_size)+b_size; ii++) {
        for(int jj = j*b_size; jj < (j*b_size)+b_size; jj++) {
            for(int kk = k*b_size; kk < (k*b_size)+b_size; kk++) {
                A[jj][kk] -= A[jj][ii] * A[ii][kk];
                num_mac++;
                num_subt++;
            }
        }
    }
}

void diag_func(float **A, int i, int size, int b_size, int &num_div, int &num_mac, int &num_subt) {
    for(int ii = i*b_size; ii < (i*b_size)+b_size-1; ii++) {
        for(int jj = ii+1; jj < (i*b_size)+b_size; jj++) {
            A[jj][ii] /= A[ii][ii];
            num_div++;

            for(int kk = ii+1; kk < (i*b_size) + b_size; kk++) {
                A[jj][kk] -= A[jj][ii] * A[ii][kk];
                num_mac++;
                num_subt++;
            }
        }
    }
}

void init(float **A, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            A[i][j] = rand() % 200 + 2;
        }
    }
}

void display(float **A, int size) {
    printf("---------------------\n");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%8.2f ", A[i][j]);
        }
        printf("\n");
    }
}

void LU_Decomposition(float **&A, int size, int b_size,
                      int &num_div, int &num_mac, int &num_subt) {
    int num_blocks = size / b_size;

    // Initialize counters
    num_div = 0;
    num_mac = 0;
    num_subt = 0;

    for(int i=0; i<num_blocks; i++) {
        diag_func(A, i, size, b_size, num_div, num_mac, num_subt);

        for(int j=i+1; j<num_blocks; j++) {
            row_func(A, i, j, size, b_size, num_mac, num_subt);
        }

        for(int j=i+1; j<num_blocks; j++) {
            col_func(A, i, j, size, b_size, num_div, num_mac, num_subt);

            for(int k=i+1; k<num_blocks; k++) {
                inner_func(A, i, j, k, size, b_size, num_mac, num_subt);
            }
        }
    }
}

// helper to allocate a 2D float array
float** allocate_2d(int size) {
    float **A = new float*[size];
    for (int i = 0; i < size; ++i) {
        A[i] = new float[size];
    }
    return A;
}

// helper to free a 2D float array
void free_2d(float **A, int size) {
    for (int i = 0; i < size; ++i) {
        delete[] A[i];
    }
    delete[] A;
}