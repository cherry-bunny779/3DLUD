#ifndef LU_DECOMP_BLOCKED_H
#define LU_DECOMP_BLOCKED_H

// LU decomposition of a square matrix A[size][size] with block size b_size
// Modifies A in-place into its LU factors
void LU_Decomposition(float **&A, int size, int b_size);

// Initializes A[size][size] with random values
void init(float **A, int size);

// Prints the matrix A[size][size] to stdout
void display(float **A, int size);

// Allocates a 2D float matrix of size [size][size]
float** allocate_2d(int size);

// Frees a 2D float matrix of size [size][size]
void free_2d(float **A, int size);

#endif // LU_DECOMP_BLOCKED_H