#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

void diagonalize_matrix(int n, double** A, double* eigenvalues, double* eigenvectors) {
    // LAPACK variables
    int lda = n;
    int ldvr = n;
    int lvd1=n;
    int lwork = 4 * n; // Size of work array
    double* wr = (double*)malloc(n * sizeof(double)); // Real parts of eigenvalues
    double* wi = (double*)malloc(n * sizeof(double)); // Imaginary parts of eigenvalues
    double* vl = NULL; // Left eigenvectors (not used here)
    double* vr = (double*)malloc(n * n * sizeof(double)); // Right eigenvectors
    double* work = (double*)malloc(lwork * sizeof(double)); // Work array
    int info;

    // Convert 2D array to 1D array required by LAPACK
    double* A_flat = (double*)malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_flat[i * n + j] = A[i][j];
        }
    }

    // Call LAPACK's dgeev to compute eigenvalues and eigenvectors
    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, A_flat, lda, wr, wi, vl, lvd1, vr, ldvr);

    if (info != 0) {
        fprintf(stderr, "dgeev failed with error code %d\n", info);
        exit(EXIT_FAILURE);
    }
   // Copy eigenvalues to the output array
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = wr[i];
        if (wi[i] != 0.0) {
            printf("Complex eigenvalue detected: %.4f + %.4fi\n", wr[i], wi[i]);
        }
    }

    // Copy eigenvectors to the output array
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            eigenvectors[i * n + j] = vr[i * n + j];
        }
    }
       // Free allocated memory
    free(wr);
    free(wi);
    free(vr);
    free(work);
    free(A_flat);
}


double** malloc_2d(int m, int n){
 double** a = malloc(m*sizeof(double*));
 if (a == NULL){
  return NULL;
 }
 a[0] = malloc(n*m*sizeof(double));
 if (a[0] == NULL) {
  free(a);
  return NULL;
 }
 for (int i=1;i<m;i++){
        a[i] = a[i-1]+n;
 }
 return a;
}

void free_2d(double** a){
free(a[0]);
a[0]=NULL;
free(a);
}


