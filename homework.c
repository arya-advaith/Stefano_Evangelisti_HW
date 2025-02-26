#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

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

void main(){

printf("Enter the number of Carbon atoms: \n");
int  n_c;
scanf("%d",&n_c);

double** huckel=malloc_2d(n_c,n_c);
printf("Enter how many beta values: \n");
int numb;
scanf("%d",&numb);

float beta[numb];
for (int i=0;i<numb;i++){
printf("Enter the beta value: \n");
scanf("%f",&beta[i]);
}


printf("Enter if Open[=1] or Cyclic[=0]: \n");
int open;
scanf("%d",&open);

if (open ==0){

	if (n_c%numb ==0){
		huckel[n_c-1][0]=beta[1];
		huckel[0][n_c-1]=beta[1];
	}

	else{
		huckel[n_c-1][0]=beta[0];
		huckel[0][n_c-1]=beta[0];
	}
}

for (int i=0;i<n_c-1;i++){

		//printf("%d,%d\t%d,%d\n",i,i+1,i+1,i);
		if (i%numb==0){
			huckel[i][i+1]=beta[0];
			huckel[i+1][i]=beta[0];
		}
		
		else{
			huckel[i][i+1]=beta[1];
                	huckel[i+1][i]=beta[1];
		}

}


printf("Successful building of Huckel Matrix: \n");

for(int i=0;i<n_c;i++){
	for(int j=0;j<n_c;j++){
	printf("%lf\t",huckel[i][j]);
	}
	printf("\n");
}

// Allocate memory for eigenvalues and eigenvectors
    double* eigenvalues = (double*)malloc(n_c * sizeof(double));
    double* eigenvectors = (double*)malloc(n_c * n_c * sizeof(double));

    if (eigenvalues == NULL || eigenvectors == NULL) {
        fprintf(stderr, "Failed to allocate memory for eigenvalues or eigenvectors.\n");
        free_2d(huckel);
    }

    // Diagonalize the matrix
    diagonalize_matrix(n_c, huckel, eigenvalues, eigenvectors);

FILE* fptr;
     fptr = fopen("eigenvalues_.txt","w");
    for (int i = 0; i < n_c; i++) {
       fprintf(fptr,"%.4f", eigenvalues[i]);
       fputs("\n",fptr);
    }
fclose(fptr);

FILE* fptr1;
fptr1=fopen("eigenvectors.txt","w");
    for (int i = 0; i < n_c; i++) {
        for (int j = 0; j < n_c; j++) {
            fprintf(fptr1,"%.4f", eigenvectors[i * n_c + j]);
	    fputs("\t",fptr1);
        }
        fputs("\n",fptr1);
    }
fclose(fptr1);

free_2d(huckel);

}
