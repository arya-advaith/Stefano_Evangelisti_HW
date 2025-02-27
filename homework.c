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

// Total Number of Electrons 
	printf("Enter the total number of atoms: \n");
	int  n_c;
	scanf("%d",&n_c);

// Allocation of the Huckel Matrix 
	double** huckel=malloc_2d(n_c,n_c);

// Asking if There are different atoms in the system
	printf("Different Atoms in the system? [Yes=1]: \n");
	int differ;
	scanf("%d",&differ);

// Adding the Alpha values if there are different atoms 
	if(differ==1){

		// Just a trial! Asking if the atoms alternate in the system or not.
			printf("Do the 2 (or more) atoms alternate? [No= 0] : \n");
			int altern;
			scanf("%d",&altern);
			
			// If it does not alternate, then manually entering the values of Alpha
				if(altern==0){
					for(int i=0;i<n_c;i++){ 
					printf("Enter the alpha value of atom %d out of %d \n: ",i,n_c);
						int alpha;
						scanf("%d",&alpha);
						huckel[i][i] = alpha;
						}
					}
				
				// If the atoms alternate in the system [2 atoms i guess] 
				else{
					printf("-----------------------------------------------------\n");
					printf("This feature only works for 2 atoms unfortunately :/\n");
					printf("-----------------------------------------------------\n");
					
					printf("Enter the alpha value of atom 1\n:");
					double alpha1;
				        scanf("%lf",&alpha1);
					
					printf("Enter the alpha value of atom 2\n:");
					double alpha2;
					scanf("%lf",&alpha2);
					
					for(int i=0;i<n_c;i++){
						
						if(i%2==0){
						huckel[i][i]=alpha1;
						}
						
						else{
						huckel[i][i]=alpha2;
						}
					}	
				}
}


printf("Enter how many beta values: \n");
int numb;
scanf("%d",&numb);

double beta[numb];
for (int i=0;i<numb;i++){
printf("Enter the beta value: \n");
scanf("%lf",&beta[i]);
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
