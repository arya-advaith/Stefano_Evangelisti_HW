#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_solve.h"

int main() {
    FILE *fp = fopen("input.txt", "r");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }
    
    int n_c;
    double *alpha = (double*)calloc(2, sizeof(double));
    double *beta = (double*)calloc(2, sizeof(double));
    int open;
    int read_count;

    if (!alpha || !beta) {
        fprintf(stderr, "Memory allocation failed\n");
        free(alpha);
        free(beta);
        fclose(fp);
        return 1;
    }

    // Try reading with one or two values for each
    read_count = fscanf(fp, "%*s %d \n %*s %lf %lf \n %*s %lf %lf \n %*s %d",
            &n_c, &alpha[0], &alpha[1], &beta[0], &beta[1], &open);
    
    if (read_count != 6) {
        rewind(fp);
        read_count = fscanf(fp, "%*s %d \n %*s %lf \n %*s %lf %lf \n %*s %d",
                &n_c, &alpha[0], &beta[0], &beta[1], &open);
        
	if (read_count != 5) {
		rewind(fp);
		read_count = fscanf(fp, "%*s %d \n %*s %lf \n %*s %lf \n %*s %d",
                &n_c, &alpha[0], &beta[0], &open);
			
			if (read_count != 4){
            			printf("Error in input file\n");
            			free(alpha);
            			free(beta);
            			fclose(fp);
            			return 1;
        		}
    	}
    }

    printf("Number of atoms: %d\n", n_c);
    printf("Alpha values: %lf", alpha[0]);
    if (read_count == 6) {
	    printf(" %lf", alpha[1]);
    }
    printf("\nBeta values: %lf", beta[0]);
    if (read_count >= 5) {
	    printf(" %lf", beta[1]);
    }
    printf("\n");
    printf("Type of System: \n");
    if(open==1){
	printf("Open or Linear\n");
     }
     else{
	printf("Cyclic\n");
     }

    fclose(fp);

// Allocation of the Huckel Matrix 
	double** huckel=malloc_2d(n_c,n_c);

// Now Allocating the Beta values:
   	int numb;
	if (read_count >= 5){
		numb = 2;
	} 
	else {
		numb = 1;
	}
	printf("num of beta %d \n",numb);

// Now Allocating the lpha values
	int numa;
	if (read_count == 6){
		numa = 2;
	}
// If the system is cyclic
	if (open ==0){
			huckel[n_c-1][0]=beta[0];
			huckel[0][n_c-1]=beta[0];
			if (numb == 2){
				if (n_c%numb ==0){
					huckel[n_c-1][0]=beta[1];
					huckel[0][n_c-1]=beta[1];
				}
			}
	}

// Allocating the general Beta values
	for (int i=0;i<n_c-1;i++){
			
			huckel[i][i+1]=beta[0];
			huckel[i+1][i]=beta[0];
		
		if (numb == 2){
			
			if (i%numb == 1 ){
			huckel[i][i+1]=beta[1];
                	huckel[i+1][i]=beta[1];
			}
		}
	}

//Allocating the general Alpha values
	for (int i=0;i<n_c;i++){

               huckel[i][i]=alpha[0];
               if (numa == 2){

                        if (i%numa == 1 ){
                        huckel[i][i]=alpha[1];
                        }
                }
        }	


	printf("Successful building of Huckel Matrix: \n");

// Printing the Huckel Matrix  
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

// Writing the eigenvalues to a file
	FILE* fptr;
     	fptr = fopen("eigenvalues_.txt","w");
    	for (int i = 0; i < n_c; i++) {
       		fprintf(fptr,"%.4f", eigenvalues[i]);
       	fputs("\n",fptr);
    	}
	fclose(fptr);

// Writing the eigenvectors to a file
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

// Freeing the allocated space
    free_2d(huckel);
    free(eigenvalues);
    free(eigenvectors);
    free(alpha);
    free(beta);
    
    return 0;
}
