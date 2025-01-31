#include<stdio.h>
#include<complex.h>
#include<stdlib.h>

#define MAX_LEN 50
#define MAX_ITER 100000

void printpoly(double complex coeffs[MAX_LEN], int deg){
    if (deg > MAX_LEN)
    {
        printf("change the value of MAX_LEN\n");
    }
    else
    {   
        printf("(%.2f + %.2f i)", creal(coeffs[0]), cimag(coeffs[0]));
        for (int i = 1; i <= deg; i++)
        {
            printf(" +(%.2f + %.2f i)*x^%d", creal(coeffs[i]), cimag(coeffs[i]), i);
        }
        printf("\n");
    }
}

// divide polynomial c[0] + c[1] * x + c[2] * x ^ 2 + ... + c[deg] * x^deg by the polynomial x - a. The new coefficients are stored in coeffs
// and the remainder in the rem variable.

void poly_div(double complex coeffs[MAX_LEN], int deg, double complex a, double complex *rem){
    double complex swap;
    *rem = coeffs[deg];
    coeffs[deg] = 0.0;
    for(int i = deg - 1; i >= 0; i--) 
    {
        swap = coeffs[i];
        coeffs[i] = *rem;
        *rem = swap + (*rem) * a;
    }
}

void first_and_second_derivative(double complex coeffs[MAX_LEN], double complex first_der[MAX_LEN], double complex second_der[MAX_LEN], int deg){
    first_der[0] = coeffs[1];
    for (int i = 2; i <= deg; i++)
    {
        first_der[i - 1] = ((double complex) i) * coeffs[i];
        second_der[i - 2] = ((double complex) (i * (i - 1))) * coeffs[i]; 
    }
}

double complex evalpol(double complex coeffs[MAX_LEN], double complex x, int deg, int *error){
    double complex res = 0.0;
    int j;
    if (deg >= MAX_LEN){
        fprintf(stderr, "increase MAX_LEN");
        *error = 1;
    }
    else
    {
    res = coeffs[j = deg];
    while (j > 0) res = res * x + coeffs[--j];
    *error = 0;
    }   
    return res;
}

int main(void){
    double complex coeffs[] = {4.0, 2.0, 1.0}, res;
    double complex first_der[MAX_LEN], second_der[MAX_LEN];
    int deg = 2, error, iter = 0, initial_deg;
    double complex root_list[deg], G, H, valp, valm, valsqrt, denom, a, rem;
    double accuracy_target = 1e-10, accuracy;
    FILE *fp;
    fp = fopen("./../data/roots.dat", "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", "./../data/roots.dat", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    fprintf(fp, "#root_real\troot_imag\titerations\taccuracy\n");

    for (int j = 0; j < deg; j++)
    {
        root_list[j] = 1.0;
    }
    
    initial_deg = deg;
    for (int j = 0; j < initial_deg; j++)
    {   
        first_and_second_derivative(coeffs, first_der, second_der, deg);
        res = evalpol(coeffs, root_list[j], deg, &error);
        accuracy = cabs(res);
        while (accuracy > accuracy_target && iter < MAX_ITER)
        {
            G = evalpol(first_der, root_list[j], deg - 1, &error) / evalpol(coeffs, root_list[j], deg, &error);
            H = G * G - evalpol(second_der, root_list[j], deg - 2, &error) / evalpol(coeffs, root_list[j], deg, &error);
            valsqrt = csqrt((deg - 1) * (deg * H - G * G));
            valp = G + valsqrt;
            valm = G - valsqrt;
            denom = (cabs(valp) > cabs(valm)) ? valp : valm;
            a = deg / denom;
            root_list[j] -= a;
            res = evalpol(coeffs, root_list[j], deg, &error);
            accuracy = cabs(res);
            iter++;
        }
        poly_div(coeffs, deg, root_list[j], &rem);
        deg--;
        fprintf(fp, "%.10e\t%.10e\t%d\t%.10e\n", creal(root_list[j]), cimag(root_list[j]), iter, accuracy);
    }
    fclose(fp);    
    return EXIT_SUCCESS;
}

