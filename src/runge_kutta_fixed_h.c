#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#

#define DIM 5
#define UNUSED(t) (void)(t)
#define N_params 13
#define MAX_ITER (int) pow(10, 7)
#define MAX_LEN_FILE 100
#define N_simul 1000

double c[6] = {0.0, 0.2, 0.3, 0.6, 1.0, 0.875};
double a[6][5] = {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.2, 0.0, 0.0, 0.0, 0.0},
                  {0.075, 0.225, 0.0, 0.0, 0.0}, {0.3, -0.9, 1.2, 0.0, 0.0},
                  {-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0},
                  {1631.0/55296.0,175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0}};
double b[6] = {37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0};
double b_ast[6] = {2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25};

double a1 = 1.0, a2 = 25.0;
double b1 = 0.6, b2 = 2.3;
double c1 = 0.646, c2 = 0.6, c34 = 0.2;
double m1 = 2.3;
double Q = 5.2e-2, q1 = 0.8, w1 = 2.4, w2 = 0.02, q54 = 0.6;

double yscale(double y[DIM], double field[DIM], double h){
    double ymin, fieldmin;
    ymin = fabs(y[0]);
    fieldmin = fabs(field[0]);
    for (int jdim = 1; jdim < DIM; jdim++)
    {
        if (fabs(y[jdim]) < ymin)
        {
            ymin = fabs(y[jdim]);
        }
        if (fabs(field[jdim]) < fieldmin)
        {
            fieldmin = fabs(field[jdim]);
        }
        
    }
    return ymin + h * fieldmin;
}

void func2d(double x[DIM], double field[DIM], double t, double params[N_params]){
    UNUSED(t);
    field[0] = params[0] * x[0] + params[1] * x[1];
    field[1] = - params[1] * x[0] + params[0] * x[1];
}

void cancer_system(double x[DIM], double field[DIM], double t, double params[N_params]){
    UNUSED(t);
    field[0] = params[0] * x[0] - params[1] * x[0] * x[0];
    field[1] = params[2] * x[1] - params[3] * x[1] * x[1] + x[0] - x[1] * x[3];
    field[2] = -params[4] * x[2] + params[5] * x[1] - params[6] * x[2] * x[3];
    field[3] = x[3] * (params[7] - params[7] * x[3] - x[1] - x[2] - x[4]);
    field[4] = params[8] - params[9] * x[4] + params[10] * x[3] / (params[11] + x[3]) * x[4] - params[12] * x[4] * x[3];
}



void rungekutta(void (*system)(double xin[DIM], double xout[DIM], double t, double params[N_params]), double params[N_params],
                double yn[DIM], double tn, double tfin, double h, double rel_acc, char filename[20]){
    UNUSED(rel_acc);
    int iter;
    double err[DIM], temp[DIM], field[DIM], maxerr=0.0, yscale_loc;
    double k[6][DIM];
    FILE *fp;

    fp = fopen(filename, "w");
    //N = (int) ((tfin - tn) / h);
    if (fp == NULL)
    {
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    fprintf(fp, "t\tx1\tx2\tx3\tx4\tx5\terr\th\tyscale\t"
                "err1\terr2\terr3\terr4\terr5\n"
                "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.10e\t%.10e\t%.10e\t"
                "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", 
                tn, yn[0], yn[1], yn[2], yn[3], yn[4], maxerr, h, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    //fprintf(fp, "#t\tx1\tx2\tabserr\th\tyscale\n"
    //            "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", 
    //            tn, yn[0], yn[1], maxerr, h, 1.0);

    iter = 0;
    while (tn < tfin && iter < MAX_ITER)    
    {   
        for (int klen = 0; klen < 6; klen++)
        {
            for (int jdim = 0; jdim < DIM; jdim++)
            {
                temp[jdim] = 0.0;
                
                for (int isum = 0; isum < klen; isum++)
                {
                    temp[jdim] += a[klen][isum] * k[isum][jdim];
                }
                temp[jdim] = yn[jdim] + h * temp[jdim];           
            }

            system(temp, field, tn + c[klen] * h, params);

            for (int jdim = 0; jdim < DIM; jdim++)
            {
                k[klen][jdim] = field[jdim];
            }
            
        }

        for (int jdim = 0; jdim < DIM; jdim++)
        {
            temp[jdim] = yn[jdim];
        }
        
        for (int j = 0; j < 6; j++)
        {
            for (int jdim = 0; jdim < DIM; jdim++)
            {
                temp[jdim] += h * b[j] * k[j][jdim];
            }
        }
        
        for (int jdim = 0; jdim < DIM; jdim++)
        {
            err[jdim] = 0.0;
            for (int j = 0; j < 6; j++)
            {
                err[jdim] += h * (b[j] - b_ast[j]) * k[j][jdim];
            } 
        }
        

        maxerr = fabs(err[0]);
        for (int jdim = 1; jdim < DIM; jdim++)
        {
            if (fabs(err[jdim]) > maxerr)
            {
                maxerr = fabs(err[jdim]);
            }
            
        }

        system(yn, field, tn, params);
        yscale_loc = yscale(yn, field, h);
        tn += h;
        for (int jdim = 0; jdim < DIM; jdim++)
        {
            yn[jdim] = temp[jdim];
        }
        fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", 
                    tn, yn[0], yn[1], yn[2], yn[3], yn[4], maxerr, h, yscale_loc, err[0], err[1], err[2], err[3], err[4]);
        iter++;
    }
    fclose(fp);
}

int main(void){
    double yin[DIM], tin, tfin, h, rel_acc;
    double params[N_params + 1];
    char filename_data[200], filename_params[200];
    int count = 0, index;
    FILE *fp;

    tin = 0.0;
    tfin = 10.0;
    h = 1.0e-6;
    rel_acc = 1e-12;
   
    strcpy(filename_params, "./../data/parameters.dat");
    fp = fopen(filename_params, "r");

    if (fp == NULL)
    {
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename_params, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    //params[0] = a1; params[1] = a2;
    //params[2] = b1; params[3] = b2;
    //params[4] = c1; params[5] = c2; params[6] = c34;
    //params[7] = m1;
    //params[8] = Q; params[9] = q1; params[10] = w1; params[11] = w2; params[12] = q54;

    //fprintf(fp, "#a1\ta2\tb1\tb2\tc1\tc2\tc34\tm1\tQ\tq1\tw1\tw2\tq54"
    //            "%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e",
    //            a1, a2, b1, b2, c1, c2, c34, m1, Q, q1, w1, w2, q54);
    while (fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                &index,
                &params[0], &params[1], &params[2], &params[3], &params[4], &params[5], &params[6], 
                &params[7], &params[8], &params[9], &params[10], &params[11], &params[12]) !=EOF && ++count < N_simul)
    {
        if (index == 12)
        {
                   
        printf("%d\n", count);
        /* let us set, as initial conditions, the asymptotic tumor-free behavior, and a small initial value for tumor cells */
        yin[0] = params[0] / params[1];
        yin[1] = 0.5 * (params[2] / params[3] + sqrt((params[2] / params[3] * params[2] / params[3]) + 4 * params[0] / params[1] / params[3]));
        yin[2] = params[5] / params[4] * yin[1];
        yin[3] = 1.0e-6;
        yin[4] = params[8] / params[9];

        snprintf(filename_data, sizeof(filename_data), "./../data/cancer_%d_test.dat", index);
        rungekutta(&cancer_system, params, yin, tin, tfin, h, rel_acc, filename_data);       
        }        
    }

    fclose(fp);

    
    
    //params[0] = -1.0;
    //for (int i = 0; i < 10; i++)
    //{   
    //    printf("%d\n", i);
    //    for (int jdim = 0; jdim < DIM; jdim++)
    //    {
    //        yin[jdim] = 1.0;
    //    }
    //    params[1] = 10.0 / (double) (i + 1);
    //    snprintf(filename, sizeof(filename), "data_%d.dat", i);
    //    tin = 0.0;
    //    tfin = 100.0;
    //    h = 1.0e-3;
    //    rel_acc = 1.0e-10;
    //    rungekutta(&func2d, params, yin, tin, tfin, h, rel_acc, filename);    
    //}
    return EXIT_SUCCESS;
}

