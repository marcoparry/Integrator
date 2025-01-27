#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define DIM 2

double func(double x){
    return x * (1 - x);
}

void func2d(double x[DIM], double y[DIM]){
    y[0] = - x[0] + 10 * x[1];
    y[1] = - 10 * x[0] - x[1];
}

double c[] = {0.0, 1/5};

int main(void){
    double x[DIM], t, tfin, dt, y[DIM];
    int i, N;
    char datafile[20];
    FILE *fp;

    strcpy(datafile, "data.dat");
    fp = fopen(datafile, "w");

    x[0] = 1.0;
    x[1] = 1.0;
    t = 0.1;
    tfin = 100.0;
    dt = 1.0e-3;

    N = (int) ((tfin - t) / dt);

    fprintf(fp, "#t\tx\ty\n%.14f\t%.14f\t%.14f\n", t, x[0], x[1]);

    for(i = 0; i < N + 1; i++){
        func2d(x, y);
        x[0] += dt * y[0];
        x[1] += dt * y[1];
        t += dt;
        fprintf(fp, "%.14f\t%.14f\t%.14f\n", t, x[0], x[1]);
    }

    fclose(fp);

    return EXIT_SUCCESS;
}
