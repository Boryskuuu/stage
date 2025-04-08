#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

#define PI 3.141592653589793

// Generates a Gaussian wavepacket
double *gauss(double sig, double mu, double a, double *x, int N) {
    double *g = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        g[i] = a * exp(-pow(x[i] - mu, 2) / (2 * sig * sig));
    }
    return g;
}

// Computes density
void density(double rc, double complex *phi, double *x, double complex *result, int N) {
    #pragma omp parallel for
    for (int j = 0; j < N; j++) {
        result[j] = 0.0;
        for (int i = 0; i < N; i++) {
            double complex norm = phi[i] * conj(phi[i]);
            result[j] += norm * cexp(-pow(x[j] - x[i], 2) / (2 * rc * rc));
        }
        result[j] *= 1.0 / (rc * sqrt(2 * PI));
    }
}

// Hamiltonian: drift towards center
void Hcentre(int N, double **H) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            H[i][j] = 0.0;
    for (int i = 0; i < N / 2; i++)
        H[i][i + 1] = 1.0;
    for (int i = N / 2; i < N - 1; i++)
        H[i][i - 1] = 1.0;
}

// Time evolution
void timeevol(double complex *phi, double rc, double gam, double *x, double dt, double dx, int N, double **H) {
    double complex *schro = malloc(N * sizeof(double complex));
    double complex *dens = malloc(N * sizeof(double complex));
    double *dW = malloc(N * sizeof(double));
    double complex *wiener = malloc(N * sizeof(double complex));
    double complex *deterministic = malloc(N * sizeof(double complex));
    double norm = 0.0;

    // Schrödinger term
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        schro[i] = 0.0;
        for (int j = 0; j < N; j++) {
            schro[i] += -I * dt * H[i][j] * phi[j];
        }
    }

    // Compute density
    density(rc, phi, x, dens, N);

    // Generate Wiener increments and compute evolution
    #pragma omp parallel for reduction(+:norm)
    for (int i = 0; i < N; i++) {
        dW[i] = sqrt(dt) * ((double) rand() / RAND_MAX * 2.0 - 1.0);
        wiener[i] = dens[i] * dW[i] * phi[i];
        deterministic[i] = -0.5 * gam * dens[i] * dens[i] * dt * phi[i];
        phi[i] += wiener[i] + deterministic[i] + schro[i];
        norm += creal(phi[i] * conj(phi[i]));
    }

    norm = sqrt(norm * dx);
    for (int i = 0; i < N; i++) {
        phi[i] /= norm;
    }

    free(schro); free(dens); free(dW); free(wiener); free(deterministic);
}

// Initializes a symmetric double Gaussian
void doublegauss(double sigma, double center, double a, double *x, double complex *phi, double dx, int N) {
    double *g1 = gauss(sigma, center, a, x, N);
    double *g2 = gauss(sigma, -center, a, x, N);
    double norm = 0.0;

    for (int i = 0; i < N; i++) {
        phi[i] = g1[i] + g2[i];
        norm += creal(phi[i] * conj(phi[i]));
    }

    norm = sqrt(norm * dx);
    for (int i = 0; i < N; i++) {
        phi[i] /= norm;
    }

    free(g1); free(g2);
}



// Run single simulation and return collapse step
int runvariance(double sigma, double center, double a, double *x, double dx, int Nt, int N, double dt, double **H, double rc, double gam) {
    double complex *phi = malloc(N * sizeof(double complex));
    doublegauss(sigma, center, a, x, phi, dx, N);

    double *lvar = calloc(Nt, sizeof(double));
    int indvar = 0, conv = 0, stop = 0;

    for (int i = 0; i < Nt; i++) {
        timeevol(phi, rc, gam, x, dt, dx, N, H);
        double sum1 = 0.0, sum2 = 0.0;
        for (int j = 0; j < N; j++) {
            double prob = creal(phi[j] * conj(phi[j]));
            sum1 += x[j] * prob;
            sum2 += x[j] * x[j] * prob;
        }
        lvar[i] = sum2 - sum1 * sum1;
        if (i > 0 && fabs(lvar[i] - lvar[i - 1]) < 1e-8)
            conv++;
        else
            conv = 0;
        if (conv == 3 && !stop) {
            indvar = i - 3;
            stop = 1;
        }
    }

    free(phi); free(lvar);
    return indvar;
}

// Top-level simulation loop
void vartotrun(int cmax, int cfreq, double sigma, double a, double *x, double dx, int Nt, int N, double dt, double **H, double rc, double gam, double *indvar) {
    #pragma omp parallel for
    for (int c = 0; c < cmax; c++) {
        double center = (double)c / cfreq;
        int nruns = 20;
        int sum = 0;
        for (int run = 0; run < nruns; run++) {
            sum += runvariance(sigma, center, a, x, dx, Nt, N, dt, H, rc, gam);
        }
        indvar[c] = (double)sum / nruns;
    }
}

#define N 1000
#define D 25
#define Nt 500
#define dt 0.01
#define cmax 60
#define cfreq 30
#define nruns 20
#define sigma 1.0
#define a 1.0





int main() {

double lambda_csl = 1e-6;
double rc = 1.0;
double dx;
double gam;
double x[N];
double **H = malloc(N * sizeof(double *));
double indvar[cmax];


for (int i = 0; i < N; i++) {
    H[i] = malloc(N * sizeof(double));
}
    
    dx = 2.0 * D / N;
    gam = lambda_csl * pow((4.0 * M_PI * rc * rc), 1.5);

    // Init grille x
    for (int i = 0; i < N; i++) {
        x[i] = -D + i * dx;
    }

    // Init Hamiltonien
    for (int i = 0;i < N; i++) {
    for (int j = 0;j < N; j++) {
    H[i][j] = 0;
    }
    }


    // Boucle principale
    #pragma omp parallel for
    for (int c = 0; c < cmax; c++) {
        double center = (double)c / cfreq;
        double sum = 0.0;
        printf("center %.2f commencé (%d/%d)\n", center, c , cmax);
        for (int r = 0; r < nruns; r++) {
            sum += runvariance(sigma, center, a, x, dx, Nt, N, dt, H, rc, gam); 
        }
        indvar[c] = sum / nruns;
    }

    // Sauvegarde des résultats
    FILE *f = fopen("resultats.txt", "w");
    if (f == NULL) {
        printf("Erreur d'ouverture de fichier\n");
        return 1;
    }

    for (int c = 0; c < cmax; c++) {
        double center = (double)c / cfreq;
        fprintf(f, "%lf %lf\n", center, indvar[c]);
    }

    fclose(f);
    printf("Simulation terminée. Résultats écrits dans resultats.txt\n");

    return 0;
}

