#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <time.h>

#define PI 3.141592653589793

double randn (double mu, double sigma)
{
// taken from : https://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/ 
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


// Generates a Gaussian
double *gauss(double sig, double mu, double a, double *x, int N) {
double *g = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        g[i] = a * exp(-pow(x[i] - mu, 2) / (2 * sig * sig));
    }
    return g;
}

// Computes density
void wiener_func(double rc, double *x, double complex *result, int N, double dx, double dt, double gamma) {
    for (int j = 0; j < N; j++) {
        result[j] = 0.0;
        for (int i = 0; i < N; i++) {
            result[j] += cexp(-pow(x[j] - x[i], 2) / (2 * rc * rc)) * randn(0,1);
        }
        result[j] *= dx * sqrt(dt * gamma) / (rc * sqrt(2 * PI));
    }
}

// Hamiltonian: drift towards center
void Hcentre(int N, double *H) {
    for (int i = 0; i < N / 2; i++)
        H[i * N + (i+1)] = 1.0;
    for (int i = N / 2; i < N - 1; i++)
        H[i * N + (i-1)] = 1.0;
}

// Time evolution
void timeevol(double complex *phi, double rc, double gam, double *x, double dt, double dx, int N, double complex *schro, double complex *wiener) { //  double *H,
    double norm = 0.0;
/*
   // Schrödinger term
    for (int i = 0; i < N; i++) {
        schro[i] = 0.0;
        for (int j = 0; j < N; j++) {
            schro[i] += -I * dt * H[i * N + j] * phi[j];
        }
    } 
*/
    // Compute integrals
    wiener_func(rc, x, wiener, N,dx,  dt , gam);
    
    // Generate Wiener increments and compute evolution

    for (int i = 0; i < N; i++) {
    	phi[i] += wiener[i] * phi[i] ; // Wiener
	phi[i] += - gam * dt * phi[i] / (4 * rc * sqrt(PI)) ; // deterministe
	// phi[i] += schro[i] ;
        norm += creal(phi[i] * conj(phi[i]));
    }

    norm = sqrt(norm);
    for (int i = 0; i < N; i++) {
        phi[i] /= norm;
    }

}

// Initializes a symmetric double Gaussian
void doublegauss(double sigma, double center, double a, double *x, double complex *phi, double dx, int N) {
    double norm = 0.0;
    double *g1 = gauss(sigma, center, a, x, N);
    double *g2 = gauss(sigma, -center, a, x, N);

    for (int i = 0; i < N; i++) {
        phi[i] = g1[i] + g2[i] ;
        norm += creal(phi[i] * conj(phi[i]));
    }

    norm = sqrt(norm * dx);
    for (int i = 0; i < N; i++) {
        phi[i] /= norm;
    }

}



// Run single simulation and return collapse step
int runvariance(double sigma, double center, double a, double *x, double dx, int Nt, int N, double dt,  double rc, double gam) { //double *H,
    double complex *phi = malloc(N * sizeof(double complex));
    doublegauss(sigma, center, a, x, phi, dx, N);
    double complex *schro = malloc(N * sizeof(double complex));
    double *lvar = calloc(Nt, sizeof(double));
    double complex *dens = malloc(N * sizeof(double complex));
    int indvar = 0, conv = 0, stop = 0;

    for (int i = 0; i < Nt; i++) {
        timeevol(phi, rc, gam, x, dt, dx, N, schro,dens); // H,
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

void vartotrun(int cmax, int cfreq, double sigma, double a, double *x, double dx, int Nt, int N, double dt, double rc, double gam, double *indvar) { //  double *H,
  int nruns = 20 ;
    #pragma omp parallel for
    for (int c = 0; c < cmax; c++) {
        double center = (double)c / cfreq;
        double sum = 0;
        printf("center %.2f commencé (%d/%d)\n", center, c , cmax);
        for (int run = 0; run < nruns; run++) {
            sum += runvariance(sigma, center, a, x, dx, Nt, N, dt, rc, gam); // H,
        }
        indvar[c] = sum / nruns;
    }
}

int singlerun(double sigma, double a, double *x, double dx, int Nt, int N, double dt, double rc, double gam){ //  double *H,

double center = 4 ;
double complex *phi = malloc(N * sizeof(double complex));
doublegauss(sigma, center, a, x, phi, dx, N);
double complex *schro = malloc(N * sizeof(double complex));
double complex *dens = malloc(N * sizeof(double complex));
int nstock = 5 ;
int etape_evol[] = {0,10,20,30,40} ;

double *phistock = malloc(N*nstock*sizeof(double)) ;

    for (int i = 0; i < Nt; i++)
    {
        timeevol(phi, rc, gam, x, dt, dx, N, schro,dens); // H,
        for(int j = 0 ; j < nstock ; j++){
            if(etape_evol[j] == i){
                for(int k = 0; k < N;k++){
                    phistock[j*N+k] = creal(phi[k] * conj(phi[k])) ;
                }
            }
        }
    }

    // Sauvegarde des résultats
    FILE *f = fopen("singlerun.txt", "w");
    if (f == NULL) {
        printf("Erreur d'ouverture de fichier\n");
        return 1;
    }
    
    fprintf(f, "x \t\t");
    for (int j = 0 ; j < nstock ; j++){
    fprintf(f, "t = %d \t\t",etape_evol[j]);
    }
    fprintf(f,"\n") ;
    
    for (int i = 0 ; i < N ; i++){
    fprintf(f, "    %lf \t",x[i]);
    for (int j = 0 ; j < nstock ; j++){
    fprintf(f, "    %lf \t",phistock[j*N + i]);
    }
    fprintf(f,"\n") ;
    }
    

    fclose(f);
    printf("Simulation terminée. Résultats écrits dans singlerun.txt\n");

free(phi); free(schro) ; free(dens) ; free(phistock) ; 

return 0 ;
}



#define N 1000
#define D 25
#define Nt 500
#define dt 0.01
#define cmax 60
#define cfreq 30
#define sigma 1.0
#define a 1.0






int main() {

double lambda_csl = 1e-16;
double rc = 1.0;
double dx;
double gam;
double x[N];
// double *H = calloc(N * N, sizeof(double)); // calloc : init à 0
double indvar[cmax];

    
    dx = 2.0 * D / N;
    gam = lambda_csl * pow((4.0 * M_PI * rc * rc), 1.5);

    // Init grille x
    for (int i = 0; i < N; i++) {
        x[i] = -D + i * dx;
    }

singlerun(sigma,a,x,dx,Nt, N, dt, rc, gam) ;// H,


/* /// PROBABLEMENT PLEIN DE CONNERIES ICI C'EST CHAT QUI A ÉCRIS CE CON
vartotrun(cmax, cfreq, sigma, a, x, dx, Nt, N, dt, rc, gam, indvar) ; // H, 

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
*/

    return 0;
}

