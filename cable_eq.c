/*-----------------------------------------------------
Cable Equation with Minimal Ventricular model
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

/*-----------------------------------------------------
Model definition
https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub
-----------------------------------------------------*/

/*-----------------------------------------------------
Possible adaptations
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052234
-----------------------------------------------------*/

#include "parameters.h"
#include "functions.h"

/*-----------------------------------------------------
Model parameters
-----------------------------------------------------*/
double D = 1.171;           // Diffusion coefficient -> cm²/s
double chi = 1400.0;           // Surface area to volume ratio -> cm^-1


/*-----------------------------------------------------
Simulation parameters
-----------------------------------------------------*/
int L = 2.0;            // Length of each side -> cm
double dx = 0.02;       // Spatial step -> cm
double T = 500.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = 1.0;          // Stimulation strength -> uA/cm^2

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 2.0;         // Stimulation duration -> ms
double s1_x_limit = 0.2;            // Stimulation x limit -> cm


/*-----------------------------------------------------
Main function
-----------------------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage: %s <num_threads> <delta_t (ms)>\n", argv[0]);
        exit(1);
    }

    int num_threads = atoi(argv[1]);
    double dt = atof(argv[2]);

    if (num_threads <= 0)
    {
        fprintf(stderr, "Number of threads must greater than 0\n");
        exit(1);
    }

    // Number of steps
    int N = (int)(L / dx);  // Number of spatial steps (square tissue)
    int M = (int)(T / dt);  // Number of time steps

    // Variables
    double *u, *v, *w, *s;

    // Allocate memory
    u = (double *)malloc(N * sizeof(double));
    v = (double *)malloc(N * sizeof(double));
    w = (double *)malloc(N * sizeof(double));
    s = (double *)malloc(N * sizeof(double));

    double du_dt, dv_dt, dw_dt, ds_dt;
    double ustep, vstep, wstep, sstep;

    int step = 0;
    double tstep = 0.0;
    double *time = (double *)malloc(M * sizeof(double));
    int n;  // Time index
    for (n = 0; n < M; n++)
    {
        time[n] = n * dt;
    }

    // Stim variables
    double I_stim = 0.0;
    int x_lim = s1_x_limit / dx;

    // Diffusion coefficient
    double phi = D * dt / (chi * dx * dx);        // For Thomas algorithm - isotropic

    // Initial conditions
    int i;                         
    for (i = 0; i < N; i++)
    {
        u[i] = 0.0;
        v[i] = 1.0;
        w[i] = 1.0;
        s[i] = 0.0;
    }

    // Prepare files to save data
    // Convert dt to string
    char s_dt[10];
    sprintf(s_dt, "%.03f", dt);

    system("mkdir -p simulation-files");
    
    // Open the file to write for complete gif
    char fname_complete[100] = "./simulation-files/mm-";
    strcat(fname_complete, "cable-eq");
    strcat(fname_complete, "-");
    strcat(fname_complete, s_dt);
    strcat(fname_complete, ".txt");
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");
    int save_rate = ceil(M / 150.0);

    // Open the file to write for times
    char fname_times[100] = "./simulation-files/sim-times-";
    strcat(fname_times, "cable-eq");
    strcat(fname_times, "-");
    strcat(fname_times, s_dt);
    strcat(fname_times, ".txt");
    FILE *fp_times = NULL;
    fp_times = fopen(fname_times, "w");

    // For velocity
    bool tag = true;
    double velocity = 0.0;
    
    // Timer
    double start, finish, elapsed = 0.0;

    start = omp_get_wtime();

    // Forward Euler
    #pragma omp parallel num_threads(num_threads) default(none) \
    private(i, I_stim, du_dt, dv_dt, dw_dt, ds_dt, ustep, vstep, wstep, sstep) \
    shared(u, v, w, s, N, M, dt, L, s1_x_limit, stim_strength, t_s1_begin, stim_duration, x_lim, \
    time, phi, T, tstep, step, \
    fp_all, velocity, save_rate, fp_times, tag)
    {
        while (step < M)
        {
            // Get time step
            tstep = time[step];

            #pragma omp for
            for (i = 1; i < N-1; i++)
            {
                // Stimulus 1
                if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration && i <= x_lim)
                {
                    I_stim = stim_strength;
                }
                else 
                {
                    I_stim = 0.0;
                }

                ustep = u[i];
                vstep = v[i];
                wstep = w[i];
                sstep = s[i];

                // Get du_dt, dv_dt, dw_dt and ds_dt
                du_dt = - reaction_u(ustep, vstep, wstep, sstep) + I_stim;
                dv_dt = dvdt(ustep, vstep);
                dw_dt = dwdt(ustep, wstep);
                ds_dt = dsdt(ustep, sstep);
                
                // Update u_aux, v, w and s
                u[i] = ustep + dt * du_dt + phi * (u[i-1] - 2.0*u[i] + u[i+1]);
                v[i] = vstep + dt * dv_dt;
                w[i] = wstep + dt * dw_dt;
                s[i] = sstep + dt * ds_dt;
            }

            // Boundary conditions
            #pragma omp for nowait
            for (i = 0; i < N; i++)
            {
                u[0] = u[1];
                u[N-1] = u[N-2];
            }

            // Save data to file
            #pragma omp master
            {
                // Write to file
                if (step % save_rate == 0)
                {
                    for (int i = 0; i < N; i++)
                    {
                        fprintf(fp_all, "%lf\n", rescale_u(u[i]));
                    }
                    fprintf(fp_times, "%lf\n", time[step]);
                }

                // Check S1 velocity
                if (rescale_u(u[N-1]) > 40 && tag)
                {
                    velocity = ((10*(L - s1_x_limit)) / (time[step]));
                    printf("S1 velocity: %lf\n", velocity);
                    tag = false;
                }
            } 
            
            // Update step
            #pragma omp master
            {
                step++;
            }
            #pragma omp barrier 
        }
    } 
    
    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;

    printf("\nElapsed time = %e seconds\n", elapsed);

    // Close files
    fclose(fp_all);
    fclose(fp_times);
    
    // Free alocated memory
    free(time);
    free(u);
    free(v);
    free(w);
    free(s);

    return 0;
}