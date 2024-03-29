/*-----------------------------------------------------
Single cell with Minimal Ventricular model
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
Simulation parameters
-----------------------------------------------------*/
double T = 600.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = 1.0;          // Stimulation strength -> uA/cm^2

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 1.0;         // Stimulation duration -> ms


/*-----------------------------------------------------
Main function
-----------------------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s <delta_t (ms)>\n", argv[0]);
        exit(1);
    }

    double dt = atof(argv[1]);

    if (dt <= 0)
    {
        fprintf(stderr, "Delta_t must greater than 0\n");
        exit(1);
    }

    // Number of steps
    int M = (int)(T / dt);  // Number of time steps

    // Variables
    double u, v, w, s;

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

    // Initial conditions
    u = 0.0;
    v = 1.0;
    w = 1.0;
    s = 0.0;

    // Prepare files to save data
    // Convert dt to string
    char s_dt[10];
    sprintf(s_dt, "%.03f", dt);

    system("mkdir -p simulation-files");

    // Open the file to write for complete gif
    char fname_complete[100] = "./simulation-files/mm-";
    strcat(fname_complete, "cell");
    strcat(fname_complete, "-");
    strcat(fname_complete, s_dt);
    strcat(fname_complete, ".txt");
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");
    
    // Open the file to write for times
    char fname_times[100] = "./simulation-files/sim-times-";
    strcat(fname_times, "cell");
    strcat(fname_times, "-");
    strcat(fname_times, s_dt);
    strcat(fname_times, ".txt");
    FILE *fp_times = NULL;
    fp_times = fopen(fname_times, "w");
    
    // Timer
    double start, finish, elapsed = 0.0;

    start = omp_get_wtime();

    // Forward Euler
    while (step < M)
    {
        // Get time step
        tstep = time[step];

        // Stimulus 1
        if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration)
        {
            I_stim = stim_strength;
        }
        else 
        {
            I_stim = 0.0;
        }

        ustep = u;
        vstep = v;
        wstep = w;
        sstep = s;

        // Get du_dt, dv_dt, dw_dt and ds_dt
        du_dt = - reaction_u(ustep, vstep, wstep, sstep) + I_stim;
        dv_dt = dvdt(ustep, vstep);
        dw_dt = dwdt(ustep, wstep);
        ds_dt = dsdt(ustep, sstep);
        
        // Update u, v, w and s
        u = ustep + dt * du_dt;
        v = vstep + dt * dv_dt;
        w = wstep + dt * dw_dt;
        s = sstep + dt * ds_dt;

        // Save data to file
        // Write to file
        fprintf(fp_all, "%lf\n", rescale_u(u));
        fprintf(fp_times, "%lf\n", time[step]);
        
        // Update step
        step++;
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

    return 0;
}