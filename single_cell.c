/*-----------------------------------------------------
Single cell with Minimal Ventricular model (Minimal model / Bueno-Orovio model)
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

/*-----------------------------------------------------
Model definition
https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub
-----------------------------------------------------*/

/*---------
Possible adaptations
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052234
-------*/

/*-----------------------------------------------------
Model parameters
Adjusted to TNNP
-----------------------------------------------------*/
double u_o = 0.0;
double u_u = 1.58;
double theta_v = 0.3;
double theta_w = 0.015;
double theta_vminus = 0.015;
double theta_o = 0.006;
double tau_v1minus = 60.0;
double tau_v2minus = 1150.0;
double tau_vplus = 1.4506;
double tau_w1minus = 70.0;
double tau_w2minus = 20.0;
double k_wminus = 65.0;
double u_wminus = 0.03;
double tau_wplus = 280.0;
double tau_fi = 0.11;
double tau_o1 = 6.0;
double tau_o2 = 6.0;
double tau_so1 = 43.0;
double tau_so2 = 0.2;
double k_so = 2.0;
double u_so = 0.65;
double tau_s1 = 2.7342;
double tau_s2 = 3.0;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 2.8723;
double tau_winf = 0.07;
double w_infstar = 0.94;
double D = 1.171e-3;           // Diffusion coefficient -> cm²/ms

/*-----------------------------------------------------
Auxiliary functions
-----------------------------------------------------*/
// Standard Heaviside function
double H(double x)
{
    if (x > 0.0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

// Convert u to voltage in mV
double rescale_u(double u)
{
    return 85.7*u - 84.0;
}


/*-----------------------------------------------------
Functions of voltage variable u
-----------------------------------------------------*/
double tau_vminus(double u)
{
    double h = H(u - theta_vminus);
    return (1.0 - h) * tau_v1minus + (h * tau_v2minus);
}

double tau_wminus(double u)
{
    return tau_w1minus + (((tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus*(u - u_wminus)))) * 0.5);
}

double tau_so(double u)
{
    return tau_so1 + (((tau_so2 - tau_so1) * (1.0 + tanh(k_so*(u - u_so)))) * 0.5);
}

double tau_s(double u)
{
    double h = H(u - theta_w);
    return (1.0 - h) * tau_s1 + (h * tau_s2);
}

double tau_o(double u)
{
    double h = H(u - theta_o);
    return (1.0 - h) * tau_o1 + (h * tau_o2);
}

// Infinity values
double v_inf_function(double u)
{
    if (u < theta_vminus)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

double w_inf_function(double u)
{
    double h = H(u - theta_o);
    return (1.0 - h) * (1.0 - (u/tau_winf)) + (h * w_infstar);
}


/*-----------------------------------------------------
Currents functions
-----------------------------------------------------*/
double J_fi(double u, double v)
{
    return -v * H(u-theta_v) * (u-theta_v) * (u_u-u) / tau_fi;
}

double J_so(double u)
{
    double h = H(u-theta_w);
    return ((u-u_o) * (1.0 - h) / tau_o(u)) + (h / tau_so(u));
}

double J_si(double u, double w, double s)
{
    return - H(u-theta_w) * w * s / tau_si;
}


/*-----------------------------------------------------
Differential equations for each variable
-----------------------------------------------------*/
double reaction_u(double u, double v, double w, double s)
{
    return (J_fi(u, v) + J_so(u) + J_si(u, w, s));
}

double dvdt(double u, double v)
{
    double h = H(u - theta_v);
    return (1.0 - h) * (v_inf_function(u) - v) / tau_vminus(u) - (h * v / tau_vplus);
}

double dwdt(double u, double w)
{
    double h = H(u - theta_w);
    return (1.0 - h) * (w_inf_function(u) - w) / tau_wminus(u) - (h * w / tau_wplus);
}

double dsdt(double u, double s)
{
    return (((1.0 + tanh(k_s*(u - u_s))) * 0.5) - s) / tau_s(u);
}


/*-----------------------------------------------------
Simulation parameters
-----------------------------------------------------*/
double T = 600.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = 1.0;          // Stimulation strength -> uA/cm^2 (???)       ~52 mV   

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 2.0;         // Stimulation duration -> ms


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
    #pragma omp parallel num_threads(num_threads) default(none) \
    private(I_stim, du_dt, dv_dt, dw_dt, ds_dt, ustep, vstep, wstep, sstep) \
    shared(u, v, w, s, M, dt, stim_strength, t_s1_begin, stim_duration, \
    time, T, tstep, step, \
    fp_all, fp_times)
    {
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
            
            // Update u_aux, v, w and s
            u = ustep + dt * du_dt;
            v = vstep + dt * dv_dt;
            w = wstep + dt * dw_dt;
            s = sstep + dt * ds_dt;

            // Save data to file
            #pragma omp master
            {
                // Write to file
                fprintf(fp_all, "%lf\n", rescale_u(u));
                fprintf(fp_times, "%lf\n", time[step]);
                
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

    return 0;
}