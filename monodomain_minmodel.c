/*-----------------------------------------------------
Monodomain with Minimal Ventricular model (Minimal model)
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
double theta_u = 0.3;
double theta_w = 0.015;
double theta_uminus = 0.015;
double theta_o = 0.006;
double tau_u1minus = 60.0;
double tau_u2minus = 1150.0;
double tau_uplus = 1.4506;
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
// double D = 0.0001;           // Diffusion coefficient -> cm²/ms
double chi = 1400.0;           // Surface area to volume ratio -> cm^-1
double Cm = 0.185;          // Cell capacitance per unit surface area -> uF/ (???)^2

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

// Adapted for 2nd order approximation
void thomas_algorithm_2nd(double *d, double *solution, unsigned long N, double alpha, double *c_, double *d_)
{   
    // Coefficients
    double a = -alpha;    // subdiagonal
    double b = 1 + 2 * alpha; // diagonal (1st and last row)
    double c = - 2 * alpha;    // superdiagonal
    
    // 1st: update auxiliary arrays
    c_[0] = c / b;
    d_[0] = d[0] / b;

    c = -alpha;
    
    for (int i = 1; i <= N - 2; i++)
    {
        c_[i] = c / (b - a * c_[i - 1]);
        d_[i] = (d[i] - a * d_[i - 1]) / (b - a * c_[i - 1]);
    }
    
    a = - 2 * alpha;
    d_[N - 1] = (d[N - 1] - a * d_[N - 2]) / (b - a * c_[N - 2]);

    a = -alpha;

    // 2nd: update solution
    solution[N - 1] = d_[N - 1];
    
    for (int i = N - 2; i >= 0; i--)
    {
        solution[i] = d_[i] - c_[i] * solution[i + 1];
    }
}

// Adapted for 2nd order approximation
double diffusion_i_2nd(int i, int j, int N, double **v)
{
    double result = 0.0;
    if (i == 0)
    {
        result = - 2.0*v[i][j] + 2.0*v[i + 1][j]; 
    }
    else if (i == N - 1)
    {
        result = 2.0*v[i - 1][j] - 2.0*v[i][j]; 
    }
    else
    {
        result = v[i - 1][j] - 2.0*v[i][j] + v[i + 1][j];
    }

    return result;
}

// Adapted for 2nd order approximation
double diffusion_j_2nd(int i, int j, int N, double **v)
{
    double result = 0.0;
    if (j == 0)
    {
        result = - 2.0*v[i][j] + 2.0*v[i][j + 1]; 
    }
    else if (j == N - 1)
    {
        result = 2.0*v[i][j - 1] - 2.0*v[i][j]; 
    }
    else
    {
        result = v[i][j - 1] - 2.0*v[i][j] + v[i][j + 1];
    }

    return result;
}

// Convert u to voltage in mV
double rescale_u(double u)
{
    return 85.7*u - 84;
}


/*-----------------------------------------------------
Functions of voltage variable u
-----------------------------------------------------*/
double tau_uminus(double u)
{
    double h = H(u - theta_uminus);
    return (1.0 - h) * tau_u1minus + (h * tau_u2minus);
}

double tau_wminus(double u)
{
    return tau_w1minus + (((tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus*(u - u_wminus)))) * 0.5);
}

double tau_so(double u)
{
    return tau_so1 + (((tau_so2 - tau_so1) * (1.0 + tanh(k_so*(u - u_o)))) * 0.5);
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
    if (u < theta_uminus)
    {
        return 1.0;
    }
    else
        return 0.0;
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
    return -v * H(u-theta_u) * (u-theta_u) * (u_u-u) / tau_fi;
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
    double h = H(u - theta_u);
    return (1.0 - h) * (v_inf_function(u) - v) / tau_uminus(u) - (h * v / tau_uplus);
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
int L = 2.0;            // Length of each side -> cm
double dx = 0.02;       // Spatial step -> cm
double dy = 0.02;       // Spatial step -> cm
double T = 200.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = 0.325;          // Stimulation strength -> uA/cm^2 (???)       ~52 mV   

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 4.0;         // Stimulation duration -> ms
double s1_x_limit = 0.2;            // Stimulation x limit -> cm

double t_s2_begin = 120.0;          // Stimulation start time -> ms
double stim2_duration = 2.0;        // Stimulation duration -> ms
double s2_x_max = 1.0;              // Stimulation x max -> cm
double s2_y_max = 1.0;              // Stimulation y limit -> cm
double s2_x_min = 0.0;              // Stimulation x min -> cm
double s2_y_min = 0.0;              // Stimulation y min -> cm


/*-----------------------------------------------------
Main function
-----------------------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: %s <num_threads> <delta_t (ms)> <method>\n", argv[0]);
        exit(1);
    }

    int num_threads = atoi(argv[1]);
    double dt = atof(argv[2]);
    char *method = argv[3];

    if (num_threads <= 0)
    {
        fprintf(stderr, "Number of threads must greater than 0\n");
        exit(1);
    }
    if (strcmp(method, "ADI2") != 0 && strcmp(method, "FE") != 0)
    {
        fprintf(stderr, "Method must be ADI2 (second order) or FE\n");
        exit(1);
    }

    // Number of steps
    int N = (int)(L / dx);  // Number of spatial steps (square tissue)
    int M = (int)(T / dt);  // Number of time steps

    // Variables
    double **u, **v, **w, **s;
    double **u_aux, **v_aux, **w_aux, **s_aux;

    // Auxiliary arrays
    double **r_u, **rightside, **solution, **c_, **d_;

    // Allocate memory
    u = (double **)malloc(N * sizeof(double *));
    v = (double **)malloc(N * sizeof(double *));
    w = (double **)malloc(N * sizeof(double *));
    s = (double **)malloc(N * sizeof(double *));

    u_aux = (double **)malloc(N * sizeof(double *));
    v_aux = (double **)malloc(N * sizeof(double *));
    w_aux = (double **)malloc(N * sizeof(double *));
    s_aux = (double **)malloc(N * sizeof(double *));

    r_u = (double **)malloc(N * sizeof(double *));
    rightside = (double **)malloc(N * sizeof(double *));
    solution = (double **)malloc(N * sizeof(double *));   
    c_ = (double **)malloc((N) * sizeof(double *));
    d_ = (double **)malloc((N) * sizeof(double *));

    for (int i = 0; i < N; i++)
    {
        u[i] = (double *)malloc(N * sizeof(double));
        v[i] = (double *)malloc(N * sizeof(double));
        w[i] = (double *)malloc(N * sizeof(double));
        s[i] = (double *)malloc(N * sizeof(double));

        u_aux[i] = (double *)malloc(N * sizeof(double));
        v_aux[i] = (double *)malloc(N * sizeof(double));
        w_aux[i] = (double *)malloc(N * sizeof(double));
        s_aux[i] = (double *)malloc(N * sizeof(double));

        r_u[i] = (double *)malloc(N * sizeof(double));
        rightside[i] = (double *)malloc(N * sizeof(double));
        solution[i] = (double *)malloc(N * sizeof(double));
        c_[i] = (double *)malloc(N * sizeof(double));
        d_[i] = (double *)malloc(N * sizeof(double));
    }

    double du_dt, dv_dt, dw_dt, ds_dt, diff_term = 0.0;
    double ustep, vstep, wstep, sstep, uaux, vaux, waux, saux;

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
    int x_max = s2_x_max / dx;
    int x_min = s2_x_min / dx;
    int y_max = N;
    int y_min = N - s2_y_max / dy;

    // Diffusion coefficient and phi for ADI
    // double D = ga / (chi * Cm);             // Diffusion coefficient - isotropic
    double phi = D * dt / (dx * dx);        // For Thomas algorithm - isotropic

    // Initial conditions
    int i, j;                               // Spatial indexes i for y-axis and j for x-axis
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            w[i][j] = 0.0;
            s[i][j] = 0.0;

            u_aux[i][j] = 0.0;
            v_aux[i][j] = 0.0;
            w_aux[i][j] = 0.0;
            s_aux[i][j] = 0.0;

            r_u[i][j] = 0.0;
            rightside[i][j] = 0.0;
            solution[i][j] = 0.0;
        }
    }

    // Prepare files to save data
    // Convert dt to string
    char s_dt[10];
    sprintf(s_dt, "%.03f", dt);

    // Open the file to write for complete gif
    char fname_complete[100] = "./simulation-files/mm-";
    strcat(fname_complete, method);
    strcat(fname_complete, "-");
    strcat(fname_complete, s_dt);
    strcat(fname_complete, ".txt");
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");
    int save_rate = ceil(M / 100.0);

    // Open the file to write for times
    char fname_times[100] = "./simulation-files/sim-times-";
    strcat(fname_times, method);
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
    double start_ode, finish_ode, elapsed_ode = 0.0;
    double start_pde, finish_pde, elapsed_pde = 0.0;

    start = omp_get_wtime();

    // ADI second order
    if (strcmp(method, "ADI2") == 0)
    {
        #pragma omp parallel num_threads(num_threads) default(none) \
        private(i, j, I_stim, du_dt, dv_dt, dw_dt, ds_dt, diff_term, ustep, vstep, wstep, sstep, uaux, vaux, waux, saux) \
        shared(u, v, w, s, N, M, dt, L, s1_x_limit, stim_strength, t_s1_begin, stim_duration, x_lim, t_s2_begin, stim2_duration, \
        x_max, y_max, x_min, y_min,  time,  c_, d_, u_aux, v_aux, w_aux, s_aux, phi, T, tstep, step, \
        fp_all, velocity, save_rate, fp_times, tag, \
        r_u, rightside, solution,  start_ode, finish_ode, elapsed_ode, start_pde, finish_pde, elapsed_pde)
        {
            while (step < M)
            {
                // Get time step
                tstep = time[step];

                // Start ode timer
                #pragma omp master
                {
                    start_ode = omp_get_wtime();
                }

                // Predict aux variables with explicit method
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        // Stimulus 1
                        if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration && j <= x_lim)
                        {
                            I_stim = stim_strength;
                        }
                        // Stimulus 2
                        else if (tstep >= t_s2_begin && tstep <= t_s2_begin + stim2_duration && j >= x_min && j <= x_max && i >= y_min && i <= y_max)
                        {
                            I_stim = stim_strength;
                        }
                        else 
                        {
                            I_stim = 0.0;
                        }

                        ustep = u[i][j];
                        vstep = v[i][j];
                        wstep = w[i][j];
                        sstep = s[i][j];

                        // Get du_dt, dv_dt, dw_dt and ds_dt
                        du_dt = - reaction_u(ustep, vstep, wstep, sstep) + I_stim;
                        dv_dt = dvdt(ustep, vstep);
                        dw_dt = dwdt(ustep, wstep);
                        ds_dt = dsdt(ustep, sstep);
                        
                        // Update u_aux, v_aux, w_aux and s_aux
                        uaux = ustep + ((phi * 0.5) * (diffusion_i_2nd(i, j, N, u) + diffusion_j_2nd(i, j, N, u))) + (dt * 0.5 * du_dt);
                        vaux = vstep + (dt * 0.5) * dv_dt;
                        waux = wstep + (dt * 0.5) * dw_dt;
                        saux = sstep + (dt * 0.5) * ds_dt;

                        // Update r_u for Thomas algorithm
                        du_dt = - reaction_u(uaux, vaux, waux, saux) + I_stim;
                        r_u[i][j] = 0.5 * dt * du_dt;
                        
                        // Update v, w, s
                        dv_dt = dvdt(uaux, vaux);
                        dw_dt = dwdt(uaux, waux);
                        ds_dt = dsdt(uaux, saux);

                        v[i][j] = vstep + dt * dv_dt;
                        w[i][j] = wstep + dt * dw_dt;
                        s[i][j] = sstep + dt * ds_dt;

                        // Update aux variables 
                        u_aux[i][j] = uaux;
                        v_aux[i][j] = vaux;
                        w_aux[i][j] = waux;
                        s_aux[i][j] = saux;
                    }
                }

                // Finish ode timer
                #pragma omp master
                {
                    finish_ode = omp_get_wtime();
                    elapsed_ode += finish_ode - start_ode;
                }
                
                // Update rightside for Thomas algorithm
                #pragma omp barrier
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        // Explicit diffusion in x-axis
                        diff_term = u[i][j] + (phi*0.5) * diffusion_j_2nd(i, j, N, u);
                        
                        // Update rightside
                        rightside[j][i] = diff_term + r_u[i][j];
                    }
                }

                // Start pde timer
                #pragma omp master
                {
                    start_pde = omp_get_wtime();
                }

                // Solve tridiagonal system for u
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    thomas_algorithm_2nd(rightside[i], solution[i], N, (phi*0.5), c_[i], d_[i]);
                    
                    // Update u
                    for (j = 0; j < N; j++)
                    {
                        u[j][i] = solution[i][j];
                    }
                }

                // Finish pde timer
                #pragma omp master
                {
                    finish_pde = omp_get_wtime();
                    elapsed_pde += finish_pde - start_pde;
                }
                
                // Update rightside for Thomas algorithm
                #pragma omp barrier
                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        // Explicit diffusion in y-axis
                        diff_term = u[i][j] + (phi*0.5) * diffusion_i_2nd(i, j, N, u);
                        
                        // Update rightside
                        rightside[i][j] = diff_term + r_u[i][j];
                    }
                }

                // Start pde timer
                #pragma omp master
                {
                    start_pde = omp_get_wtime();
                }
                
                // Solve tridiagonal system for u
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    thomas_algorithm_2nd(rightside[i], u[i], N, (phi*0.5), c_[i], d_[i]);
                }

                // Finish pde timer
                #pragma omp master
                {
                    finish_pde = omp_get_wtime();
                    elapsed_pde += finish_pde - start_pde;
                }

                // Save data to file
                #pragma omp master
                {
                    // Write to file
                    if (step % save_rate == 0)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            for (int j = 0; j < N; j++)
                            {
                                fprintf(fp_all, "%lf\n", rescale_u(u[i][j]));
                            }
                        }
                        fprintf(fp_times, "%lf\n", time[step]);
                    }

                    // Check S1 velocity
                    if (rescale_u(u[0][N-1]) > 20 && tag)
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
    }

    // Forward Euler
    else if (strcmp(method, "FE") == 0)
    {
        #pragma omp parallel num_threads(num_threads) default(none) \
        private(i, j, I_stim, du_dt, dv_dt, dw_dt, ds_dt, diff_term, ustep, vstep, wstep, sstep) \
        shared(u, v, w, s, N, M, dt, L, s1_x_limit, stim_strength, t_s1_begin, stim_duration, x_lim, t_s2_begin, stim2_duration, \
        x_max, y_max, x_min, y_min,  time,  c_, d_, u_aux, v_aux, w_aux, s_aux, phi, T, tstep, step, \
        fp_all, velocity, save_rate, fp_times, tag, \
        r_u, rightside, solution,  start_ode, finish_ode, elapsed_ode, start_pde, finish_pde, elapsed_pde)
        {
            while (step < M)
            {
                // Get time step
                tstep = time[step];

                // Start ode timer
                #pragma omp master
                {
                    start_ode = omp_get_wtime();
                }

                #pragma omp for collapse(2)
                for (i = 1; i < N-1; i++)
                {
                    for (j = 1; j < N-1; j++)
                    {
                        // Stimulus 1
                        if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration && j <= x_lim)
                        {
                            I_stim = stim_strength;
                        }
                        // Stimulus 2
                        else if (tstep >= t_s2_begin && tstep <= t_s2_begin + stim2_duration && j >= x_min && j <= x_max && i >= y_min && i <= y_max)
                        {
                            I_stim = stim_strength;
                        }
                        else 
                        {
                            I_stim = 0.0;
                        }

                        ustep = u[i][j];
                        vstep = v[i][j];
                        wstep = w[i][j];
                        sstep = s[i][j];

                        // Get du_dt, dv_dt, dw_dt and ds_dt
                        du_dt = - reaction_u(ustep, vstep, wstep, sstep) + I_stim;
                        dv_dt = dvdt(ustep, vstep);
                        dw_dt = dwdt(ustep, wstep);
                        ds_dt = dsdt(ustep, sstep);
                        
                        // Update u_aux, v, w and s
                        u_aux[i][j] = ustep + dt * du_dt;
                        v[i][j] = vstep + dt * dv_dt;
                        w[i][j] = wstep + dt * dw_dt;
                        s[i][j] = sstep + dt * ds_dt;
                    }
                }

                // Finish ode timer
                #pragma omp master
                {
                    finish_ode = omp_get_wtime();
                    elapsed_ode += finish_ode - start_ode;
                }

                // Boundary conditions
                #pragma omp for
                for (i = 0; i < N; i++)
                {
                    u_aux[i][0] = u_aux[i][1];
                    u_aux[i][N-1] = u_aux[i][N-2];
                    u_aux[0][i] = u_aux[1][i];
                    u_aux[N-1][i] = u_aux[N-2][i];
                }

                // Start pde timer
                #pragma omp master
                {
                    start_pde = omp_get_wtime();
                }
                
                // Diffusion
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 1; i < N-1; i++)
                {
                    for (j = 1; j < N-1; j++)
                    {   
                        u[i][j] = u_aux[i][j] + phi * (u_aux[i-1][j] - 2.0*u_aux[i][j] + u_aux[i+1][j]);
                        u[i][j] += phi * (u_aux[i][j-1] - 2.0*u_aux[i][j] + u_aux[i][j+1]);
                    }
                }

                // Finish pde timer
                #pragma omp master
                {
                    finish_pde = omp_get_wtime();
                    elapsed_pde += finish_pde - start_pde;
                }

                // Boundary conditions
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    u[i][0] = u[i][1];
                    u[i][N-1] = u[i][N-2];
                    u[0][i] = u[1][i];
                    u[N-1][i] = u[N-2][i];
                }

                // Save data to file
                #pragma omp master
                {
                    // Write to file
                    if (step % save_rate == 0)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            for (int j = 0; j < N; j++)
                            {
                                fprintf(fp_all, "%lf\n", rescale_u(u[i][j]));
                            }
                        }
                        fprintf(fp_times, "%lf\n", time[step]);
                    }

                    // Check S1 velocity
                    if (rescale_u(u[0][N-1]) > 20 && tag)
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
    }

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;

    printf("\nElapsed time = %e seconds\n", elapsed);

    // Comparison file
    FILE *fp_comp = NULL;
    fp_comp = fopen("comparison.txt", "a");
    fprintf(fp_comp, "%s  \t|\t%d threads\t|\t%.3f ms\t|\t%lf m/s\t|\t%e seconds\n", method, num_threads, dt, velocity, elapsed);

    // ODE/PDE times file
    FILE *fp_all_times = NULL;
    fp_all_times = fopen("times.txt", "a");
    fprintf(fp_all_times, "%s  \t|\t%d threads\t|\t%.3f ms\t|\t%e seconds\t|\t%e seconds\t|\t%e seconds\n", method, num_threads, dt, elapsed_ode, elapsed_pde, elapsed);

    // Close files
    fclose(fp_all);
    fclose(fp_times);
    fclose(fp_comp);
    fclose(fp_all_times);
    
    // Free alocated memory
    free(time);
    free(u);
    free(v);
    free(w);
    free(s);
    free(u_aux);
    free(v_aux);
    free(w_aux);
    free(s_aux);
    free(r_u);
    free(rightside);
    free(solution);
    free(c_);
    free(d_);

    return 0;
}