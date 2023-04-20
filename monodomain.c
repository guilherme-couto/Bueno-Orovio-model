/*-----------------------------------------------------
Monodomain with Minimal Ventricular model
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
double D = 1.171;              // Diffusion coefficient -> cm²/s
double chi = 1400.0;           // Surface area to volume ratio -> cm^-1


/*-----------------------------------------------------
Simulation parameters
-----------------------------------------------------*/
int L = 2.0;            // Length of each side -> cm
double dx = 0.02;       // Spatial step -> cm
double dy = 0.02;       // Spatial step -> cm
double T = 400.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = 38;          // Stimulation strength -> uA/cm^2 

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 2.0;         // Stimulation duration -> ms
double s1_x_limit = 0.04;            // Stimulation x limit -> cm

double t_s2_begin = 310.0;          // Stimulation start time -> ms
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

    // Auxiliar arrays
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

    // Phi
    double phi = D * dt / (chi * dx * dx);        // For Thomas algorithm - isotropic

    // Initial conditions
    int i, j;                               // Spatial indexes i for y-axis and j for x-axis
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 1.0;
            w[i][j] = 1.0;
            s[i][j] = 0.0;

            u_aux[i][j] = 0.0;
            v_aux[i][j] = 1.0;
            w_aux[i][j] = 1.0;
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

    system("mkdir -p simulation-files");

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
                    if (rescale_u(u[0][N-1]) > 40.0 && tag)
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
                    if (rescale_u(u[0][N-1]) > 40.0 && tag)
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