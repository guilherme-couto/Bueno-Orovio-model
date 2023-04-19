#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

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


#endif // FUNCTIONS_H