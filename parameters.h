#ifndef PARAMETERS_H
#define PARAMETERS_H

/*-----------------------------------------------------
Model definition
https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub
-----------------------------------------------------*/

/*-----------------------------------------------------
Options for define cell type:
-----------------------------------------------------*/

// #define EPI
// #define ENDO
// #define M
// #define PB
#define TNNP


/*-----------------------------------------------------
EPI parameters
-----------------------------------------------------*/
#ifdef EPI
double u_o = 0.0;
double u_u = 1.55;
double theta_v = 0.3;
double theta_w = 0.13;
double theta_vminus = 0.006;
double theta_o = 0.006;
double tau_v1minus = 60.0;
double tau_v2minus = 1150.0;
double tau_vplus = 1.4506;
double tau_w1minus = 60.0;
double tau_w2minus = 15.0;
double k_wminus = 65.0;
double u_wminus = 0.03;
double tau_wplus = 200.0;
double tau_fi = 0.11;
double tau_o1 = 400.0;
double tau_o2 = 6.0;
double tau_so1 = 30.0181;
double tau_so2 = 0.9957;
double k_so = 2.0458;
double u_so = 0.65;
double tau_s1 = 2.7342;
double tau_s2 = 16.0;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 1.8875;
double tau_winf = 0.07;
double w_infstar = 0.94;
#endif


/*-----------------------------------------------------
ENDO parameters
-----------------------------------------------------*/
#ifdef ENDO
double u_o = 0.0;
double u_u = 1.56;
double theta_v = 0.3;
double theta_w = 0.13;
double theta_vminus = 0.2;
double theta_o = 0.006;
double tau_v1minus = 75.0;
double tau_v2minus = 10.0;
double tau_vplus = 1.4506;
double tau_w1minus = 6.0;
double tau_w2minus = 140.0;
double k_wminus = 200.0;
double u_wminus = 0.016;
double tau_wplus = 280.0;
double tau_fi = 0.1;
double tau_o1 = 470.0;
double tau_o2 = 6.0;
double tau_so1 = 40.0;
double tau_so2 = 1.2;
double k_so = 2.0;
double u_so = 0.65;
double tau_s1 = 2.7342;
double tau_s2 = 2.0;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 2.9013;
double tau_winf = 0.0273;
double w_infstar = 0.78;
#endif


/*-----------------------------------------------------
M parameters
-----------------------------------------------------*/
#ifdef M
double u_o = 0.0;
double u_u = 1.61;
double theta_v = 0.3;
double theta_w = 0.13;
double theta_vminus = 0.1;
double theta_o = 0.005;
double tau_v1minus = 80.0;
double tau_v2minus = 1.4506;
double tau_vplus = 1.4506;
double tau_w1minus = 70.0;
double tau_w2minus = 8.0;
double k_wminus = 200.0;
double u_wminus = 0.016;
double tau_wplus = 280.0;
double tau_fi = 0.078;
double tau_o1 = 410.0;
double tau_o2 = 7.0;
double tau_so1 = 91.0;
double tau_so2 = 0.8;
double k_so = 2.1;
double u_so = 0.6;
double tau_s1 = 2.7342;
double tau_s2 = 4.0;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 3.3849;
double tau_winf = 0.01;
double w_infstar = 0.5;
#endif


/*-----------------------------------------------------
PB parameters
-----------------------------------------------------*/
#ifdef PB
double u_o = 0.0;
double u_u = 1.45;
double theta_v = 0.35;
double theta_w = 0.13;
double theta_vminus = 0.175;
double theta_o = 0.006;
double tau_v1minus = 10.0;
double tau_v2minus = 1150.0;
double tau_vplus = 1.4506;
double tau_w1minus = 140.0;
double tau_w2minus = 6.25;
double k_wminus = 65.0;
double u_wminus = 0.015;
double tau_wplus = 326.0;
double tau_fi = 0.105;
double tau_o1 = 400.0;
double tau_o2 = 6.0;
double tau_so1 = 30.0181;
double tau_so2 = 0.9957;
double k_so = 2.0458;
double u_so = 0.65;
double tau_s1 = 2.7342;
double tau_s2 = 16.0;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 1.8875;
double tau_winf = 0.175;
double w_infstar = 0.9;
#endif


/*-----------------------------------------------------
TNNP parameters
-----------------------------------------------------*/
#ifdef TNNP
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
#endif

#endif // PARAMETERS_H