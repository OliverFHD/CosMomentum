#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstring>


namespace constants{

/***** Enums *****/

enum SHOT_NOISE_MODUS {BA0A1, BR, INVALID_SHOT_NOISE};
enum PDF_MODUS {BERNARDEAU, LOGNORMAL, LOGNORMAL_FIX_D0, GAUSSIAN, INVALID_PDF}; 
enum BINNING {LIN, LOG};
enum ON_OFF {ON, OFF};
enum INITIALISATION {INITIALISED, UNINITIALISED};

const char *SHOT_NOISE_types[] = {"BA0A1", "BR", "INVALID_SHOTNOISE"};
const char *PDF_types[] = {"BERNARDEAU", "LOGNORMAL", "LOGNORMAL_FIX_D0", "GAUSSIAN", "INVALID_PDF"};
const char *BINNING_types[] = {"LIN", "LOG"};
const char *ON_OFF_types[] = {"ON", "OFF"};
const char *INITIALISATION_types[] = {"INITIALISED", "UNINITIALISED"};
  

/***** General constants *****/

static double eulers_constant = exp(1.0);
static double pi = 3.14159265;
static double pi2 = 2.0*pi;
static double sqrtPi_times_2 = 2.0*sqrt(pi);
static double pi_sq = pi*pi;
static double one_over_2_pi_sq = 1.0/(2.0*pi_sq); // = 1.0/(2.0*pi^2)
static double sqrt2 = 1.414213562373095;
static double sqrt2_plus_1 = 2.414213562373095;
static double rad = pi/(3600*180);
static double c_in_si =  2.99792458e8;
static double c_over_e5 = c_in_si/1.0e5;
static double c_over_e5_sq = c_over_e5*c_over_e5;

  /*
    TOM WAS HERE
    - c_over_e5 is divided into many things. To get rid of devisions we can calculate the inverse ahead of time and just multiply this on.
   */

static double inverse_c_over_e5 = 1.0/c_over_e5;
static double arcmin = 2.908882087e-4;

/***** epsilon to judge whether something is 0 (zero) *****/

static double eps_from_0 = 1.0e-10;

/***** Variables controlling the precision finding roots by Newton's method or my nested intervals *****/

static double eps_Newton = 1.0e-5;
static double eps_nested = 1.0e-5;


/***** Max ell cut off for 2D power spectra *****/

static int ell_max = 20000;


/***** Variable controlling the order of interpolation polynomials *****/

//static int order_of_interpolation = 6;
static int order_of_interpolation = 4;


/***** Variable controlling maximal stepsize when integration expansion history *****/

static double maximal_da = 0.001;


/***** Variable controlling the order of interpolating polynomials for generating functions *****/

static int generating_function_coeff_order = 19;
//static int generating_function_coeff_order = 15;
//static int generating_function_coeff_order = 10;
//static int generating_function_coeff_order = 7;

/***** Variable controlling for how many delta values the PDF is computed *****/

static int N_delta_values_for_PDFs = 600;
//static int N_delta_values_for_PDFs = 200;
static double max_contrast = 10.0;
//static int N_delta_values_for_PDFs = 2000;


/***** Variable controlling range over which PDF(kappa_CMB) is computed *****/
static double kappa_CMB_min = -0.15;
static double kappa_CMB_max = 0.17;


/***** Variables controlling the precision of gsl integrators *****/

static double gsl_hstart_relative = 1e-6; // Would correspond to 1000000 steps of integration.
static double gsl_eps_absolute = 1e-6; // Hard to define a reasonable default here... Try to use gsl_eps_relative * some expectation of the dynamic range instead
static double gsl_eps_relative = 1e-10; // gsl_eps_relative/gsl_hstart_relative \approx 0.00001% .

/***** Variables controlling the precision of all integrals over the 3D-power spectrum in Matter.cpp *****/

//static int number_of_k = 2048;
// values of wave numbers in h/Mpc
//static double minimal_wave_number = 3.336e-5;
//static double maximal_wave_number = 336.0;

static int number_of_k = 4096;
//static double minimal_wave_number = 1.0e-3;
//static double maximal_wave_number = 1.0e3;
//static double minimal_wave_number = 0.158671E-04;
static double minimal_wave_number = 0.0001;
//static double maximal_wave_number = 0.518697E+02;
//static double maximal_wave_number = 3360.0;
static double maximal_wave_number = 4962.49;
static double high_k_cutoff = maximal_wave_number;

static double log_minimal_wave_number = log(minimal_wave_number);
static double log_maximal_wave_number = log(maximal_wave_number);

//static double product_of_kmax_and_R = 1.0e2;//1.0e3;
static double product_of_kmax_and_R = 1.0e3;//1.0e3;

// values of wave numbers in H_0/c
static double high_k_cutoff_in_H0_units = high_k_cutoff*c_over_e5;
static double minimal_wave_number_in_H0_units = minimal_wave_number*c_over_e5;
static double maximal_wave_number_in_H0_units = maximal_wave_number*c_over_e5;

static double log_minimal_wave_number_in_H0_units = log(minimal_wave_number_in_H0_units);
static double log_maximal_wave_number_in_H0_units = log(maximal_wave_number_in_H0_units);


/***** Variable controlling how close the variance of delta must be to 1.0 when finding k_NL (for halofit)  *****/

static double sig_sq_precision = 0.0001;

/***** Variable controlling the maximal relative bin width when integrating kappa(theta) to get gamma(theta) *****/

static double max_relative_bin_width = 0.02;


/***** minimum number of steps when performing projections along the LOS (too prevent input histograms being to broadly binned) *****/
/***** NOTE: this should be changed in the future, to integrate redshift histograms with variable stepsizes! ****/

static int minimal_n_w = 300;
static double maximal_dw = 0.01;
static double z_last_scattering = 1090.30; // from https://arxiv.org/pdf/1807.06209.pdf , TT-only

};





namespace overloads{

  std::string to_string(int val) {
    return std::to_string((long long) (val));
  }

  std::string to_string(double val) {
    return std::to_string((long double) (val));
  }

}
