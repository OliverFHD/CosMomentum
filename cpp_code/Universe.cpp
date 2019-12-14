#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "cosmology_utils.cpp"
#include "Universe.h"
#include "Universe_integrands.cpp"
#include "Universe_analytic_special_cases.cpp"
#include "Universe_initial_conditions.cpp"
 
using namespace std;
using namespace constants;


/* _____ Content Universe.cpp
 * 
 * _____ 1. Initialization
 * 
 * ..... 1.1 Constructor
 * ..... 1.2 Destructor
 * ..... 1.3 estimate_initial_step_size
 * ..... 1.4 set_background_cosmology
 * 
 * _____ 2. Cosmological Functions
 * 
 * ..... 2.1 hubble_from_Friedmann
 * ..... 2.2 hubble_prime_from_Friedmann
 * 
 * _____ 3. Return Values
 * 
 * ..... 3.1  return_eta    
 * ..... 3.3  return_cosmology             
 * ..... 3.4  f_k                          
 * ..... 3.5  rho_m_of_a                   
 * ..... 3.6  rho_r_of_a                   
 * ..... 3.7  rho_L_of_a                   
 * ..... 3.8  w_L_of_a                     
 * ..... 3.9  a_at_eta                     
 * ..... 3.10 H_at_eta                     
 * ..... 3.11 H_prime_at_eta               
 * ..... 3.12 eta_at_a                     
 * 
 * _____ 4. Output and Checks
 * 
 * ..... 4.1 print_background_cosmology
 * 
 */


/*******************************************
 *******************************************
 **__________ 1. INITIALIZATION __________**
 *******************************************
 ******************************************* 
 *                                         *
 * ..... 1.1 Constructor                   *
 * ..... 1.2 Destructor                    *
 * ..... 1.3 estimate_initial_step_size    *
 * ..... 1.4 set_background_cosmology      *
 *                                         *
 *******************************************
 *******************************************/


/*******************************************************************************************************************************************************
 * 1.1 Constructor
 * Description:
 *  - Constructor of class Universe
 * Arguments:
 *  - cosmo: struct containing cosmological parameters
 *  - t_start: initial value of time parameter
 *  - a_min: minimal scale faction (= initial scale factor for expanding universe)
 *  - a_max: maximal scale faction
 *  - expand_or_collapse: 1 == expanding universe, 0 == collapsing universe
 * 
*******************************************************************************************************************************************************/


Universe::Universe(cosmological_model cosmo, double a_min, double a_max, int expand_or_collapse){
  
  this->expansion_or_collapse = expand_or_collapse;
  
  {
    using namespace error_handling;
    string error_message;
    
    error_message = "Omega_r = " + to_string(cosmo.Omega_r)+" should be >= 0.0 .";
    test_value_of_double(cosmo.Omega_r, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "theta_27 = " + to_string(cosmo.theta_27)+" should be > 0.0 .";
    test_value_of_double(cosmo.theta_27, 0.0, LARGER, error_message);
    
    error_message = "Omega_m = " + to_string(cosmo.Omega_m)+" should be >= 0.0 .";
    test_value_of_double(cosmo.Omega_m, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "Omega_L = " + to_string(cosmo.Omega_L)+" should be >= 0.0 .";
    test_value_of_double(cosmo.Omega_L, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "Omega_b = " + to_string(cosmo.Omega_b)+" should be >= 0.0 .";
    test_value_of_double(cosmo.Omega_b, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "Omega_m = " + to_string(cosmo.Omega_m)+" should be greater or equal to Omega_b = " + to_string(cosmo.Omega_b) +" .";
    test_value_of_double(cosmo.Omega_m, cosmo.Omega_b, LARGER_EQUAL, error_message);
    
    error_message = "h_100 = " + to_string(cosmo.h_100)+" should be > 0.0 .";
    test_value_of_double(cosmo.h_100, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "sigma_8 = " + to_string(cosmo.sigma_8)+" should be > 0.0 .";
    test_value_of_double(cosmo.sigma_8, 0.0, LARGER_EQUAL, error_message);
  }

  this->cosmology = cosmo;
  
  if(expand_or_collapse == 1){
    this->a_initial = a_min;
    this->a_final = a_max;
  }
  else{
    this->a_initial = a_max;
    this->a_final = a_min;
  }
  
  this->set_initial_conditions();
  this->set_background_cosmology();
  this->print_background_cosmology("expansion_history.dat");

}

/*******************************************************************************************************************************************************
 * 1.2 Destructor
 * Description:
 * 
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

Universe::~Universe(){
  //cout << "\nThe Universe was destroyed - probably by a protouniverse intruding our universe through a gap in space time.\n(see Star Trek: Deep Space 9)\n\n";
}

/*******************************************************************************************************************************************************
 * 1.4 set_background_cosmology
 * Description:
 *  - Integrates Friedmann equations to give
 *    a(eta)
 *    H(eta)
 *    H'(eta)
 *    t(eta)
 *    where eta is conformal time and H is the conformal expansion rate a'/a.
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Universe::set_background_cosmology(){
  
  integration_parameters_Universe params;
  params.pointer_to_Universe = this;
  integration_parameters_Universe * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {scale_factor_gsl, scale_factor_gsl_jac, 3, (void *) pointer_to_params};

  // NOTE: this function is supposed to hide the blinding mechanism uppon release
  this->set_number_of_time_steps(30000);
  
  this->a.resize(this->number_of_time_steps);
  this->H.resize(this->number_of_time_steps);
  this->eta.resize(this->number_of_time_steps);
  this->t.resize(this->number_of_time_steps);
  this->H_prime.resize(this->number_of_time_steps);
  
  this->a[0] = this->a_initial;
  this->eta[0] = this->eta_initial;
  this->H[0] = this->hubble_from_Friedmann(this->a_initial);
  this->H_prime[0] = this->hubble_prime_from_Friedmann(this->a_initial);
  this->t[0] = this->t_initial;
  
  
  double da = (this->a_final - this->a_initial)/double(this->number_of_time_steps - 1);
  double a_i = a_initial;
  double a_f = a_initial + da;
  double y[3] = { this->eta[0], this->H[0] , this->t[0]};
  
  /* 
   * 
   * from https://www.gnu.org/software/gsl/manual/html_node/Stepping-Functions.html :
   * 
   * Step Type: gsl_odeiv2_step_rkf45
   * Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator.
   * 
   * Step Type: gsl_odeiv2_step_rk8pd
   * Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.
   * 
   */
  
  // Most quantities are of order 1.0 at their maximum. Only H can be bigger in the early universe.
  double hstart = 1.0*constants::gsl_hstart_relative;
  double eps_absolute = constants::gsl_eps_relative*max(1.0, this->H[0]);
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, eps_absolute, constants::gsl_eps_relative);

  for (int i = 1; i < this->number_of_time_steps; i++){
    this->a[i] = a_f;
    int status = gsl_odeiv2_driver_apply(d, &a_i, a_f, y);
    
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }

    this->eta[i] = y[0];
    this->H[i] = y[1];
    this->t[i] = y[2];
    this->H[i] = this->hubble_from_Friedmann(this->a[i]);
    this->H_prime[i] = this->hubble_prime_from_Friedmann(this->a[i]);
    a_f += da;
  }

  gsl_odeiv2_driver_free(d);
}


/***********************************************
 ***********************************************
 **________ 2. COSMOLOGICAL FUNCTIONS ________**
 ***********************************************
 *********************************************** 
 *                                             *
 * ..... 2.1 hubble_from_Friedmann             *
 * ..... 2.2 hubble_prime_from_Friedmann                      *
 *                                             *
 ***********************************************
 ***********************************************/

/*******************************************************************************************************************************************************
 * 2.1 hubble_from_Friedmann
 * Description:
 *  - Gives conformal expansion rate H at given scale factor. Within the code, this function is mainly used
 *    to set the initial value of H when integrating the Friedmann equations.
 * Arguments:
 *  - a_start: scale factor at which H shall be computed
 * 
*******************************************************************************************************************************************************/

 double Universe::hubble_from_Friedmann(double a_start){
  
  double rho_m = rho_m_of_a(a_start);
  double rho_l = rho_L_of_a(a_start);
  double rho_r = rho_r_of_a(a_start);

  return sqrt(a_start*a_start*(rho_m+rho_r+rho_l)+this->cosmology.Omega_k);
  
}

/*******************************************************************************************************************************************************
 * 2.2 hubble_prime_from_Friedmann
 *  - Gives derivative of conformal expansion rate, H', at given scale factor. Within the code, this function
 *    is mainly used to set the initial value of H' when integrating the Friedmann equations.
 * Arguments:
 *  - a_start: scale factor at which H' shall be computed
 * 
*******************************************************************************************************************************************************/

double Universe::hubble_prime_from_Friedmann(double scale){
  
  double rho_m = rho_m_of_a(scale);
  double rho_l = rho_L_of_a(scale);
  double rho_r = rho_r_of_a(scale);
  double w = this->w_L_of_a(scale);

  return -0.5*(rho_m + 2*rho_r+(1.0+3.0*w)*rho_l)*scale*scale;
  
}


/*******************************************
 *******************************************
 **__________ 3. RETURN VALUES ___________**
 *******************************************
 ******************************************* 
 *                                         *
 * ..... 3.1  return_eta_initial           *
 * ..... 3.2  return_eta_final             *
 * ..... 3.3  return_cosmology             *
 * ..... 3.4  f_k                          *
 * ..... 3.5  rho_m_of_a                   *
 * ..... 3.6  rho_r_of_a                   *
 * ..... 3.7  rho_L_of_a                   *
 * ..... 3.8  w_L_of_a                     *
 * ..... 3.9  a_at_eta                     *
 * ..... 3.10 H_at_eta                     *
 * ..... 3.11 H_prime_at_eta               *
 * ..... 3.12 eta_at_a                     *
 *                                         *
 *******************************************
 *******************************************/


/*******************************************************************************************************************************************************
 * 3.1 return_eta_initial
 * Description:
 *  Returns i-th element of the array containing the values of conformal time.
 * Arguments:
 *  - i: index, at which array element shall be returned
 * 
*******************************************************************************************************************************************************/

double Universe::return_eta_initial(){
    
  return this->eta[0];
  
}

double Universe::return_a_initial(){
    
  return this->a[0];
  
}


/*******************************************************************************************************************************************************
 * 3.2 return_eta_final
 * Description:
 *  Returns i-th element of the array containing the values of conformal time.
 * Arguments:
 *  - i: index, at which array element shall be returned
 * 
*******************************************************************************************************************************************************/

double Universe::return_eta_final(){
    
  return this->eta[this->eta.size()-1];
  
}

double Universe::return_a_final(){
    
  return this->a[this->a.size()-1];
  
}

/*******************************************************************************************************************************************************
 * 3.3 return_cosmology
 * Description:
 *  - Returns struct with cosmological parameters.
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

cosmological_model Universe::return_cosmology(){
  
  return this->cosmology;
  
}

/*******************************************************************************************************************************************************
 * 3.4 f_k
 * Description:
 *  - Co-moving angular diameter distance.
 * Arguments:
 *  - w: co-moving distance at which f_k shall be evaluated
 * 
*******************************************************************************************************************************************************/

double Universe::f_k(double w){
  double omega_0 = 1.0-this->cosmology.Omega_k;
  
  if(abs(1.0-omega_0) < eps_from_0)
    return w;

  if(omega_0 > 1.0)
    return sin(w);
  
  return sinh(w);

  
}


/*******************************************************************************************************************************************************
 * 2.5 rho_m_of_a
 * Description:
 *  - Matter density at scale factor a in units of present day critical density.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/

double Universe::rho_m_of_a(double a){
  
  return this->cosmology.Omega_m/pow(a, 3.0);
  
}



/*******************************************************************************************************************************************************
 * 2.6 rho_r_of_a
 * Description:
 *  - Radiation density at scale factor a in units of present day critical density.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/


double Universe::rho_r_of_a(double a){
  
  return this->cosmology.Omega_r/pow(a, 4.0);
  
}


/*******************************************************************************************************************************************************
 * 2.7 rho_L_of_a
 * Description:
 *  - Dark energy density at scale factor a in units of present day critical density.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/


double Universe::rho_L_of_a(double a){
  
  double w = this->w_L_of_a(a);
  
  return this->cosmology.Omega_L/pow(a, 3.0*(1.0+w));
  
}


/*******************************************************************************************************************************************************
 * 2.8 w_L_of_a
 * Description:
 * - Equation of State parameter of Dark energy. 
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/

double Universe::w_L_of_a(double a){
  
  double w1 = this->cosmology.w1;
  if(w1 == 0.0)
    return this->cosmology.w0;
  
  return this->cosmology.w0 + this->cosmology.w1*(1.0-a); // Chevallier-Polarski-Linde
  //return this->cosmology.w0 + this->cosmology.w1*(1.0-a)/a; // Linear-redshif
  //return this->cosmology.w0 + this->cosmology.w1*(1.0-a)*a; // Jassal-Bagla-Padmanabhan
  
}




/*******************************************************************************************************************************************************
 * 3.9 a_at_eta
 * Description:
 *  - Gives value of scale factor at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double Universe::a_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->a, constants::order_of_interpolation);
  
}



/*******************************************************************************************************************************************************
 * 3.9 a_at_eta
 * Description:
 *  - Gives value of scale factor at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double Universe::t_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->t, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 3.10 H_at_eta
 * Description:
 *  - Gives value of conformal expansion rate at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double Universe::H_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->H, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 3.11 H_prime_at_eta
 * Description:
 *  - Gives derivative of conformal expansion rate at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double Universe::H_prime_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->H_prime, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 3.12 eta_at_a
 * Description:
 *  - Gives value of conformal time at given scale factor.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/

double Universe::eta_at_a(double a){
  
  return interpolate_neville_aitken(a, &this->a, &this->eta, constants::order_of_interpolation);
  
}


/*******************************************************************************************************************************************************
 * 3.13 eta_at_a
 * Description:
 *  - Gives value of conformal time at given scale factor.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/


vector<vector<double> > Universe::return_background_expansion(){
  
  int n_column = 5;
  int time_steps = this->a.size();
  vector<vector<double> > background_expansion(n_column, vector<double>(time_steps, 0.0));
  
  background_expansion[0] = this->t;
  background_expansion[1] = this->eta;
  background_expansion[2] = this->a;
  background_expansion[3] = this->H;
  background_expansion[4] = this->H_prime;
  
  return background_expansion;
  
}

vector<vector<double> > Universe::return_background_expansion(int conformal_time_steps){
  
  int n_column = 5;
  vector<vector<double> > background_expansion(n_column, vector<double>(conformal_time_steps, 0.0));
  
  double a_max = this->a[this->a.size()-1];
  double a_min = this->a[0];
  double da = (a_max - a_min)/double(conformal_time_steps-1);
  double a;
  double e;
  
  for(int i = 0; i < conformal_time_steps; i++){
    a = a_min + double(i)*da;
    e = eta_at_a(a);
    background_expansion[0][i] = this->t_at_eta(e);
    background_expansion[1][i] = e;
    background_expansion[2][i] = a;
    background_expansion[3][i] = this->H_at_eta(e);
    background_expansion[4][i] = this->H_prime_at_eta(e);
  }
  
  return background_expansion;
  
}


/*******************************************
 *******************************************
 **________ 4. OUTPUT AND CHECKS _________**
 *******************************************
 ******************************************* 
 *                                         *
 * ..... 4.1 print_background_cosmology    *
 *                                         *
 *******************************************
 *******************************************/

/*******************************************************************************************************************************************************
 * 4.1 print_background_cosmology
 * Description:
 *  - Prints table with expansion history of the universe.
 * Arguments:
 *  - file_name: file name to which the table shall be output
 * 
*******************************************************************************************************************************************************/

void Universe::print_background_cosmology(string file_name){
  
  double Om_m = this->cosmology.Omega_m;
  double Om_l = this->cosmology.Omega_L;
  double Om_r = this->cosmology.Omega_r;
  fstream cosmo_stream;
  
  remove(file_name.c_str());
  FILE * F = fopen(file_name.c_str(), "w");
  fclose(F);
  cosmo_stream.open(file_name.c_str());
  
  cosmo_stream << setprecision(10) << scientific;
  cosmo_stream << "#Cosmological Parameters: Om_m = " << Om_m << ", Om_l = " << Om_l << ", Om_r = " << Om_r << '\n';
  cosmo_stream << "#a_max = " << this->a_final << '\n';
  cosmo_stream << "#t" << setw(20) << "eta(t)" << setw(20) << "w(eta)" << setw(20) << "z(eta)" << setw(20) << "a(eta)" << setw(20) << "H(eta)" << setw(20) << "H_prime(eta)\n";
  for(int i = 0; i < this->a.size(); i++){
    cosmo_stream << this->t[i];
    cosmo_stream << setw(20) << this->eta[i];
    cosmo_stream << setw(20) << this->eta_at_a(1.0)-this->eta[i];
    cosmo_stream << setw(20) << 1.0/this->a[i] - 1.0;
    cosmo_stream << setw(20) << this->a[i];
    cosmo_stream << setw(20) << this->H[i];
    cosmo_stream << setw(20) << this->H_prime[i] << '\n';
  }
    
  cosmo_stream.close();
  
}




  




