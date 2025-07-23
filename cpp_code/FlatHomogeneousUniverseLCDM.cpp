#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "FlatHomogeneousUniverseLCDM_integrands.cpp"
#include "FlatHomogeneousUniverseLCDM_analytic_special_cases.cpp"
#include "FlatHomogeneousUniverseLCDM_initial_conditions.cpp"
 
using namespace std;
using namespace constants;


/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::FlatHomogeneousUniverseLCDM
 * 
 * Constructor of class FlatHomogeneousUniverseLCDM
 * 
*******************************************************************************************************************************************************/


FlatHomogeneousUniverseLCDM::FlatHomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max) : FriedmannUniverse(){
  
  { using namespace error_handling;
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
    
    error_message = "Omega_m + Omega_r + Omega_L = " + to_string(cosmo.Omega_m+cosmo.Omega_r+cosmo.Omega_L)+" should be equal to 1 in flat universe .";
    test_value_of_double(constants::eps_from_0, abs(cosmo.Omega_m+cosmo.Omega_r+cosmo.Omega_L-1.0), LARGER_EQUAL, error_message);
    
    error_message = "h_100 = " + to_string(cosmo.h_100)+" should be > 0.0 .";
    test_value_of_double(cosmo.h_100, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "sigma_8 = " + to_string(cosmo.sigma_8)+" should be > 0.0 .";
    test_value_of_double(cosmo.sigma_8, 0.0, LARGER_EQUAL, error_message);
    
    error_message = "a_max = " + to_string(a_max)+" should be > a_min = " + to_string(a_min)+" .";
    test_value_of_double(a_max, a_min, LARGER_EQUAL, error_message);
  }

  this->a_initial = a_min;
  this->a_final = a_max;
  this->set_spatial_curvature(0.0);
  this->cosmology = cosmo;
  this->set_initial_conditions();
  this->set_expansion_history();

}


/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::~FlatHomogeneousUniverseLCDM
 * 
 * Destructor of class FlatHomogeneousUniverseLCDM
 * 
*******************************************************************************************************************************************************/

FlatHomogeneousUniverseLCDM::~FlatHomogeneousUniverseLCDM(){
}

/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::set_expansion_history
 * 
 * computing the expansion history of the universe 
 * 
*******************************************************************************************************************************************************/

void FlatHomogeneousUniverseLCDM::set_expansion_history(){
  
  integration_parameters_Universe params;
  params.pointer_to_Universe = this;
  integration_parameters_Universe * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {scale_factor_gsl, scale_factor_gsl_jac, 3, (void *) pointer_to_params};

  
  double a_i = this->a_initial;
  double a_f;
  double eta_i, d_eta;
  
  vector<double> a_dummy(1, 0.0);
  vector<double> eta_dummy(1, 0.0);
  vector<double> t_dummy(1, 0.0);
  double y[3];
  a_dummy[0] = this->a_initial;
  eta_dummy[0] = this->eta_initial;
  t_dummy[0] = this->t_initial;
  
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
  int status;
  double hstart = 1.0*constants::gsl_hstart_relative;
  double eps_absolute = constants::gsl_eps_relative*max(1.0, this->H_initial);
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, eps_absolute, constants::gsl_eps_relative);

  double da_attempt;
  int step = 0;
  while(a_i < this->a_final){
    step++;
    da_attempt = constants::maximal_dw*this->hubble_from_Friedmann(a_dummy[step-1])*a_dummy[step-1];
    da_attempt = min(constants::maximal_da, da_attempt);
    da_attempt = min(this->a_final - a_dummy[step-1], da_attempt);
    do{
      a_i = a_dummy[step-1];
      a_f = a_i + da_attempt;
      y[0] = eta_dummy[step-1];
      y[1] = this->hubble_from_Friedmann(a_i);
      y[2] = t_dummy[step-1];
      status = gsl_odeiv2_driver_apply(d, &a_i, a_f, y);
      d_eta = y[0]-eta_dummy[step-1];
      da_attempt *= 0.9;
    }while(d_eta > constants::maximal_dw);
    
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }

    eta_dummy.push_back(y[0]);
    t_dummy.push_back(y[2]);
    a_dummy.push_back(a_f);
  }
  
  gsl_odeiv2_driver_free(d);
  
  int number_of_time_steps = a_dummy.size();
  
  this->a = a_dummy;
  this->eta = eta_dummy;
  this->t = t_dummy;
  this->H.resize(number_of_time_steps);
  this->H_prime.resize(number_of_time_steps);
  
  for (int i = 0; i < number_of_time_steps; i++){
    this->H[i] = this->hubble_from_Friedmann(this->a[i]);
    this->H_prime[i] = this->hubble_prime_from_Friedmann(this->a[i]);
  }
}

/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::hubble_from_Friedmann
 * Gives conformal expansion rate H at given scale factor. Within the code, this function is mainly used
 * to set the initial value of H when integrating the Friedmann equations.
 * 
*******************************************************************************************************************************************************/

 double FlatHomogeneousUniverseLCDM::hubble_from_Friedmann(double a_start){
  
  double rho_m = rho_m_of_a(a_start);
  double rho_l = rho_L_of_a(a_start);
  double rho_r = rho_r_of_a(a_start);

  return sqrt(a_start*a_start*(rho_m+rho_r+rho_l)+this->cosmology.Omega_k);
  
}

/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::hubble_prime_from_Friedmann
 * Gives derivative of conformal expansion rate, H', at given scale factor. Within the code, this function
 * is mainly used to set the initial value of H' when integrating the Friedmann equations.
 * 
*******************************************************************************************************************************************************/

double FlatHomogeneousUniverseLCDM::hubble_prime_from_Friedmann(double scale){
  
  double rho_m = rho_m_of_a(scale);
  double rho_l = rho_L_of_a(scale);
  double rho_r = rho_r_of_a(scale);
  double w = this->w_L_of_a(scale);

  return -0.5*(rho_m + 2*rho_r+(1.0+3.0*w)*rho_l)*scale*scale;
  
}


/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::return_cosmology
 * Description:
 *  - Returns struct with cosmological parameters.
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

cosmological_model FlatHomogeneousUniverseLCDM::return_cosmology(){
  
  return this->cosmology;
  
}

/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::rho_m_of_a
 * Description:
 *  - Matter density at scale factor a in units of present day critical density.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/

double FlatHomogeneousUniverseLCDM::rho_m_of_a(double a){
  
  return this->cosmology.Omega_m/pow(a, 3.0);
  
}



/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::rho_r_of_a
 * Description:
 *  - Radiation density at scale factor a in units of present day critical density.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/


double FlatHomogeneousUniverseLCDM::rho_r_of_a(double a){
  
  return this->cosmology.Omega_r/pow(a, 4.0);
  
}


/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::rho_L_of_a
 * Description:
 *  - Dark energy density at scale factor a in units of present day critical density.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/


double FlatHomogeneousUniverseLCDM::rho_L_of_a(double a){
  
  double w = this->w_L_of_a(a);
  
  return this->cosmology.Omega_L/pow(a, 3.0*(1.0+w));
  
}


/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::w_L_of_a
 * Description:
 * - Equation of State parameter of Dark energy. 
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/

double FlatHomogeneousUniverseLCDM::w_L_of_a(double a){
  
  double w1 = this->cosmology.w1;
  if(w1 == 0.0)
    return this->cosmology.w0;
  
  return this->cosmology.w0 + this->cosmology.w1*(1.0-a); // Chevallier-Polarski-Linde
  //return this->cosmology.w0 + this->cosmology.w1*(1.0-a)/a; // Linear-redshift
  //return this->cosmology.w0 + this->cosmology.w1*(1.0-a)*a; // Jassal-Bagla-Padmanabhan
  
}
  
  
/*******************************************************************************************************************************************************
 * FlatHomogeneousUniverseLCDM::Mass_within_R_in_Mpc_over_h
 * Description:
 * - Mass within a given radius according to critical density [in M_solar/h]
 * Arguments:
 *  - R_in_Mpc_over_h: radius
 * 
*******************************************************************************************************************************************************/

double FlatHomogeneousUniverseLCDM::Mass_within_R_in_Mpc_over_h(double R_in_Mpc_over_h){
  double R = R_in_Mpc_over_h/constants::c_over_e5; //changing to units of c/H_0
  double M0_crit = constants::M0_crit_in_Hubble_radius_in_Msol_over_h*pow(R,3);
  double M0 = M0_crit*this->cosmology.Omega_m;
  
  return M0;
  
}


double FlatHomogeneousUniverseLCDM::R_in_Mpc_over_h_enclosing_M(double M_in_Msol_over_h){
  double R = pow(M_in_Msol_over_h/(constants::M0_crit_in_Hubble_radius_in_Msol_over_h*this->cosmology.Omega_m), 1.0/3.0);
  double R_in_Mpc_over_h = R*constants::c_over_e5;
  
  return R_in_Mpc_over_h;
}




