 
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

struct integration_parameters_Universe{
  
  cosmological_model cosmo;
  
  Universe* pointer_to_Universe;
  
};

int scale_factor_gsl(double scale, const double y[], double dfda[], void *params){
  
  integration_parameters_Universe *integration_params = (integration_parameters_Universe *) params;
  Universe* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double hubble = y[1];
  
  
  double rho_m = pointer_to_Universe->rho_m_of_a(scale);
  double rho_l = pointer_to_Universe->rho_L_of_a(scale);
  double rho_r = pointer_to_Universe->rho_r_of_a(scale);
  double w = pointer_to_Universe->w_L_of_a(scale);

  dfda[0] = 1.0/(scale*hubble); // derivative: conformal time as a function of scale factor
  dfda[1] = -0.5*(rho_m + 2*rho_r+(1.0+3.0*w)*rho_l)*scale/hubble; // derivative: conformal expansion rate as a function of scale factor
  dfda[2] = 1.0/hubble; // derivative: physical time as a function of scale factor
  return GSL_SUCCESS;
}


int scale_factor_gsl_jac(double scale, const double y[], double *dfdy, double dfde[], void *params){

  integration_parameters_Universe *integration_params = (integration_parameters_Universe *) params;
  Universe* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double hubble = y[1];
  
  double rho_m = pointer_to_Universe->rho_m_of_a(scale);
  double rho_l = pointer_to_Universe->rho_L_of_a(scale);
  double rho_r = pointer_to_Universe->rho_r_of_a(scale);
  double w = pointer_to_Universe->w_L_of_a(scale);
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix;
  
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, -1.0/(scale*hubble*hubble));
  gsl_matrix_set (m, 0, 2, 0.0);
  
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, 0.5*(rho_m + 2*rho_r+(1.0+3.0*w)*rho_l)*scale/pow(hubble, 2));
  gsl_matrix_set (m, 1, 2, 0.0);
  
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, -1.0/(hubble*hubble));
  gsl_matrix_set (m, 2, 2, 0.0);

  return GSL_SUCCESS;
}




void Universe::set_number_of_time_steps(int n_entries){

  /*
  
  Planck 1-sigma uncertainties:
  Omega_m: 0.013
  Omega_b: 0.0005
  sigma_8: 0.014
      n_s: 0.0062
    h_100: 0.0096
      
  3-sigma uncertainties on Omega_b and h_100 by Olive et al. (2014) and Riess et al. (2016):
  
  Omega_b: 0.0042
    h_100: 0.054
  
  Fisher forecast of 1-sigma uncertainties for non-tomographic 20arcmin trough lensing:
  Omega_m: 0.0690696011991
  sigma_8: 0.0725324281146
      n_s: 0.14309020171
  
  */
  
  this->number_of_time_steps = n_entries;
  
  
  double dummy;
  /*
  dummy = 3.0*sin(10000.0);
  this->cosmology.Omega_m += dummy*0.0690696011991;
  this->cosmology.Omega_L = 1.0 - this->cosmology.Omega_m;
  dummy = 3.0*sin(20000.0);
  this->cosmology.Omega_b += dummy*0.0042;
  dummy = 3.0*sin(30000.0);
  this->cosmology.n_s += dummy*0.14309020171;
  dummy = 3.0*sin(40000.0);
  this->cosmology.h_100 += dummy*0.054;
  dummy = 3.0*sin(50000.0);
  this->cosmology.sigma_8 += dummy*0.0725324281146;
  */
  
}
