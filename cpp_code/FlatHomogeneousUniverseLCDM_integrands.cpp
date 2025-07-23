 
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

struct integration_parameters_Universe{
  
  cosmological_model cosmo;
  
  FlatHomogeneousUniverseLCDM* pointer_to_Universe;
  
};

int scale_factor_gsl(double scale, const double y[], double dfda[], void *params){
  
  integration_parameters_Universe *integration_params = (integration_parameters_Universe *) params;
  FlatHomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
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
  FlatHomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
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
