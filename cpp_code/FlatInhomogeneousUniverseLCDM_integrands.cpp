#include <stdio.h>
#include <gsl/gsl_randist.h>


struct integration_parameters{
  
  double top_hat_radius;
  double top_hat_length;
  double second_top_hat_radius;
  double k;
  double k_1, k_2;
  double lnk;
  double r;
  double n_s;
  double Omega_m;
  
  double WR_1, WR_2, WR_3;
  double dWR_1, dWR_2, dWR_3;
  double d2WR_1, d2WR_2, d2WR_3;
  
  double alpha_1, alpha_2, alpha_3;
  
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe;
  
};





double int_gsl_integrate_high_precision(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  int workspace_size = 10000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = func;
  F.params  = arg;
  int key = 6;
  gsl_integration_qag(&F,a,b,0, 1e-7, workspace_size, key, w, &res, &err);
  if(NULL!=error)
    *error=err;
  gsl_integration_workspace_free(w);
  return res;
}

double int_gsl_integrate_medium_precision(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
  gsl_integration_cquad(&F,a,b,0,1.0e-6,w,&res,&err,0);
  if(NULL!=error)
    *error=err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}


double int_gsl_integrate_low_precision(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
  //gsl_integration_cquad(&F,a,b,0,1.0e-4,w,&res,&err,0);
  gsl_integration_cquad(&F,a,b,0,1.0e-3,w,&res,&err,0);
  // --> ISSUE: precision of 1.0e-3 may not be good enough in skewness integral for all PNG types.
  if(NULL!=error)
    *error=err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}


double norm_derivs_before_norm_was_determined_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 3.0+integration_params->n_s;
  double WR = w_R(k, integration_params->top_hat_radius);
  double T_sq = pointer_to_Universe->transfer_function_at(k); T_sq *= T_sq;
  
  return one_over_2_pi_sq*pow(k,index)*T_sq*WR*WR;
  
}


double norm_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 3.0;
  double WR = w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return one_over_2_pi_sq*pow(k,index)*P*WR*WR;
  
}

double norm_NL_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 3.0;
  double WR = w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_NL_at(lnk);
  
  return one_over_2_pi_sq*pow(k,index)*P*WR*WR;
  
}

double dvar_dR_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 3.0;
  double WR = w_R(k, integration_params->top_hat_radius);
  double dWR = deriv_of_w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return 2.0*one_over_2_pi_sq*pow(k,index)*P*dWR*WR;
  
}



double var_derivs_2D_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 2.0;
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k,index)*P*WR*WR;
  
}

double cov_derivs_2D_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 2.0;
  double WR1 = w_R_2D(k, integration_params->top_hat_radius);
  double WR2 = w_R_2D(k, integration_params->second_top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k,index)*P*WR1*WR2;
  
}

double dcov_dR2_derivs_2D_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 2.0;
  double WR1 = w_R_2D(k, integration_params->top_hat_radius);
  double dWR2_dR2 = deriv_of_w_R_2D(k, integration_params->second_top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k,index)*P*WR1*dWR2_dR2;
  
}

double squared_saddle_point_derivs_2D_gsl(double R, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double R_lin = integration_params->top_hat_radius;
  integration_params->second_top_hat_radius = R;
  double Cov = int_gsl_integrate_medium_precision(cov_derivs_2D_gsl,params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  double dCov_dR = int_gsl_integrate_medium_precision(dcov_dR2_derivs_2D_gsl,params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  
  return R*pow(Cov + 0.5*R*dCov_dR, 2);
  
}

double var_NL_derivs_2D_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 2.0;
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_NL_at(lnk);
  
  return pow(k,index)*P*WR*WR;
  
}

double var_NL_derivs_2D_2_gsl(double lnk_2D, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k_2D = exp(lnk_2D);
  double lnk = 0.5*log(k_2D*k_2D+pow(integration_params->k,2));
  double index = 2.0;
  double WR = w_R_2D(k_2D, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_NL_at(lnk);
  
  return pow(k_2D,index)*P*WR*WR;
  
}

double var_NL_derivs_2D_1_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 1.0;
  double WL = w_L_1D(k, integration_params->top_hat_length);
  integration_params->k = k;
  
  return pow(k,index)*WL*WL*int_gsl_integrate_medium_precision(var_NL_derivs_2D_2_gsl,params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  
}

double var_derivs_2D_2_gsl(double lnk_2D, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k_2D = exp(lnk_2D);
  double lnk = 0.5*log(k_2D*k_2D+pow(integration_params->k,2));
  double index = 2.0;
  double WR = w_R_2D(k_2D, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k_2D,index)*P*WR*WR;
  
}

double var_derivs_2D_1_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 1.0;
  double WL = w_L_1D(k, integration_params->top_hat_length);
  integration_params->k = k;
  
  return pow(k,index)*WL*WL*int_gsl_integrate_medium_precision(var_derivs_2D_2_gsl,params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);;
  
}

double dvar_dR_derivs_2D_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 2.0;
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double dWR_dR = deriv_of_w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k,index)*P*2.0*dWR_dR*WR;
  
}

double d2var_dR2_derivs_2D_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 2.0;
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double dWR_dR = deriv_of_w_R_2D(k, integration_params->top_hat_radius);
  double d2WR_dR2 = second_deriv_of_w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k,index)*P*2.0*(d2WR_dR2*WR + dWR_dR*dWR_dR);
  
}

double dvar_dR_derivs_2D_2_gsl(double lnk_2D, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k_2D = exp(lnk_2D);
  double lnk = 0.5*log(k_2D*k_2D+pow(integration_params->k,2));
  double index = 2.0;
  double WR = w_R_2D(k_2D, integration_params->top_hat_radius);
  double dWR_dR = deriv_of_w_R_2D(k_2D, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k_2D,index)*P*2.0*WR*dWR_dR;
  
}

double dvar_dR_derivs_2D_1_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 1.0;
  double WL = w_L_1D(k, integration_params->top_hat_length);
  integration_params->k = k;
  
  return pow(k,index)*WL*WL*int_gsl_integrate_medium_precision(dvar_dR_derivs_2D_2_gsl,params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  
}

double P_to_P2D_derivs_gsl(double lnk_parallel, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k_parallel = exp(lnk_parallel);
  double lnk = 0.5*log(k_parallel*k_parallel+pow(integration_params->k,2));
  double index = 1.0;
  double WL = w_L_1D(k_parallel, integration_params->top_hat_length);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return pow(k_parallel,index)*WL*WL*P;
}


double skewness_integral_3_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double WR = w_R(k, integration_params->top_hat_radius);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double alpha = integration_params->alpha_3;
  
  return pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*WR*pow(P, alpha);
  
}

double skewness_integral_2_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double WR = w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double alpha = integration_params->alpha_2;
  
  double q_min;
  if(k>integration_params->k){
    q_min = max(k-integration_params->k, minimal_wave_number_in_H0_units);
  }
  else{
    q_min = max(integration_params->k-k, minimal_wave_number_in_H0_units);
  }
  double q_max = min(k+integration_params->k, maximal_wave_number_in_H0_units);
  
  double integral = int_gsl_integrate_low_precision(skewness_integral_3_derivs_gsl,params,log(q_min),log(q_max),NULL,1000);
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*WR*pow(P, alpha);
  
}
 

double skewness_integral_1_derivs_gsl(double lnk, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->k = k;
  double WR = w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P_phi = P;
  
  double log_k_max = min(constants::product_of_kmax_and_R/integration_params->top_hat_radius, maximal_wave_number_in_H0_units);
  log_k_max = log(log_k_max);
  
  double integral = int_gsl_integrate_low_precision(skewness_integral_2_derivs_gsl,params,log_minimal_wave_number_in_H0_units,log_k_max,NULL,1000);
  
  double result = 0.0;
  double alpha = integration_params->alpha_1;
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*WR*pow(P, alpha);
  
}




double dskewness_integral_3_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->WR_3 = w_R(k, integration_params->top_hat_radius);
  integration_params->dWR_3 = deriv_of_w_R(k, integration_params->top_hat_radius);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double alpha = integration_params->alpha_3;
  
  double W_dervis = 0.0;
  W_dervis += integration_params->dWR_1*integration_params->WR_2*integration_params->WR_3;
  W_dervis += integration_params->WR_1*integration_params->dWR_2*integration_params->WR_3;
  W_dervis += integration_params->WR_1*integration_params->WR_2*integration_params->dWR_3;
  
  return pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*pow(P, alpha)*W_dervis;
  
}

double dskewness_integral_2_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->WR_2 = w_R(k, integration_params->top_hat_radius);
  integration_params->dWR_2 = deriv_of_w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double alpha = integration_params->alpha_2;
  
  double q_min;
  if(k>integration_params->k){
    q_min = max(k-integration_params->k, minimal_wave_number_in_H0_units);
  }
  else{
    q_min = max(integration_params->k-k, minimal_wave_number_in_H0_units);
  }
  double q_max = min(k+integration_params->k, maximal_wave_number_in_H0_units);
  
  double integral = int_gsl_integrate_low_precision(dskewness_integral_3_derivs_gsl,params,log(q_min),log(q_max),NULL,1000);
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*pow(P, alpha);
  
}
 

double dskewness_integral_1_derivs_gsl(double lnk, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->k = k;
  integration_params->WR_1 = w_R(k, integration_params->top_hat_radius);
  integration_params->dWR_1 = deriv_of_w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P_phi = P;
  
  double log_k_max = min(constants::product_of_kmax_and_R/integration_params->top_hat_radius, maximal_wave_number_in_H0_units);
  log_k_max = log(log_k_max);
  
  double integral = int_gsl_integrate_low_precision(dskewness_integral_2_derivs_gsl,params,log_minimal_wave_number_in_H0_units,log_k_max,NULL,1000);
  
  double result = 0.0;
  double alpha = integration_params->alpha_1;
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*pow(P, alpha);
  
}


/**************************************
 * 
 * 2D skewness integrands start here.
 * 
 * ************************************/




double skewness_integral_2D_3_derivs_gsl(double phi, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k_1 = integration_params->k_1;
  double k_2 = integration_params->k_2;
  double cosphi = cos(phi);
  double k;
  
  if (cosphi >= 1.0)
    k = k_1 + k_2;
  else if (cosphi <= -1.0)
    k = abs(k_1 - k_2);
  else
    k = sqrt(k_1*k_1 + k_2*k_2 + 2.0*k_1*k_2*cosphi);
  
  if(k > maximal_wave_number_in_H0_units || k < minimal_wave_number_in_H0_units)
    return 0.0;
  
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P = pointer_to_Universe->current_P_L_at(log(k));
  double alpha = integration_params->alpha_3;
  
  return pow(k, 2.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*WR*pow(P, alpha);
  
}

double skewness_integral_2D_2_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->k_2 = k;
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double alpha = integration_params->alpha_2;
  
  double integral = int_gsl_integrate_low_precision(skewness_integral_2D_3_derivs_gsl,params,0.0,constants::pi,NULL,1000);
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*WR*pow(P, alpha);
  
}
 

double skewness_integral_2D_1_derivs_gsl(double lnk, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->k_1 = k;
  double WR = w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P_phi = P;
  
  double log_k_max = min(constants::product_of_kmax_and_R/integration_params->top_hat_radius, maximal_wave_number_in_H0_units);
  log_k_max = log(log_k_max);
  
  double integral = int_gsl_integrate_low_precision(skewness_integral_2D_2_derivs_gsl,params,log_minimal_wave_number_in_H0_units,log_k_max,NULL,1000);
  
  double result = 0.0;
  double alpha = integration_params->alpha_1;
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*WR*pow(P, alpha);
  
}




double dskewness_integral_2D_3_derivs_gsl(double phi, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k_1 = integration_params->k_1;
  double k_2 = integration_params->k_2;
  double cosphi = cos(phi);
  double k;
  
  if (cosphi >= 1.0)
    k = k_1 + k_2;
  else if (cosphi <= -1.0)
    k = abs(k_1 - k_2);
  else
    k = sqrt(k_1*k_1 + k_2*k_2 + 2.0*k_1*k_2*cosphi);
  
  if(k > maximal_wave_number_in_H0_units || k < minimal_wave_number_in_H0_units)
    return 0.0;
  
  
  integration_params->WR_3 = w_R_2D(k, integration_params->top_hat_radius);
  integration_params->dWR_3 = deriv_of_w_R_2D(k, integration_params->top_hat_radius);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P = pointer_to_Universe->current_P_L_at(log(k));
  double alpha = integration_params->alpha_3;
  
  double W_dervis = 0.0;
  W_dervis += integration_params->dWR_1*integration_params->WR_2*integration_params->WR_3;
  W_dervis += integration_params->WR_1*integration_params->dWR_2*integration_params->WR_3;
  W_dervis += integration_params->WR_1*integration_params->WR_2*integration_params->dWR_3;
  
  return pow(k, 2.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*W_dervis*pow(P, alpha);
  
}

double dskewness_integral_2D_2_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->k_2 = k;
  integration_params->WR_2 = w_R_2D(k, integration_params->top_hat_radius);
  integration_params->dWR_2 = deriv_of_w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double alpha = integration_params->alpha_2;
  
  double integral = int_gsl_integrate_low_precision(dskewness_integral_2D_3_derivs_gsl,params,0.0,constants::pi,NULL,1000);
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*pow(P, alpha);
  
}
 

double dskewness_integral_2D_1_derivs_gsl(double lnk, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  integration_params->k_1 = k;
  integration_params->WR_1 = w_R_2D(k, integration_params->top_hat_radius);
  integration_params->dWR_1 = deriv_of_w_R_2D(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  double T = pointer_to_Universe->transfer_function_at(k);
  double P_phi = P;
  
  double log_k_max = min(constants::product_of_kmax_and_R/integration_params->top_hat_radius, maximal_wave_number_in_H0_units);
  log_k_max = log(log_k_max);
  
  double integral = int_gsl_integrate_low_precision(dskewness_integral_2D_2_derivs_gsl,params,log_minimal_wave_number_in_H0_units,log_k_max,NULL,1000);
  
  double result = 0.0;
  double alpha = integration_params->alpha_1;
  
  return integral*pow(k, 4.0-4.0*alpha)*pow(T, 1.0-2.0*alpha)*pow(P, alpha);
  
}




/**************************************
 * 
 * 2D skewness integrands end here.
 * 
 * ************************************/

/*
double norm_derivs_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double k = exp(lnk);
  double index = 3.0;
  double WR = w_R(k, integration_params->top_hat_radius);
  double P = pointer_to_Universe->current_P_L_at(lnk);
  
  return one_over_2_pi_sq*pow(k,index)*P*WR*WR;
  
}
*/


double halofit_sig_sq_gsl(double lnk, void *params){
 
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double index = 3.0+integration_params->n_s;
  
  double k = exp(lnk);
  double y_sq = pow(integration_params->top_hat_radius*k, 2);
  double T_sq = pointer_to_Universe->transfer_function_at(k); T_sq *= T_sq;
  

  return pow(k,index)*T_sq*exp(-y_sq);
  
}


double halofit_C_gsl(double lnk, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double index = 3.0+integration_params->n_s;
  double k = exp(lnk);
  double k_sq = k*k;
  double y_sq = pow(integration_params->top_hat_radius*k, 2);
  double T_sq = pointer_to_Universe->transfer_function_at(k); T_sq *= T_sq;
  double integrand = pow(k,index)*T_sq*exp(-y_sq)*y_sq;
    
  return integrand;
  
}


double halofit_n_gsl(double lnk, void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double index = 3.0+integration_params->n_s;
  double k = exp(lnk);
  double k_sq = k*k;
  double y_sq = pow(integration_params->top_hat_radius*k, 2);
  double T_sq = pointer_to_Universe->transfer_function_at(k); T_sq *= T_sq;
  double integrand = pow(k,index)*T_sq*exp(-y_sq)*y_sq;
    
  return (1-y_sq)*integrand;
  
}

int growth_factor_gsl(double e, const double D[], double dDde[], void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double Om_m = integration_params->Omega_m;
  
  double D0 = D[0];
  double D1 = D[1];

  double scale = pointer_to_Universe->a_at_eta(e);
  double H  = pointer_to_Universe->H_at_eta(e);

  dDde[0] = D1;
  dDde[1] = -H*D1+3.0/2.0*Om_m/scale*D0;
  
  return GSL_SUCCESS;
}

int growth_factor_to_second_order_gsl(double e, const double D[], double dDde[], void *params){
  
  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double Om_m = integration_params->Omega_m;
  
  double D0 = D[0];
  double D1 = D[1];

  double scale = pointer_to_Universe->a_at_eta(e);
  double H = pointer_to_Universe->H_at_eta(e);
  double D_lin_prime = pointer_to_Universe->return_D_prime_of_eta(e);

  dDde[0] = D1;
  dDde[1] = -H*D1+3.0/2.0*Om_m/scale*D0 + D_lin_prime*D_lin_prime;
  
  return GSL_SUCCESS;
  
}


int growth_factor_gsl_jac(double e, const double D[], double *dfdD, double dDde[], void *params){

  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double Om_m = integration_params->Omega_m;
  
  double D0 = D[0];
  double D1 = D[1];

  double scale = pointer_to_Universe->a_at_eta(e);
  double H = pointer_to_Universe->H_at_eta(e);

  dDde[0] = D1;
  dDde[1] = -H*D1+3.0/2.0*Om_m/scale*D0;
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdD, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  
  gsl_matrix_set (m, 1, 0, 3.0/2.0*Om_m/scale);
  gsl_matrix_set (m, 1, 1, -H);

  return GSL_SUCCESS;
}


int growth_factor_to_second_order_gsl_jac(double e, const double D[], double *dfdD, double dDde[], void *params){

  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double Om_m = integration_params->Omega_m;
  
  double D0 = D[0];
  double D1 = D[1];

  double scale = pointer_to_Universe->a_at_eta(e);
  double H = pointer_to_Universe->H_at_eta(e);

  dDde[0] = D1;
  dDde[1] = -H*D1+3.0/2.0*Om_m/scale*D0;
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdD, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  
  gsl_matrix_set (m, 1, 0, 3.0/2.0*Om_m/scale);
  gsl_matrix_set (m, 1, 1, -H);

  return GSL_SUCCESS;
}





/*
 * 
 * Spherical collapse
 * 
 * 
 * 
 */


int F_dF_ddF_spherical_wrt_delta_gsl(double e, const double y[], double dfde[], void *params){

  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double F = y[0];
  double F_prime = y[1];
  double dF_ddelta = y[2];
  double dF_ddelta_prime = y[3];
  double d2F_ddelta2 = y[4];
  double d2F_ddelta2_prime = y[5];
  double a = pointer_to_Universe->a_at_eta(e);
  double hubble = pointer_to_Universe->H_at_eta(e);
  double Omega_m = integration_params->Omega_m;

  dfde[0] = F_prime;
  dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 4.0/3.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  dfde[2] = dF_ddelta_prime;
  dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 8.0/3.0*F_prime*dF_ddelta_prime/(1+F) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  dfde[4] = d2F_ddelta2_prime;
  dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 8.0/3.0*pow(dF_ddelta_prime, 2)/(1+F) + 8.0/3.0*F_prime*d2F_ddelta2_prime/(1+F) - 16.0/3.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 8.0/3.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;

  
  return GSL_SUCCESS;
}

int F_dF_ddF_spherical_wrt_delta_jac(double e, const double y[], double *dfdy, double dfde[], void *params){

  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;

  double F = y[0];
  double F_prime = y[1];
  double dF_ddelta = y[2];
  double dF_ddelta_prime = y[3];
  double d2F_ddelta2 = y[4];
  double d2F_ddelta2_prime = y[5];
  double a = pointer_to_Universe->a_at_eta(e);
  double hubble = pointer_to_Universe->H_at_eta(e);
  double Omega_m = integration_params->Omega_m;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
  gsl_matrix * m = &dfdy_mat.matrix; 
  
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 0, 4, 0.0);
  gsl_matrix_set (m, 0, 5, 0.0);
  
  //dfde[0] = F_prime;
  //dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 4.0/3.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  //dfde[2] = dF_ddelta_prime;
  //dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 8.0/3.0*F_prime*dF_ddelta_prime/(1+F) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  //dfde[4] = d2F_ddelta2_prime;
  //dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 8.0/3.0*pow(dF_ddelta_prime, 2)/(1+F) + 8.0/3.0*F_prime*d2F_ddelta2_prime/(1+F) - 16.0/3.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 8.0/3.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;
  
  gsl_matrix_set (m, 1, 0, 1.5*Omega_m/a*(1.0+2.0*F) - 4.0/3.0*pow(F_prime/(1.0+F), 2.0));
  gsl_matrix_set (m, 1, 1, 8.0/3.0*F_prime/(1+F) - hubble);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 1, 4, 0.0);
  gsl_matrix_set (m, 1, 5, 0.0);
  
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 3, 1.0);
  gsl_matrix_set (m, 2, 4, 0.0);
  gsl_matrix_set (m, 2, 5, 0.0);
  
  //dfde[0] = F_prime;
  //dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 4.0/3.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  //dfde[2] = dF_ddelta_prime;
  //dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 8.0/3.0*F_prime*dF_ddelta_prime/(1+F) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  //dfde[4] = d2F_ddelta2_prime;
  //dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 8.0/3.0*pow(dF_ddelta_prime, 2)/(1+F) + 8.0/3.0*F_prime*d2F_ddelta2_prime/(1+F) - 16.0/3.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 8.0/3.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;
  
  gsl_matrix_set (m, 3, 0, 3.0*Omega_m/a*dF_ddelta - 8.0/3.0*F_prime*dF_ddelta_prime/pow(1.0+F, 2.0) + 8.0/3.0*pow(F_prime, 2)/pow(1.0+F, 3)*dF_ddelta);
  gsl_matrix_set (m, 3, 1, 8.0/3.0*dF_ddelta_prime/(1+F) - 8.0/3.0*F_prime/pow(1+F, 2)*dF_ddelta);
  gsl_matrix_set (m, 3, 2, 1.5*Omega_m/a*(1.0+2.0*F) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2));
  gsl_matrix_set (m, 3, 3, 8.0/3.0*F_prime/(1+F) - hubble);
  gsl_matrix_set (m, 3, 4, 0.0);
  gsl_matrix_set (m, 3, 5, 0.0);
  
  gsl_matrix_set (m, 4, 0, 0.0);
  gsl_matrix_set (m, 4, 1, 0.0);
  gsl_matrix_set (m, 4, 2, 0.0);
  gsl_matrix_set (m, 4, 3, 0.0);
  gsl_matrix_set (m, 4, 4, 0.0);
  gsl_matrix_set (m, 4, 5, 1.0);
  
  //dfde[0] = F_prime;
  //dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 4.0/3.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  //dfde[2] = dF_ddelta_prime;
  //dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 8.0/3.0*F_prime*dF_ddelta_prime/(1+F) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  //dfde[4] = d2F_ddelta2_prime;
  //dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 8.0/3.0*pow(dF_ddelta_prime, 2)/(1+F) + 8.0/3.0*F_prime*d2F_ddelta2_prime/(1+F) - 16.0/3.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 8.0/3.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;
  
  gsl_matrix_set (m, 5, 0, 3.0*Omega_m/a*d2F_ddelta2 - 8.0/3.0*pow(dF_ddelta_prime, 2)/pow(1.0+F,2.0) - 8.0/3.0*F_prime*d2F_ddelta2_prime/pow(1.0+F,2.0) + 32.0/3.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 3.0) + 8.0/3.0*pow(F_prime, 2)/pow(1+F, 3)*d2F_ddelta2 - 24.0/3.0*pow(F_prime, 2)/pow(1+F, 4)*pow(dF_ddelta, 2.0));
  gsl_matrix_set (m, 5, 1, 8.0/3.0*d2F_ddelta2_prime/(1+F) - 16.0/3.0*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 8.0/3.0*F_prime/pow(1+F, 2)*d2F_ddelta2 + 16.0/3.0*F_prime/pow(1+F, 3)*pow(dF_ddelta, 2.0));
  gsl_matrix_set (m, 5, 2, 6.0*Omega_m/a*dF_ddelta - 16.0/3.0*F_prime*dF_ddelta_prime/pow(1+F, 2.0) + 16.0/3.0*pow(F_prime, 2)/pow(1+F, 3)*dF_ddelta);
  gsl_matrix_set (m, 5, 3, 16.0/3.0*dF_ddelta_prime/(1+F) - 16.0/3.0*F_prime*dF_ddelta/pow(1+F, 2.0));
  gsl_matrix_set (m, 5, 4, 1.5*Omega_m/a*(1.0+2.0*F) - 4.0/3.0*pow(F_prime, 2)/pow(1+F, 2));
  gsl_matrix_set (m, 5, 5, 8.0/3.0*F_prime/(1+F) - hubble);
  
  return GSL_SUCCESS;
}


// <--- You got until here with checking the equations!






/*
 * 
 * 2D collapse
 * 
 * 
 * 
 */


int F_dF_ddF_cylindrical_wrt_delta_gsl(double e, const double y[], double dfde[], void *params){

  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;
  
  double F = y[0];
  double F_prime = y[1];
  double dF_ddelta = y[2];
  double dF_ddelta_prime = y[3];
  double d2F_ddelta2 = y[4];
  double d2F_ddelta2_prime = y[5];
  double a = pointer_to_Universe->a_at_eta(e);
  double hubble = pointer_to_Universe->H_at_eta(e);
  double Omega_m = integration_params->Omega_m;

  dfde[0] = F_prime;
  dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 3.0/2.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  dfde[2] = dF_ddelta_prime;
  dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 6.0/2.0*F_prime*dF_ddelta_prime/(1+F) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  dfde[4] = d2F_ddelta2_prime;
  dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 6.0/2.0*pow(dF_ddelta_prime, 2)/(1+F) + 6.0/2.0*F_prime*d2F_ddelta2_prime/(1+F) - 12.0/2.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 6.0/2.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;

  
  return GSL_SUCCESS;
}

int F_dF_ddF_cylindrical_wrt_delta_jac(double e, const double y[], double *dfdy, double dfde[], void *params){

  integration_parameters *integration_params = (integration_parameters *) params;
  FlatInhomogeneousUniverseLCDM* pointer_to_Universe = integration_params->pointer_to_Universe;

  double F = y[0];
  double F_prime = y[1];
  double dF_ddelta = y[2];
  double dF_ddelta_prime = y[3];
  double d2F_ddelta2 = y[4];
  double d2F_ddelta2_prime = y[5];
  double a = pointer_to_Universe->a_at_eta(e);
  double hubble = pointer_to_Universe->H_at_eta(e);
  double Omega_m = integration_params->Omega_m;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
  gsl_matrix * m = &dfdy_mat.matrix; 
  
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 0, 4, 0.0);
  gsl_matrix_set (m, 0, 5, 0.0);
  
  //dfde[0] = F_prime;
  //dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 3.0/2.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  //dfde[2] = dF_ddelta_prime;
  //dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 6.0/2.0*F_prime*dF_ddelta_prime/(1+F) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  //dfde[4] = d2F_ddelta2_prime;
  //dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 6.0/2.0*pow(dF_ddelta_prime, 2)/(1+F) + 6.0/2.0*F_prime*d2F_ddelta2_prime/(1+F) - 12.0/2.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 6.0/2.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;
  
  gsl_matrix_set (m, 1, 0, 1.5*Omega_m/a*(1.0+2.0*F) - 3.0/2.0*pow(F_prime/(1.0+F), 2.0));
  gsl_matrix_set (m, 1, 1, 6.0/2.0*F_prime/(1+F) - hubble);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 1, 4, 0.0);
  gsl_matrix_set (m, 1, 5, 0.0);
  
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 3, 1.0);
  gsl_matrix_set (m, 2, 4, 0.0);
  gsl_matrix_set (m, 2, 5, 0.0);
  
  //dfde[0] = F_prime;
  //dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 3.0/2.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  //dfde[2] = dF_ddelta_prime;
  //dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 6.0/2.0*F_prime*dF_ddelta_prime/(1+F) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  //dfde[4] = d2F_ddelta2_prime;
  //dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 6.0/2.0*pow(dF_ddelta_prime, 2)/(1+F) + 6.0/2.0*F_prime*d2F_ddelta2_prime/(1+F) - 12.0/2.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 6.0/2.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;
  
  gsl_matrix_set (m, 3, 0, 3.0*Omega_m/a*dF_ddelta - 6.0/2.0*F_prime*dF_ddelta_prime/pow(1.0+F, 2.0) + 6.0/2.0*pow(F_prime, 2)/pow(1.0+F, 3)*dF_ddelta);
  gsl_matrix_set (m, 3, 1, 6.0/2.0*dF_ddelta_prime/(1+F) - 6.0/2.0*F_prime/pow(1+F, 2)*dF_ddelta);
  gsl_matrix_set (m, 3, 2, 1.5*Omega_m/a*(1.0+2.0*F) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2));
  gsl_matrix_set (m, 3, 3, 6.0/2.0*F_prime/(1+F) - hubble);
  gsl_matrix_set (m, 3, 4, 0.0);
  gsl_matrix_set (m, 3, 5, 0.0);
  
  gsl_matrix_set (m, 4, 0, 0.0);
  gsl_matrix_set (m, 4, 1, 0.0);
  gsl_matrix_set (m, 4, 2, 0.0);
  gsl_matrix_set (m, 4, 3, 0.0);
  gsl_matrix_set (m, 4, 4, 0.0);
  gsl_matrix_set (m, 4, 5, 1.0);
  
  //dfde[0] = F_prime;
  //dfde[1] = 1.5*Omega_m/a*F*(1.0+F) + 3.0/2.0*pow(F_prime, 2)/(1+F) - hubble*F_prime;
  //dfde[2] = dF_ddelta_prime;
  //dfde[3] = 1.5*Omega_m/a*(dF_ddelta*(1.0+2.0*F)) + 6.0/2.0*F_prime*dF_ddelta_prime/(1+F) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*dF_ddelta - hubble*dF_ddelta_prime;
  //dfde[4] = d2F_ddelta2_prime;
  //dfde[5] = 1.5*Omega_m/a*(d2F_ddelta2*(1.0+2.0*F) + 2.0*pow(dF_ddelta, 2)) + 6.0/2.0*pow(dF_ddelta_prime, 2)/(1+F) + 6.0/2.0*F_prime*d2F_ddelta2_prime/(1+F) - 12.0/2.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2)*d2F_ddelta2 + 6.0/2.0*pow(F_prime, 2)/pow(1+F, 3)*pow(dF_ddelta, 2.0) - hubble*d2F_ddelta2_prime;
  
  gsl_matrix_set (m, 5, 0, 3.0*Omega_m/a*d2F_ddelta2 - 6.0/2.0*pow(dF_ddelta_prime, 2)/pow(1.0+F,2.0) - 6.0/2.0*F_prime*d2F_ddelta2_prime/pow(1.0+F,2.0) + 24.0/2.0*F_prime*dF_ddelta*dF_ddelta_prime/pow(1+F, 3.0) + 6.0/2.0*pow(F_prime, 2)/pow(1+F, 3)*d2F_ddelta2 - 18.0/2.0*pow(F_prime, 2)/pow(1+F, 4)*pow(dF_ddelta, 2.0));
  gsl_matrix_set (m, 5, 1, 6.0/2.0*d2F_ddelta2_prime/(1+F) - 12.0/2.0*dF_ddelta*dF_ddelta_prime/pow(1+F, 2.0) - 6.0/2.0*F_prime/pow(1+F, 2)*d2F_ddelta2 + 12.0/2.0*F_prime/pow(1+F, 3)*pow(dF_ddelta, 2.0));
  gsl_matrix_set (m, 5, 2, 6.0*Omega_m/a*dF_ddelta - 12.0/2.0*F_prime*dF_ddelta_prime/pow(1+F, 2.0) + 12.0/2.0*pow(F_prime, 2)/pow(1+F, 3)*dF_ddelta);
  gsl_matrix_set (m, 5, 3, 12.0/2.0*dF_ddelta_prime/(1+F) - 12.0/2.0*F_prime*dF_ddelta/pow(1+F, 2.0));
  gsl_matrix_set (m, 5, 4, 1.5*Omega_m/a*(1.0+2.0*F) - 3.0/2.0*pow(F_prime, 2)/pow(1+F, 2));
  gsl_matrix_set (m, 5, 5, 6.0/2.0*F_prime/(1+F) - hubble);
  
  return GSL_SUCCESS;
}
