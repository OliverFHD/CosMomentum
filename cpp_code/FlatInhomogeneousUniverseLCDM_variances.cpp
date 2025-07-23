
#include <thread>




void skewness_of_matter_within_R_for_multithread(double R, double n_s, FlatInhomogeneousUniverseLCDM* pointer_to_Universe, double f_NL_rescaling_factor, double alpha_1, double alpha_2, double alpha_3, double* skew){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = n_s;
  params.pointer_to_Universe = pointer_to_Universe;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  (*skew) = f_NL_rescaling_factor*int_gsl_integrate_low_precision(skewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}



void dskewness_of_matter_within_R_dR_for_multithread(double R, double n_s, FlatInhomogeneousUniverseLCDM* pointer_to_Universe, double f_NL_rescaling_factor, double alpha_1, double alpha_2, double alpha_3, double* dskew){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = n_s;
  params.pointer_to_Universe = pointer_to_Universe;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  (*dskew) = f_NL_rescaling_factor*int_gsl_integrate_low_precision(dskewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}


void skewness_of_matter_within_R_2D_for_multithread(double R, double n_s, FlatInhomogeneousUniverseLCDM* pointer_to_Universe, double f_NL_rescaling_factor, double alpha_1, double alpha_2, double alpha_3, double* skew){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = n_s;
  params.pointer_to_Universe = pointer_to_Universe;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  (*skew) = f_NL_rescaling_factor*int_gsl_integrate_low_precision(skewness_integral_2D_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}



void dskewness_of_matter_within_R_dR_2D_for_multithread(double R, double n_s, FlatInhomogeneousUniverseLCDM* pointer_to_Universe, double f_NL_rescaling_factor, double alpha_1, double alpha_2, double alpha_3, double* dskew){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = n_s;
  params.pointer_to_Universe = pointer_to_Universe;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  (*dskew) = f_NL_rescaling_factor*int_gsl_integrate_low_precision(dskewness_integral_2D_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}

/*******************************************************************************************************************************************************
 * 1.4 set_sphere_variances
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_sphere_variances(){
  
  int n = 512;
  
  double R_min_in_Mpc_over_h = 0.1;
  double R_max_in_Mpc_over_h = 150.0;
  vector<double> radii_in_Mpc_over_h(0, 0.0);
  log_binning(R_min_in_Mpc_over_h, R_max_in_Mpc_over_h, n-1, &radii_in_Mpc_over_h);
  
  this->log_top_hat_radii.resize(n, 0.0);
  this->top_hat_sphere_variances.resize(n, 0.0);
  this->dtop_hat_sphere_variances_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  for(int i = 0; i < n; i++){
    this->log_top_hat_radii[i] = log(radii_in_Mpc_over_h[i]/constants::c_over_e5);
    this->top_hat_sphere_variances[i] = this->variance_of_matter_within_R(radii_in_Mpc_over_h[i]/constants::c_over_e5);
    this->dtop_hat_sphere_variances_dR[i] = this->dvariance_of_matter_within_R_dR(radii_in_Mpc_over_h[i]/constants::c_over_e5);
  }
  
  this->log_top_hat_radii_for_skewnesses = this->log_top_hat_radii;
  
}

/*******************************************************************************************************************************************************
 * 1.5 set_sphere_skewnesses
 * Description:
 *
 * Arguments:
 * - int PNG_modus: 1 == local, 2 == equilateral, 3 == orthogonal
 * - double f_NL: amplitude of bispectrum
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_sphere_skewnesses(int PNG_modus){
  
  int n = this->log_top_hat_radii.size();
  this->log_top_hat_radii_for_skewnesses = this->log_top_hat_radii;
  
  this->top_hat_sphere_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  // f_NL is now set in the ipython interface
  double prefactor = -3.0*this->return_Omega_m()/pow(constants::pi2, 4);//*this->f_NL_rescaling_factor;
  double skew_1, dskew_1;
  double skew_2, dskew_2;
  double skew_3, dskew_3;
  double R;
  
  if(PNG_modus < 1 || PNG_modus > 3){
    cerr << "PNG_modus in set_sphere_skewnesses must take values 1, 2 or 3.";
    exit(1);
  }
    
  
  for(int i = 0; i < n; i++){
    R = exp(this->log_top_hat_radii[i]);
    
    std::thread thread_1(skewness_of_matter_within_R_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &skew_1);
    std::thread thread_4(dskewness_of_matter_within_R_dR_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &dskew_1);
    
    if(PNG_modus != 1){
      std::thread thread_2(skewness_of_matter_within_R_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &skew_2);
      std::thread thread_5(dskewness_of_matter_within_R_dR_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &dskew_2);
      
      std::thread thread_3(skewness_of_matter_within_R_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &skew_3);
      std::thread thread_6(dskewness_of_matter_within_R_dR_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &dskew_3);
      
      thread_2.join();
      thread_5.join();
      
      thread_3.join();
      thread_6.join();
    }
    
    thread_1.join();
    thread_4.join();
    
    if(PNG_modus == 1){
      this->top_hat_sphere_skewnesses[i] = 6.0*prefactor*skew_1;
      this->dtop_hat_sphere_skewnesses_dR[i] = 6.0*prefactor*dskew_1;
    }
    if(PNG_modus == 2){
      this->top_hat_sphere_skewnesses[i]  = -18.0*prefactor*skew_1;
      this->top_hat_sphere_skewnesses[i] += -12.0*prefactor*skew_2;
      this->top_hat_sphere_skewnesses[i] +=  36.0*prefactor*skew_3;
      this->dtop_hat_sphere_skewnesses_dR[i]  = -18.0*prefactor*dskew_1;
      this->dtop_hat_sphere_skewnesses_dR[i] += -12.0*prefactor*dskew_2;
      this->dtop_hat_sphere_skewnesses_dR[i] +=  36.0*prefactor*dskew_3;
    }
    if(PNG_modus == 3){
      this->top_hat_sphere_skewnesses[i]  = -54.0*prefactor*skew_1;
      this->top_hat_sphere_skewnesses[i] += -48.0*prefactor*skew_2;
      this->top_hat_sphere_skewnesses[i] +=  108.0*prefactor*skew_3;
      this->dtop_hat_sphere_skewnesses_dR[i]  = -54.0*prefactor*dskew_1;
      this->dtop_hat_sphere_skewnesses_dR[i] += -48.0*prefactor*dskew_2;
      this->dtop_hat_sphere_skewnesses_dR[i] +=  108.0*prefactor*dskew_3;
    }
    
  }
  
  
  this->skewness_initialisation = INITIALISED;
    
}

/*******************************************************************************************************************************************************
 * 1.6 set_cylinder_skewnesses
 * Description:
 *
 * Arguments:
 * - int PNG_modus: 1 == local, 2 == equilateral, 3 == orthogonal
 * - double f_NL: amplitude of bispectrum
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_cylinder_skewnesses(int PNG_modus){
  
  int n = this->log_top_hat_radii.size();
  this->log_top_hat_cylinder_radii_for_skewnesses = this->log_top_hat_cylinder_radii;
  
  this->top_hat_cylinder_skewnesses.resize(n, 0.0);
  this->dtop_hat_cylinder_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  // f_NL is now set in the ipython interface
  double prefactor = -3.0*this->return_Omega_m()/pow(constants::pi2, 5);//*this->f_NL_rescaling_factor;
  double skew_1, dskew_1;
  double skew_2, dskew_2;
  double skew_3, dskew_3;
  double R;
  
  if(PNG_modus < 1 || PNG_modus > 3){
    cerr << "PNG_modus in set_cylinder_skewnesses must take values 1, 2 or 3.";
    exit(1);
  }
    
  
  for(int i = 0; i < n; i++){
    R = exp(this->log_top_hat_cylinder_radii[i]);
    
    std::thread thread_1(skewness_of_matter_within_R_2D_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &skew_1);
    std::thread thread_4(dskewness_of_matter_within_R_dR_2D_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &dskew_1);
    
    if(PNG_modus != 1){
      std::thread thread_2(skewness_of_matter_within_R_2D_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &skew_2);
      std::thread thread_5(dskewness_of_matter_within_R_dR_2D_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &dskew_2);
      
      std::thread thread_3(skewness_of_matter_within_R_2D_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &skew_3);
      std::thread thread_6(dskewness_of_matter_within_R_dR_2D_for_multithread, R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &dskew_3);
      
      thread_2.join();
      thread_5.join();
      
      thread_3.join();
      thread_6.join();
    }
    
    thread_1.join();
    thread_4.join();
    
    if(PNG_modus == 1){
      this->top_hat_cylinder_skewnesses[i] = 6.0*prefactor*skew_1;
      this->dtop_hat_cylinder_skewnesses_dR[i] = 6.0*prefactor*dskew_1;
    }
    if(PNG_modus == 2){
      this->top_hat_cylinder_skewnesses[i]  = -18.0*prefactor*skew_1;
      this->top_hat_cylinder_skewnesses[i] += -12.0*prefactor*skew_2;
      this->top_hat_cylinder_skewnesses[i] +=  36.0*prefactor*skew_3;
      this->dtop_hat_cylinder_skewnesses_dR[i]  = -18.0*prefactor*dskew_1;
      this->dtop_hat_cylinder_skewnesses_dR[i] += -12.0*prefactor*dskew_2;
      this->dtop_hat_cylinder_skewnesses_dR[i] +=  36.0*prefactor*dskew_3;
    }
    if(PNG_modus == 3){
      this->top_hat_cylinder_skewnesses[i]  = -54.0*prefactor*skew_1;
      this->top_hat_cylinder_skewnesses[i] += -48.0*prefactor*skew_2;
      this->top_hat_cylinder_skewnesses[i] +=  108.0*prefactor*skew_3;
      this->dtop_hat_cylinder_skewnesses_dR[i]  = -54.0*prefactor*dskew_1;
      this->dtop_hat_cylinder_skewnesses_dR[i] += -48.0*prefactor*dskew_2;
      this->dtop_hat_cylinder_skewnesses_dR[i] +=  108.0*prefactor*dskew_3;
    }
    
  }
  
  
  this->skewness_initialisation_2D = INITIALISED;
    
}


void FlatInhomogeneousUniverseLCDM::compute_cylinder_skewnesses_for_unit_L_and_unit_fNL(int PNG_modus, double R_in_Mpc_over_h, double *skew, double *dskew_dR){
  
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  double prefactor = -3.0*this->return_Omega_m()/pow(constants::pi2, 5);
  double skew_1, dskew_1;
  double skew_2, dskew_2;
  double skew_3, dskew_3;
  
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  skewness_of_matter_within_R_2D_for_multithread(R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &skew_1);
  dskewness_of_matter_within_R_dR_2D_for_multithread(R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &dskew_1);
  
  if(PNG_modus != 1){
    skewness_of_matter_within_R_2D_for_multithread(R, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &skew_2);
    dskewness_of_matter_within_R_dR_2D_for_multithread(R, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &dskew_2);
    
    skewness_of_matter_within_R_2D_for_multithread(R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &skew_3);
    dskewness_of_matter_within_R_dR_2D_for_multithread(R, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &dskew_3);
  }
    
  if(PNG_modus == 1){
    (*skew) = 6.0*prefactor*skew_1;
    (*dskew_dR) = 6.0*prefactor*dskew_1;
  }
  if(PNG_modus == 2){
    (*skew)  = -18.0*prefactor*skew_1;
    (*skew) += -12.0*prefactor*skew_2;
    (*skew) +=  36.0*prefactor*skew_3;
    (*dskew_dR)  = -18.0*prefactor*dskew_1;
    (*dskew_dR) += -12.0*prefactor*dskew_2;
    (*dskew_dR) +=  36.0*prefactor*dskew_3;
  }
  if(PNG_modus == 3){
    (*skew)  = -54.0*prefactor*skew_1;
    (*skew) += -48.0*prefactor*skew_2;
    (*skew) +=  108.0*prefactor*skew_3;
    (*dskew_dR)  = -54.0*prefactor*dskew_1;
    (*dskew_dR) += -48.0*prefactor*dskew_2;
    (*dskew_dR) +=  108.0*prefactor*dskew_3;
  }

}


/*******************************************************************************************************************************************************
 * 1.7 set_sphere_skewnesses_from_eps3_powerlaw_approximation
 * Description:
 * - this function sets the skewnesses of the linear density field by approximating
 *   eps3(R) = skew(R)/sigma(R)^3
 *   as a power law in R, i.e.
 *   eps3(R) = A_eps3*(R/R_0)^n_eps3
 *
 * Arguments:
 * - int PNG_modus: 1 == local, 2 == equilateral, 3 == orthogonal
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_sphere_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h){
  
  int n = this->log_top_hat_radii.size();
  this->log_top_hat_radii_for_skewnesses = this->log_top_hat_radii;
  
  this->top_hat_sphere_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  // f_NL is now set in the ipython interface
  double prefactor = -3.0*this->return_Omega_m()/pow(constants::pi2, 4);//*this->f_NL_rescaling_factor;
  double skew_1, dskew_1;
  double skew_2, dskew_2;
  double skew_3, dskew_3;
  double R_0 = R_0_in_Mpc_over_h/constants::c_over_e5;
  double R;
  
  if(PNG_modus < 1 || PNG_modus > 3){
    cerr << "PNG_modus in set_sphere_skewnesses_from_eps3_powerlaw_approximation must take values 1, 2 or 3.";
    exit(1);
  }
  
  std::thread thread_1(skewness_of_matter_within_R_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &skew_1);
  std::thread thread_4(dskewness_of_matter_within_R_dR_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &dskew_1);
  
  if(PNG_modus != 1){
    std::thread thread_2(skewness_of_matter_within_R_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &skew_2);
    std::thread thread_5(dskewness_of_matter_within_R_dR_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &dskew_2);
    
    std::thread thread_3(skewness_of_matter_within_R_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &skew_3);
    std::thread thread_6(dskewness_of_matter_within_R_dR_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &dskew_3);
    
    thread_2.join();
    thread_5.join();
    
    thread_3.join();
    thread_6.join();
  }
  
  thread_1.join();
  thread_4.join();
  
  double skewness, dskewness_dR;
  
  if(PNG_modus == 1){
    skewness = 6.0*prefactor*skew_1;
    dskewness_dR = 6.0*prefactor*dskew_1;
  }
  if(PNG_modus == 2){
    skewness  = -18.0*prefactor*skew_1;
    skewness += -12.0*prefactor*skew_2;
    skewness +=  36.0*prefactor*skew_3;
    dskewness_dR  = -18.0*prefactor*dskew_1;
    dskewness_dR += -12.0*prefactor*dskew_2;
    dskewness_dR +=  36.0*prefactor*dskew_3;
  }
  if(PNG_modus == 3){
    skewness  = -54.0*prefactor*skew_1;
    skewness += -48.0*prefactor*skew_2;
    skewness +=  108.0*prefactor*skew_3;
    dskewness_dR  = -54.0*prefactor*dskew_1;
    dskewness_dR += -48.0*prefactor*dskew_2;
    dskewness_dR +=  108.0*prefactor*dskew_3;
  }
  
  double var, dvar_dR;
  this->return_2nd_moment_and_derivative(R_0, &var, &dvar_dR);
  
  double A_eps3, n_eps3;
  A_eps3 = skewness/pow(var, 1.5);
  n_eps3 = R_0*(dskewness_dR/skewness-1.5*dvar_dR/var);
  
  for(int i = 0; i < n; i++){
    R = exp(this->log_top_hat_radii[i]);
    
    this->top_hat_sphere_skewnesses[i]  = A_eps3*pow(R/R_0, n_eps3);
    this->top_hat_sphere_skewnesses[i] *= pow(this->top_hat_sphere_variances[i],1.5);
    
    this->dtop_hat_sphere_skewnesses_dR[i]  = this->top_hat_sphere_skewnesses[i]*(n_eps3/R + 1.5*this->dtop_hat_sphere_variances_dR[i]/this->top_hat_sphere_variances[i]);
    
  }
  
  this->skewness_initialisation = INITIALISED;
    
}




/*******************************************************************************************************************************************************
 * 1.7 set_cylinder_skewnesses_from_eps3_powerlaw_approximation
 * Description:
 * - this function sets the skewnesses of the linear density field by approximating
 *   eps3(R) = skew(R)/sigma(R)^3
 *   as a power law in R, i.e.
 *   eps3(R) = A_eps3*(R/R_0)^n_eps3
 *
 * Arguments:
 * - int PNG_modus: 1 == local, 2 == equilateral, 3 == orthogonal
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_cylinder_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h){
  
  int n = this->log_top_hat_cylinder_radii.size();
  this->log_top_hat_cylinder_radii_for_skewnesses = this->log_top_hat_cylinder_radii;
  
  this->top_hat_cylinder_skewnesses.resize(n, 0.0);
  this->dtop_hat_cylinder_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  // f_NL is now set in the ipython interface
  double prefactor = -3.0*this->return_Omega_m()/pow(constants::pi2, 3);//*this->f_NL_rescaling_factor;
  double skew_1, dskew_1;
  double skew_2, dskew_2;
  double skew_3, dskew_3;
  double R_0 = R_0_in_Mpc_over_h/constants::c_over_e5;
  double R;
  
  if(PNG_modus < 1 || PNG_modus > 3){
    cerr << "PNG_modus in set_cylinder_skewnesses_from_eps3_powerlaw_approximation must take values 1, 2 or 3.";
    exit(1);
  }
  
  cout << "Starting jobs.\n";
  std::thread thread_1(skewness_of_matter_within_R_2D_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &skew_1);
  std::thread thread_4(dskewness_of_matter_within_R_dR_2D_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &dskew_1);
  
  if(PNG_modus != 1){
    std::thread thread_2(skewness_of_matter_within_R_2D_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &skew_2);
    std::thread thread_5(dskewness_of_matter_within_R_dR_2D_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &dskew_2);
    
    std::thread thread_3(skewness_of_matter_within_R_2D_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &skew_3);
    std::thread thread_6(dskewness_of_matter_within_R_dR_2D_for_multithread, R_0, this->return_n_s(), this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &dskew_3);
    
    thread_2.join();
    thread_5.join();
    
    thread_3.join();
    thread_6.join();
  }
  cout << "Jobs started.\n";
  
  thread_1.join();
  thread_4.join();
  
  cout << "Jobs re-joined.\n";
  
  double skewness, dskewness_dR;
  
  if(PNG_modus == 1){
    skewness = 6.0*prefactor*skew_1;
    dskewness_dR = 6.0*prefactor*dskew_1;
  }
  if(PNG_modus == 2){
    skewness  = -18.0*prefactor*skew_1;
    skewness += -12.0*prefactor*skew_2;
    skewness +=  36.0*prefactor*skew_3;
    dskewness_dR  = -18.0*prefactor*dskew_1;
    dskewness_dR += -12.0*prefactor*dskew_2;
    dskewness_dR +=  36.0*prefactor*dskew_3;
  }
  if(PNG_modus == 3){
    skewness  = -54.0*prefactor*skew_1;
    skewness += -48.0*prefactor*skew_2;
    skewness +=  108.0*prefactor*skew_3;
    dskewness_dR  = -54.0*prefactor*dskew_1;
    dskewness_dR += -48.0*prefactor*dskew_2;
    dskewness_dR +=  108.0*prefactor*dskew_3;
  }
  cout << "Computing variance.\n";
  
  double var, dvar_dR;
  this->return_2nd_moment_and_derivative_2D(R_0, &var, &dvar_dR);
    
  cout << "Variance computed.\n";
  
  double A_eps3, n_eps3;
  double alpha = 2.0;
  A_eps3 = skewness/pow(var, alpha);
  n_eps3 = R_0*(dskewness_dR/skewness-alpha*dvar_dR/var);
  
  cout << "Setting arrays.\n";
  for(int i = 0; i < n; i++){
    R = exp(this->log_top_hat_cylinder_radii[i]);
    
    this->top_hat_cylinder_skewnesses[i]  = A_eps3*pow(R/R_0, n_eps3);
    this->top_hat_cylinder_skewnesses[i] *= pow(this->top_hat_cylinder_variances[i],alpha);
    this->dtop_hat_cylinder_skewnesses_dR[i]  = this->top_hat_cylinder_skewnesses[i]*(n_eps3/R + alpha*this->dtop_hat_cylinder_variances_dR[i]/this->top_hat_cylinder_variances[i]);
  }
  cout << "Arrays set.\n";
  
  this->skewness_initialisation_2D = INITIALISED;
    
}


void FlatInhomogeneousUniverseLCDM::set_sphere_skewnesses_from_file(string file_name){
  
  fstream input_1, input_2;
  input_1.open(file_name);
  
  int n = 0;
  double dummy_double;
  double R_in_Mpc_over_h;
  string dummy_string;
  while ( getline(input_1, dummy_string) ) ++n;
  n -= 1;
  
  input_1.close();
  input_2.open(file_name);
  getline(input_2, dummy_string);
  
  this->log_top_hat_radii_for_skewnesses.resize(n, 0.0);
  
  this->top_hat_sphere_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  for(int i = 0; i < n; i++){
    input_2 >> dummy_double;
    input_2 >> R_in_Mpc_over_h; this->log_top_hat_radii_for_skewnesses[i] = log(R_in_Mpc_over_h/constants::c_over_e5);
    
    input_2 >> top_hat_sphere_skewnesses[i];
    input_2 >> dtop_hat_sphere_skewnesses_dR[i];
  }
  
  input_2.close();
  
  this->skewness_initialisation = INITIALISED;
  
}



/*******************************************************************************************************************************************************
 * 1.7 set_cylinder_variances
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/
/*
void FlatInhomogeneousUniverseLCDM::set_cylinder_variances(){
  
  //int n = 128;
  int n = 512;
  // ISSUE: this should not be hardcoded locally
  
  //double R_min_in_Mpc_over_h = 0.1;
  //double R_max_in_Mpc_over_h = 150.0;
  double R_min_in_Mpc_over_h = 0.01;
  double R_max_in_Mpc_over_h = 500.0;
  //double R_min_in_Mpc_over_h = 1.0/0.518697E+02;
  //double R_max_in_Mpc_over_h = 1.0/0.158671E-04;
  vector<double> radii_in_Mpc_over_h(0, 0.0);
  log_binning(R_min_in_Mpc_over_h, R_max_in_Mpc_over_h, n-1, &radii_in_Mpc_over_h);
  
  this->log_top_hat_cylinder_radii.resize(n, 0.0);
  this->top_hat_cylinder_variances.resize(n, 0.0);
  this->dtop_hat_cylinder_variances_dR.resize(n, 0.0);
  this->d2top_hat_cylinder_variances_dR2.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  for(int i = 0; i < n; i++){
    this->log_top_hat_cylinder_radii[i] = log(radii_in_Mpc_over_h[i]/constants::c_over_e5);
    this->top_hat_cylinder_variances[i] = this->variance_of_matter_within_R_2D(radii_in_Mpc_over_h[i]/constants::c_over_e5);
    this->dtop_hat_cylinder_variances_dR[i] = this->dvariance_of_matter_within_R_dR_2D(radii_in_Mpc_over_h[i]/constants::c_over_e5);
  }
  
}
*/


void FlatInhomogeneousUniverseLCDM::set_cylinder_variances(){
  
  int n = this->wave_numbers.size()/8;
  
  this->log_top_hat_cylinder_radii.resize(n, 0.0);
  this->top_hat_cylinder_variances.resize(n, 0.0);
  this->average_of_squared_cylinder_saddle_point.resize(n, 0.0);
  this->dtop_hat_cylinder_variances_dR.resize(n, 0.0);
  this->d2top_hat_cylinder_variances_dR2.resize(n, 0.0);
  this->current_P_L = this->P_L(this->eta_at_a(1.0));
  
  vector<double> R_values(n, 0.0);
  
  
  for(int i = 0; i < n; i++){
    R_values[i] = 1.0/wave_numbers[8*(n-1-i)];
    this->log_top_hat_cylinder_radii[i] = -log(wave_numbers[8*(n-1-i)]);
    this->top_hat_cylinder_variances[i] = this->variance_of_matter_within_R_2D(1.0/wave_numbers[8*(n-1-i)]);
    //this->average_of_squared_cylinder_saddle_point[i] = this->average_of_squared_saddle_point_within_R_2D(1.0/wave_numbers[8*(n-1-i)]);
    //this->average_of_squared_cylinder_saddle_point[i] /= pow(this->top_hat_cylinder_variances[i],2);
    this->dtop_hat_cylinder_variances_dR[i] = this->dvariance_of_matter_within_R_dR_2D(1.0/wave_numbers[8*(n-1-i)]);
  }
  /*
  
  vector<double> P_L_2D = this->current_P_L;
  int n_k = P_L_2D.size();
  double L = 150.0/constants::c_over_e5;
  
  for(int i = 0; i < n_k; i++){
    double prefactors = L/constants::pi;
    integration_parameters params;
    params.top_hat_length = L;
    params.k = wave_numbers[i];
    params.n_s = this->return_n_s();
    params.pointer_to_Universe = this;
    integration_parameters * pointer_to_params = &params;
    P_L_2D[i] = prefactors*int_gsl_integrate_medium_precision(P_to_P2D_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
    cout << i << "  " << wave_numbers[i]/constants::c_over_e5 << "  " <<  P_L_2D[i] << "  " <<  this->current_P_L[i] << '\n';
  }
  
  this->current_P_L = P_L_2D;
  
  for(int i = 0; i < n; i++){
    R_values[i] = 1.0/wave_numbers[8*(n-1-i)];
    this->log_top_hat_cylinder_radii[i] = -log(wave_numbers[8*(n-1-i)]);
    this->top_hat_cylinder_variances[i] = this->variance_of_matter_within_R_2D(1.0/wave_numbers[8*(n-1-i)]);
    this->dtop_hat_cylinder_variances_dR[i] = this->dvariance_of_matter_within_R_dR_2D(1.0/wave_numbers[8*(n-1-i)]);
    cout << i << "  " << R_values[i]*constants::c_over_e5 << "  " <<  this->top_hat_cylinder_variances[i] << "  " <<  this->dtop_hat_cylinder_variances_dR[i] << '\n';
  }
  */
  
  for(int i = 0; i < n; i++){
    this->d2top_hat_cylinder_variances_dR2[i] = interpolate_neville_aitken_derivative(R_values[i], &R_values, &this->dtop_hat_cylinder_variances_dR, constants::order_of_interpolation);
  }
  
  
  this->log_top_hat_cylinder_radii_for_skewnesses = this->log_top_hat_cylinder_radii;
  
}


/*******************************************************************************************************************************************************
 * 3.3 variance_of_matter_within_R
 * Description:
 *  - computes variance of density field when smoothed with tophat filter. 
 * Arguments:
 *  - R_in_Mpc_over_h: radius of tophat filter in Mpc/h
 * 
 * 
 * NOTE: x * Mpc/h = x * Mpc/H_0/Mpc*100*km/s
 *                 = x * 10^5 m/s/H_0
 *                 = x * 10^5 m/s/c * c/H_0
 *                 = x / c_over_e5 * c/H_0
 * 
*******************************************************************************************************************************************************/


double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R_before_norm_was_determined(double R_in_Mpc_over_h){
 
  double s_sq = 0;
  double prefactors;
  
  integration_parameters params;
  params.top_hat_radius = R_in_Mpc_over_h/c_over_e5;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(norm_derivs_before_norm_was_determined_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R(double R){
 
  double s_sq = 0;
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(norm_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R_NL(double R){
 
  double s_sq = 0;
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(norm_NL_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::dvariance_of_matter_within_R_dR(double R){
 
  double s_sq = 0;
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(dvar_dR_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}


/*******************************************************************************************************************************************************
 * 3.5 variance_of_matter_within_R_2D
 * Description:
 *  - computes variance of density field when averaged over cylinders of radius R and length L = 1. 
 * Arguments:
 *  - R: radius of tophat filter in Mpc/h
 * 
*******************************************************************************************************************************************************/

double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(var_derivs_2D_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::average_of_squared_saddle_point_within_R_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 2.0/pow(R,2)/pow(2.0*constants::pi, 2);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(squared_saddle_point_derivs_2D_gsl,(void*)pointer_to_params,1.0/maximal_wave_number_in_H0_units,R,NULL,1000);

  return prefactors*s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R_NL_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(var_NL_derivs_2D_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R_NL_2D(double R, double L){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi_sq); // now pi_sq instead of pi, because of how sinc^2 approximates the delta function.
  integration_parameters params;
  params.top_hat_radius = R;
  params.top_hat_length = L;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(var_NL_derivs_2D_1_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::variance_of_matter_within_R_2D(double R, double L){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi_sq); // now pi_sq instead of pi, because of how sinc^2 approximates the delta function.
  integration_parameters params;
  params.top_hat_radius = R;
  params.top_hat_length = L;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(var_derivs_2D_1_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}

double FlatInhomogeneousUniverseLCDM::dvariance_of_matter_within_R_dR_2D(double R, double L){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi_sq); // now pi_sq instead of pi, because of how sinc^2 approximates the delta function.
  integration_parameters params;
  params.top_hat_radius = R;
  params.top_hat_length = L;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(dvar_dR_derivs_2D_1_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}


double FlatInhomogeneousUniverseLCDM::dvariance_of_matter_within_R_dR_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(dvar_dR_derivs_2D_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}


double FlatInhomogeneousUniverseLCDM::d2variance_of_matter_within_R_dR2_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->return_n_s();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(d2var_dR2_derivs_2D_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}
