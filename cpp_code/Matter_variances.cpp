
#include <thread>




void skewness_of_matter_within_R_for_multithread(double R, double n_s, Matter* pointer_to_Matter, double f_NL_rescaling_factor, double alpha_1, double alpha_2, double alpha_3, double* skew){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = n_s;
  params.pointer_to_Matter = pointer_to_Matter;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  (*skew) = f_NL_rescaling_factor*int_gsl_integrate_low_precision(skewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}



void dskewness_of_matter_within_R_dR_for_multithread(double R, double n_s, Matter* pointer_to_Matter, double f_NL_rescaling_factor, double alpha_1, double alpha_2, double alpha_3, double* dskew){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = n_s;
  params.pointer_to_Matter = pointer_to_Matter;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  (*dskew) = f_NL_rescaling_factor*int_gsl_integrate_low_precision(dskewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}

/*******************************************************************************************************************************************************
 * 1.4 set_sphere_variances
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_sphere_variances(){
  
  //int reduction = 8;
  //int n = this->wave_numbers.size()/reduction;
  //int n = 32;
  int n = 512;
  
  // Say we consider the PDF at R=15Mpc/h down to density contrasts of delta_NL = -0.9
  // Then we need R_max = R*(1-0.9)**(1/3) \approx 7Mpc/h --> use 5Mpc/h to have factor of 3.
  double R_min_in_Mpc_over_h = 0.1;//5.0;
  //double R_min_in_Mpc_over_h = 5.0;//5.0;
  // Say we consider the PDF at R=15Mpc/h up to density contrasts of delta_NL = 26
  // Then we need R_max = R*(1+26)**(1/3) = 45Mpc/h.
  double R_max_in_Mpc_over_h = 150.0;//60.0;
  //double R_max_in_Mpc_over_h = 45.0;//60.0;
  vector<double> radii_in_Mpc_over_h(0, 0.0);
  log_binning(R_min_in_Mpc_over_h, R_max_in_Mpc_over_h, n-1, &radii_in_Mpc_over_h);
  
  this->log_top_hat_radii.resize(n, 0.0);
  this->top_hat_sphere_variances.resize(n, 0.0);
  this->dtop_hat_sphere_variances_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  
  for(int i = 0; i < n; i++){
    //this->log_top_hat_radii[i] = log(wave_numbers[reduction*(n-1-i)]);
    this->log_top_hat_radii[i] = log(radii_in_Mpc_over_h[i]/constants::c_over_e5);
    this->top_hat_sphere_variances[i] = this->variance_of_matter_within_R(radii_in_Mpc_over_h[i]/constants::c_over_e5);
    this->dtop_hat_sphere_variances_dR[i] = this->dvariance_of_matter_within_R_dR(radii_in_Mpc_over_h[i]/constants::c_over_e5);
  }
  
}

/*******************************************************************************************************************************************************
 * 1.5 set_sphere_skewnesses
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_sphere_skewnesses(){
  
  int n = this->log_top_hat_radii.size();
  this->log_top_hat_radii_for_skewnesses = this->log_top_hat_radii;
  
  this->top_hat_sphere_local_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_local_skewnesses_dR.resize(n, 0.0);
  this->top_hat_sphere_equilateral_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_equilateral_skewnesses_dR.resize(n, 0.0);
  this->top_hat_sphere_orthogonal_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_orthogonal_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  
  // f_NL is now set in the ipython interface
  double f_NL_local = 1.0;//-0.9-5.1;
  double f_NL_equilateral = 1.0;//-26.0-47.0;
  double f_NL_orthoginal = 1.0;//-38.0-24.0;
  double prefactor = -3.0*this->cosmology.Omega_m/pow(constants::pi2, 4);
  double skew_1, skew_2, skew_3;
  double dskew_1, dskew_2, dskew_3;
  double R;
  
  //FILE *f = fopen("../data_for_PNG_paper/Final_plots/primordial_moments_Nishimishi.dat", "w");
  //fclose(f);
  
  fstream output;
  //output.open("../data_for_PNG_paper/Final_plots/primordial_moments_Nishimishi.dat");
  cout << scientific << setprecision(10);
  for(int i = 0; i < -n; i++){
  //for(int i = n-1; i > -1; i--){
    R = exp(this->log_top_hat_radii[i]);
    cout << i << '\t';
    cout << c_over_e5*R << '\t';
    cout << this->top_hat_sphere_variances[i] << '\t';
    cout << this->dtop_hat_sphere_variances_dR[i] << '\t';
    cout.flush();
    
    std::thread thread_1(skewness_of_matter_within_R_for_multithread, R, this->cosmology.n_s, this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &skew_1);
    //std::thread thread_2(skewness_of_matter_within_R_for_multithread, R, this->cosmology.n_s, this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &skew_2);
    //std::thread thread_3(skewness_of_matter_within_R_for_multithread, R, this->cosmology.n_s, this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &skew_3);
    
    std::thread thread_4(dskewness_of_matter_within_R_dR_for_multithread, R, this->cosmology.n_s, this, this->f_NL_rescaling_factor, 1.0, 1.0, 0.0, &dskew_1);
    //std::thread thread_5(dskewness_of_matter_within_R_dR_for_multithread, R, this->cosmology.n_s, this, this->f_NL_rescaling_factor, 2.0/3.0, 2.0/3.0, 2.0/3.0, &dskew_2);
    //std::thread thread_6(dskewness_of_matter_within_R_dR_for_multithread, R, this->cosmology.n_s, this, this->f_NL_rescaling_factor, 1.0, 1.0/3.0, 2.0/3.0, &dskew_3);
    
    /*
    skew_1 = skewness_of_matter_within_R(R, 1.0, 1.0, 0.0);
    dskew_1 = dskewness_of_matter_within_R_dR(R, 1.0, 1.0, 0.0);
    skew_2 = skewness_of_matter_within_R(R, 2.0/3.0, 2.0/3.0, 2.0/3.0);
    dskew_2 = dskewness_of_matter_within_R_dR(R, 2.0/3.0, 2.0/3.0, 2.0/3.0);
    skew_3 = skewness_of_matter_within_R(R, 1.0, 1.0/3.0, 2.0/3.0);
    dskew_3 = dskewness_of_matter_within_R_dR(R, 1.0, 1.0/3.0, 2.0/3.0);
    */
    
    thread_1.join();
    //thread_2.join();
    //thread_3.join();
    
    thread_4.join();
    //thread_5.join();
    //thread_6.join();
    
    
    this->top_hat_sphere_local_skewnesses[i] = 6.0*f_NL_local*prefactor*skew_1;
    this->dtop_hat_sphere_local_skewnesses_dR[i] = 6.0*f_NL_local*prefactor*dskew_1;
    /*
    this->top_hat_sphere_equilateral_skewnesses[i]  = -18.0*f_NL_equilateral*prefactor*skew_1;
    this->top_hat_sphere_equilateral_skewnesses[i] += -12.0*f_NL_equilateral*prefactor*skew_2;
    this->top_hat_sphere_equilateral_skewnesses[i] +=  36.0*f_NL_equilateral*prefactor*skew_3;
    this->dtop_hat_sphere_equilateral_skewnesses_dR[i]  = -18.0*f_NL_equilateral*prefactor*dskew_1;
    this->dtop_hat_sphere_equilateral_skewnesses_dR[i] += -12.0*f_NL_equilateral*prefactor*dskew_2;
    this->dtop_hat_sphere_equilateral_skewnesses_dR[i] +=  36.0*f_NL_equilateral*prefactor*dskew_3;
    
    this->top_hat_sphere_orthogonal_skewnesses[i]  = -54.0*f_NL_orthoginal*prefactor*skew_1;
    this->top_hat_sphere_orthogonal_skewnesses[i] += -48.0*f_NL_orthoginal*prefactor*skew_2;
    this->top_hat_sphere_orthogonal_skewnesses[i] +=  108.0*f_NL_orthoginal*prefactor*skew_3;
    this->dtop_hat_sphere_orthogonal_skewnesses_dR[i]  = -54.0*f_NL_orthoginal*prefactor*dskew_1;
    this->dtop_hat_sphere_orthogonal_skewnesses_dR[i] += -48.0*f_NL_orthoginal*prefactor*dskew_2;
    this->dtop_hat_sphere_orthogonal_skewnesses_dR[i] +=  108.0*f_NL_orthoginal*prefactor*dskew_3;
    */
    cout << this->top_hat_sphere_local_skewnesses[i] << '\t';
    //cout << this->top_hat_sphere_equilateral_skewnesses[i] << '\t';
    //cout << this->top_hat_sphere_orthogonal_skewnesses[i] << '\t';
    cout << this->dtop_hat_sphere_local_skewnesses_dR[i] << '\t';
    //cout << this->dtop_hat_sphere_equilateral_skewnesses_dR[i] << '\t';
    //cout << this->dtop_hat_sphere_orthogonal_skewnesses_dR[i] << '\t';
    //cout << get_delta0(this->top_hat_sphere_variances[i], this->top_hat_sphere_local_skewnesses[i]) << '\t';
    cout << '\n';
    
    
    output << scientific << setprecision(10);
    
    output << i << setw(20);
    output << c_over_e5*R << setw(20);
    output << this->top_hat_sphere_variances[i] << setw(20);
    output << this->dtop_hat_sphere_variances_dR[i] << setw(20);
    output << this->top_hat_sphere_local_skewnesses[i] << setw(20);
    output << this->top_hat_sphere_equilateral_skewnesses[i] << setw(20);
    output << this->top_hat_sphere_orthogonal_skewnesses[i] << setw(20);
    output << this->dtop_hat_sphere_local_skewnesses_dR[i] << setw(20);
    output << this->dtop_hat_sphere_equilateral_skewnesses_dR[i] << setw(20);
    output << this->dtop_hat_sphere_orthogonal_skewnesses_dR[i] << setw(20);
    output << '\n';
    output.flush();
    
    //std::this_thread::sleep_for(std::chrono::milliseconds(1500));
    
  }
  
  output.close();
  
}


void Matter::set_sphere_skewnesses_from_file(string file_name){
  
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
  
  this->top_hat_sphere_local_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_local_skewnesses_dR.resize(n, 0.0);
  this->top_hat_sphere_equilateral_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_equilateral_skewnesses_dR.resize(n, 0.0);
  this->top_hat_sphere_orthogonal_skewnesses.resize(n, 0.0);
  this->dtop_hat_sphere_orthogonal_skewnesses_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  
  for(int i = 0; i < n; i++){
    
    input_2 >> dummy_double;
    input_2 >> R_in_Mpc_over_h; this->log_top_hat_radii_for_skewnesses[i] = log(R_in_Mpc_over_h/constants::c_over_e5);
    
    input_2 >> dummy_double;
    input_2 >> dummy_double;
    
    input_2 >> top_hat_sphere_local_skewnesses[i];
    input_2 >> top_hat_sphere_equilateral_skewnesses[i];
    input_2 >> top_hat_sphere_orthogonal_skewnesses[i];
    
    input_2 >> dtop_hat_sphere_local_skewnesses_dR[i];
    input_2 >> dtop_hat_sphere_equilateral_skewnesses_dR[i];
    input_2 >> dtop_hat_sphere_orthogonal_skewnesses_dR[i];
    
    cout << dummy_double << '\t';
    cout << R_in_Mpc_over_h << '\t';
    
    cout << top_hat_sphere_variances[i] << '\t';
    cout << dtop_hat_sphere_variances_dR[i] << '\t';
    
    cout << top_hat_sphere_local_skewnesses[i] << '\t';
    cout << top_hat_sphere_equilateral_skewnesses[i] << '\t';
    cout << top_hat_sphere_orthogonal_skewnesses[i] << '\t';
    
    cout << dtop_hat_sphere_local_skewnesses_dR[i] << '\t';
    cout << dtop_hat_sphere_equilateral_skewnesses_dR[i] << '\t';
    cout << dtop_hat_sphere_orthogonal_skewnesses_dR[i] << '\n';
    
  }
  
  input_2.close();
  
}

/*******************************************************************************************************************************************************
 * 1.6 set_cylinder_variances
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_cylinder_variances(){
  
  int n = this->log_top_hat_radii.size();
  double R;
  
  this->log_top_hat_radii.resize(n, 0.0);
  this->top_hat_cylinder_variances.resize(n, 0.0);
  this->dtop_hat_cylinder_variances_dR.resize(n, 0.0);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  
  for(int i = 0; i < n; i++){
    R = exp(this->log_top_hat_radii[i]);
    this->top_hat_cylinder_variances[i] = this->variance_of_matter_within_R_2D(R);
    this->dtop_hat_cylinder_variances_dR[i] = this->dvariance_of_matter_within_R_dR_2D(R);
  }
  
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


double Matter::variance_of_matter_within_R_before_norm_was_determined(double R_in_Mpc_over_h){
 
  double s_sq = 0;
  double prefactors;
  
  integration_parameters params;
  params.top_hat_radius = R_in_Mpc_over_h/c_over_e5;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(norm_derivs_before_norm_was_determined_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double Matter::variance_of_matter_within_R(double R){
 
  double s_sq = 0;
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(norm_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double Matter::variance_of_matter_within_R_NL(double R){
 
  double s_sq = 0;
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(norm_NL_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double Matter::dvariance_of_matter_within_R_dR(double R){
 
  double s_sq = 0;
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  s_sq = int_gsl_integrate_medium_precision(dvar_dR_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(high_k_cutoff_in_H0_units),NULL,1000);
  return s_sq;
  
}

double Matter::skewness_of_matter_within_R(double R, double alpha_1, double alpha_2, double alpha_3){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  
  return this->f_NL_rescaling_factor*int_gsl_integrate_low_precision(skewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}

double Matter::dskewness_of_matter_within_R_dR(double R, double alpha_1, double alpha_2, double alpha_3){
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  params.alpha_1 = alpha_1;
  params.alpha_2 = alpha_2;
  params.alpha_3 = alpha_3;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  return this->f_NL_rescaling_factor*int_gsl_integrate_low_precision(dskewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,1000);
  
}

double Matter::local_skewness_of_matter_within_R(double R){
 
  double skew = 0;
  double f_NL = -0.9;
  double prefactor = -18.0*f_NL*this->cosmology.Omega_m/pow(constants::pi2, 4);
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  skew = prefactor*int_gsl_integrate_low_precision(local_skewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,200);
  return skew;
  
}

double Matter::equilateral_skewness_of_matter_within_R(double R){
 
  double skew = 0;
  double f_NL = -26.0;
  double prefactor = -18.0*f_NL*this->cosmology.Omega_m/pow(constants::pi2, 4);
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  skew = prefactor*int_gsl_integrate_low_precision(equilateral_skewness_integral_1_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,200);
  return skew;
  
}

double Matter::orthogonal_skewness_of_matter_within_R(double R){
 
  double skew = 0;
  double f_NL = -38.0;
  double prefactor = -18.0*f_NL*this->cosmology.Omega_m/pow(constants::pi2, 4);
  
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;
  
  double k_max = min(constants::product_of_kmax_and_R/R, maximal_wave_number_in_H0_units);
  skew = prefactor*int_gsl_integrate_low_precision(orthogonal_skewness_integral_derivs_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(k_max),NULL,200);
  return skew;
  
}


/*******************************************************************************************************************************************************
 * 3.5 variance_of_matter_within_R_2D
 * Description:
 *  - computes variance of density field when averaged over cylinders of radius R and length L = 1. 
 * Arguments:
 *  - R: radius of tophat filter in Mpc/h
 * 
*******************************************************************************************************************************************************/

double Matter::variance_of_matter_within_R_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(norm_derivs_2D_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}

double Matter::dvariance_of_matter_within_R_dR_2D(double R){
 
  double s_sq = 0.0;
  double prefactors;
  
  prefactors = 1.0/(2.0*constants::pi);
  integration_parameters params;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  integration_parameters * pointer_to_params = &params;

  s_sq = int_gsl_integrate_medium_precision(dnorm_dR_derivs_2D_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);

  return prefactors*s_sq;
  
}


double Matter::variance_of_matter_within_R_2D(){
 
  
  double variance = 0.0;

  for(int ell = 1; ell < ell_max; ell++){
    variance += this->current_P_delta_L_in_trough_format[ell]*pow(this->theta_trough_Legendres[ell], 2);
  }

  return variance;
  
}

double Matter::variance_of_matter_within_R_2D_NL(){
 
  
  double variance = 0.0;

  for(int ell = 1; ell < ell_max; ell++){
    variance += this->current_P_delta_NL_in_trough_format[ell]*pow(this->theta_trough_Legendres[ell], 2);
  }

  return variance;
  
}

double Matter::covariance_of_matter_within_R_2D(){

  double covariance = 0.0;

  for(int ell = 1; ell < ell_max; ell++){
    covariance += this->current_P_delta_L_in_trough_format[ell]*this->current_theta1_linear_Legendres[ell]*this->current_theta2_linear_Legendres[ell];
  }

  return covariance;
  
}

double Matter::covariance_of_matter_within_R_2D_NL(){

  double covariance = 0.0;

  for(int ell = 1; ell < ell_max; ell++){
    covariance += this->current_P_delta_NL_in_trough_format[ell]*this->current_theta1_linear_Legendres[ell]*this->current_theta2_linear_Legendres[ell];
  }

  return covariance;
  
}

double Matter::covariance_of_matter_within_R_2D_NL(int bin){

  double covariance = 0.0;

  for(int ell = 1; ell < ell_max; ell++){
    covariance += this->current_P_delta_NL_in_trough_format[ell]*this->bin_Legendres[bin][ell]*this->theta_trough_Legendres[ell];
  }

  return covariance;
  
}


