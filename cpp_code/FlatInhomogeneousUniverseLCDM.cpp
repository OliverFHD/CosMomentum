
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_legendre.h>

#include "Eisenstein_and_Hu.h"
#include "generating_function_utils.h"
#include "FlatInhomogeneousUniverseLCDM_integrands.cpp"
#include "FlatInhomogeneousUniverseLCDM_initial_conditions.cpp"
#include "FlatInhomogeneousUniverseLCDM_analytic_special_cases.cpp"
#include "FlatInhomogeneousUniverseLCDM_halofit.cpp"
#include "FlatInhomogeneousUniverseLCDM_variances.cpp"
#include "FlatInhomogeneousUniverseLCDM_output.cpp"


/*******************************************************************************************************************************************************
 * 1.1
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

FlatInhomogeneousUniverseLCDM::FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max) : FlatHomogeneousUniverseLCDM(cosmo, a_min, a_max){
  
  this->set_initial_conditions_for_growth();
  this->initialize_linear_growth_factor_of_delta();
  //this->print_background_cosmology("expansion_history.dat");
  this->initialize_up_to_second_order_growth_factor_of_delta(2.0/7.0*this->D_initial*this->D_initial, 4.0/7.0*this->D_initial*this->D_prime_initial);
  
  // ISSUE: terrible that this is hard coded!!! At least put it into constants.h
  double z_of_f_NL_rescale = 500.0;
  double a_of_f_NL_rescale = 1.0/(1.0+z_of_f_NL_rescale);
  double eta_of_f_NL_rescale = this->eta_at_a(a_of_f_NL_rescale);
  this->f_NL_rescaling_factor = this->return_D_of_eta(eta_of_f_NL_rescale)/a_of_f_NL_rescale;
  
  this->set_wave_numbers();
  if(this->return_Omega_b() != 0.0){
    this->set_transfer_function_Eisenstein_and_Hu();
  }
  else{
    this->set_transfer_function_Bond_and_Efstathiou();
  }
  //this->set_transfer_function_from_file("./Data/Aniks_new_Nbody_sims/transfer_function_z0.dat");
  //cout << "Set Anik's tranfer function.\n";
  this->norm = this->variance_of_matter_within_R_before_norm_was_determined(8.0);
  cout << "Norm = " << this->norm << '\n';
  cout << "Norm/c_over_e5 = " << this->norm/constants::c_over_e5 << '\n';
  
  cout << "Setting sphere variances.\n";
  this->set_sphere_variances();
  cout << "Setting cylinder variances.\n";
  this->set_cylinder_variances();
  this->skewness_initialisation = UNINITIALISED;
  this->skewness_initialisation_2D = UNINITIALISED;
  
  cout << "Setting spherical collapse evolution.\n";
  this->set_spherical_collapse_evolution_of_delta();
  cout << "Setting cylindrical collapse evolution.\n";
  this->set_cylindrical_collapse_evolution_of_delta();
  cout << "Done.\n";
  
  
}


FlatInhomogeneousUniverseLCDM::FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max, int reduced_computation) : FlatHomogeneousUniverseLCDM(cosmo, a_min, a_max){
  
  this->set_initial_conditions_for_growth();
  this->initialize_linear_growth_factor_of_delta();
  this->initialize_up_to_second_order_growth_factor_of_delta(2.0/7.0*this->D_initial*this->D_initial, 4.0/7.0*this->D_initial*this->D_prime_initial);
  
  // ISSUE: terrible that this is hard coded!!! At least put it into constants.h
  double z_of_f_NL_rescale = 500.0;
  double a_of_f_NL_rescale = 1.0/(1.0+z_of_f_NL_rescale);
  double eta_of_f_NL_rescale = this->eta_at_a(a_of_f_NL_rescale);
  this->f_NL_rescaling_factor = this->return_D_of_eta(eta_of_f_NL_rescale)/a_of_f_NL_rescale;
  
  this->set_wave_numbers();
  if(this->return_Omega_b() != 0.0){
    this->set_transfer_function_Eisenstein_and_Hu();
  }
  else{
    this->set_transfer_function_Bond_and_Efstathiou();
  }
  this->norm = this->variance_of_matter_within_R_before_norm_was_determined(8.0);
  
  if(reduced_computation == 0){
    cout << "Setting sphere variances.\n";
    this->set_sphere_variances();
    cout << "Setting cylinder variances.\n";
    this->set_cylinder_variances();
    this->skewness_initialisation = UNINITIALISED;
    this->skewness_initialisation_2D = UNINITIALISED;
    
    cout << "Setting spherical collapse evolution.\n";
    this->set_spherical_collapse_evolution_of_delta();
    cout << "Setting cylindrical collapse evolution.\n";
    this->set_cylindrical_collapse_evolution_of_delta();
    cout << "Done.\n";
  }
  
  
}

FlatInhomogeneousUniverseLCDM::FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max, string file_for_transfer_function) : FlatHomogeneousUniverseLCDM(cosmo, a_min, a_max){
  
  this->set_initial_conditions_for_growth();
  this->initialize_linear_growth_factor_of_delta();
  this->print_growth_history("linear_growth_history.dat");
  this->initialize_up_to_second_order_growth_factor_of_delta(2.0/7.0*this->D_initial*this->D_initial, 4.0/7.0*this->D_initial*this->D_prime_initial);
  
  // ISSUE: terrible that this is hard coded!!! At least put it into constants.h
  double z_of_f_NL_rescale = 500.0;
  double a_of_f_NL_rescale = 1.0/(1.0+z_of_f_NL_rescale);
  double eta_of_f_NL_rescale = this->eta_at_a(a_of_f_NL_rescale);
  this->f_NL_rescaling_factor = this->return_D_of_eta(eta_of_f_NL_rescale)/a_of_f_NL_rescale;
  
  this->set_wave_numbers();
  this->set_transfer_function_from_file(file_for_transfer_function);
  this->norm = this->variance_of_matter_within_R_before_norm_was_determined(8.0);
  
  cout << "Setting sphere variances.\n";
  this->set_sphere_variances();
  cout << "Setting cylinder variances.\n";
  this->set_cylinder_variances();
  this->skewness_initialisation = UNINITIALISED;
  this->skewness_initialisation_2D = UNINITIALISED;
  
  cout << "Setting spherical collapse evolution.\n";
  this->set_spherical_collapse_evolution_of_delta();
  cout << "Setting cylindrical collapse evolution.\n";
  this->set_cylindrical_collapse_evolution_of_delta();
  cout << "Done.\n";
  
}

/*******************************************************************************************************************************************************
 * 1.2 Destructor
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

FlatInhomogeneousUniverseLCDM::~FlatInhomogeneousUniverseLCDM(){
}


/*******************************************************************************************************************************************************
 * 1.3 set_wave_numbers
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_wave_numbers(){
  
  log_binning(minimal_wave_number_in_H0_units, maximal_wave_number_in_H0_units, constants::number_of_k-1, &this->wave_numbers);

  int n = this->wave_numbers.size();
  this->log_wave_numbers.resize(n);
  for(int i = 0; i<n; i++){
    this->log_wave_numbers[i] = log(this->wave_numbers[i]);
  }

}

/*******************************************************************************************************************************************************
 * 1.4 set_transfer_function_Eisenstein_and_Hu
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_transfer_function_Eisenstein_and_Hu(){
  
  tranfer_function_Eisenstein_and_Hu(&this->wave_numbers, &this->transfer_function, this->return_Omega_m(), this->return_Omega_b(), this->return_h_100(), this->return_theta_27());
  
}

/*******************************************************************************************************************************************************
 * 1.5 set_transfer_function_Eisenstein_and_Hu
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_transfer_function_Bond_and_Efstathiou(){
  
  tranfer_function_Bond_and_Efstathiou(&this->wave_numbers, &this->transfer_function, this->return_Omega_m(), this->return_h_100());
  
}

/*******************************************************************************************************************************************************
 * 1.6 set_transfer_function_from_file
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_transfer_function_from_file(string file){
  
  fstream input;
  input.open(file);
  double k, Tk, Tk_0;
  vector<double> logk_values(0,0.0);
  vector<double> Tk_values(0,0.0);
  
  this->transfer_function = this->wave_numbers;
  
  
  while(input.good()){
    input >> k;
    k *= constants::c_over_e5;
    if(input.good()){
      input >> Tk;
      if(logk_values.size() == 0) Tk_0 = Tk;
      logk_values.push_back(log(k));
      Tk_values.push_back(Tk/Tk_0);
    }
  }
  
  for(int i = 0; i < this->wave_numbers.size(); i++){
    this->transfer_function[i] = interpolate_neville_aitken(this->log_wave_numbers[i], &logk_values, &Tk_values, constants::order_of_interpolation);
  }
  
  input.close();
}


/*************************************************************
 *************************************************************
 **__________ 3. GROWTH FACTORS AND POWER SPECTRA __________**
 *************************************************************
 ************************************************************* 
 *                                                           *
 * ..... 3.1 initialize_linear_growth_factor_of_delta        *
 * ..... 3.2 Newtonian_linear_power_spectrum                 *
 * ..... 3.3 variance_of_matter_within_R                     *
 * ..... 3.4 derivative_variance_of_matter_within_R          *
 * ..... more ......                                         *
 *                                                           *
 *************************************************************
 *************************************************************/


/*******************************************************************************************************************************************************
 * 3.1 initialize_linear_growth_factor_of_delta
 * Description:
 *
 * Arguments:
 * 
 * Comments:
 * - In this method there are still NO radiation inhomogeneities includes.
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::initialize_linear_growth_factor_of_delta(){

  integration_parameters params;
  params.pointer_to_Universe = this;
  params.Omega_m = this->return_Omega_m();
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {growth_factor_gsl, growth_factor_gsl_jac, 2, (void *) pointer_to_params};
  
  int number_of_time_steps = this->return_number_of_time_steps();
  this->Newtonian_growth_factor_of_delta.resize(number_of_time_steps);
  this->Newtonian_growth_factor_of_delta_prime.resize(number_of_time_steps);
  
  this->Newtonian_growth_factor_of_delta[0] = this->D_initial;
  this->Newtonian_growth_factor_of_delta_prime[0] = this->D_prime_initial;
  
  double e_i = this->eta_initial;
  double e_f;
  double y[2] = { this->Newtonian_growth_factor_of_delta[0], this->Newtonian_growth_factor_of_delta_prime[0]};
  
  // All quantities are of order 1.0 at their maximum.
  double hstart = 1.0*constants::gsl_hstart_relative;
  double eps_absolute = 1.0*constants::gsl_eps_relative;
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, eps_absolute, constants::gsl_eps_relative);

  for (int i = 1; i < number_of_time_steps; i++){
    
    e_f = this->eta[i];
    int status = gsl_odeiv2_driver_apply(d, &e_i, e_f, y);
    
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    
    this->Newtonian_growth_factor_of_delta[i] = y[0];
    this->Newtonian_growth_factor_of_delta_prime[i] = y[1];
  }

  gsl_odeiv2_driver_free(d);
  
  double e0 = this->eta_at_a(1.0);
  double D0 = interpolate_neville_aitken(e0, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
  this->D_initial /= D0;
  this->D_prime_initial /= D0;
    
  for (int i = 0; i < number_of_time_steps; i++){
    this->Newtonian_growth_factor_of_delta[i] /= D0;
    this->Newtonian_growth_factor_of_delta_prime[i] /= D0;
  }

}


/*******************************************************************************************************************************************************
 * 3.1 initialize_up_to_second_order_growth_factor_of_delta
 * Description:
 *
 * Arguments:
 * 
 * Comments:
 * - In this method there are still NO radiation inhomogeneities includes. This would also require to take into account
 *   entropy perturbation!
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::initialize_up_to_second_order_growth_factor_of_delta(double D, double D_prime){

  integration_parameters params;
  params.pointer_to_Universe = this;
  params.Omega_m = this->return_Omega_m();
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {growth_factor_to_second_order_gsl, growth_factor_to_second_order_gsl_jac, 2, (void *) pointer_to_params};

  int number_of_time_steps = this->return_number_of_time_steps();
  this->Newtonian_growth_factor_second_order.resize(number_of_time_steps);
  this->Newtonian_growth_factor_second_order_prime.resize(number_of_time_steps);
  
  this->Newtonian_growth_factor_second_order[0] = D;
  this->Newtonian_growth_factor_second_order_prime[0] = D_prime;
  
  double e;
  double e_i = this->eta[0];
  double e_f;
  double y[2] = { this->Newtonian_growth_factor_second_order[0], this->Newtonian_growth_factor_second_order_prime[0]};
    
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

  for (int i = 1; i < number_of_time_steps; i++){

    e_f = this->eta[i];
    int status = gsl_odeiv2_driver_apply(d, &e_i, e_f, y);
    
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    
    this->Newtonian_growth_factor_second_order[i] = y[0];
    this->Newtonian_growth_factor_second_order_prime[i] = y[1];
    
  }

  gsl_odeiv2_driver_free(d);
}

/*******************************************************************************************************************************************************
 * 3.2 Newtonian_linear_power_spectrum
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double FlatInhomogeneousUniverseLCDM::Newtonian_linear_power_spectrum(double k, double e){

  
  double D = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double transfer_sq = this->transfer_function_at(k);
  transfer_sq *= transfer_sq;
  return D*D*pow(k, this->return_n_s())*transfer_sq*pow(this->return_sigma_8(), 2)/this->norm;
    
}

/*******************************************************************************************************************************************************
 * 5.13 set_spherical_collapse_evolution_of_delta
 * Description:
 * - for a number of initial values for delta, this function computes their evolution in the
 *   spherical collapse model
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_spherical_collapse_evolution_of_delta(){

  double delta_min = -10.0;
  //double delta_max = 1.67;
  //double delta_max = 5.0;
  double delta_max = 15.0;
  double ddelta = 0.01;
  double delta = delta_min-ddelta;
  
  /* ISSUE: THIS HARDCODING IS REALLY BAD!! */
  /* CHANGE REQUIRED                        */
  double a_i = this->a[0];
  double D_i = this->Newtonian_growth_factor_of_delta[0];
  double e_i = this->eta[0];
  double D_prime_i = this->Newtonian_growth_factor_of_delta_prime[0];
  if(this->Newtonian_growth_factor_of_delta[0] > 0.001){
    a_i = a_i*0.001/D_i;
    D_i = 0.001;
    double t_i, H_i, H_prime_i;
    this->expansion_in_flat_matter_dominated_universe(a_i, this->return_Omega_m(), &t_i, &e_i, &H_i, &H_prime_i);
    D_prime_i = D_i*H_i;
  }
  
  int n_delta = 1+int((delta_max - delta_min)/ddelta);
  int number_of_time_steps = this->return_number_of_time_steps();
  this->delta_values_for_spherical_collapse.resize(n_delta, 0.0);
  this->spherical_collapse_evolution_of_delta.resize(n_delta, vector<double>(number_of_time_steps, 0.0));
  this->spherical_collapse_evolution_of_delta_ddelta.resize(n_delta, vector<double>(number_of_time_steps, 0.0));
  this->spherical_collapse_evolution_of_delta_ddelta2.resize(n_delta, vector<double>(number_of_time_steps, 0.0));
  
  integration_parameters params;
  params.Omega_m = this->return_Omega_m();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;
  
  
  gsl_odeiv2_system sys = {F_dF_ddF_spherical_wrt_delta_gsl, F_dF_ddF_spherical_wrt_delta_jac, 6, (void *) pointer_to_params};
  for(int i = 0; i < n_delta; i++){
    
    delta += ddelta;
    delta_values_for_spherical_collapse[i] = delta;
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    double e_start = e_i;
    double y[6] = { delta*D_i, delta*D_prime_i, D_i, D_prime_i, 0.0, 0.0};
    for (int t = 0; t < number_of_time_steps; t++){
      double e_f = this->eta[t];
      if(t > 0 && this->spherical_collapse_evolution_of_delta[i][t-1] < constants::max_contrast){
        int status = gsl_odeiv2_driver_apply(d, &e_start, e_f, y);
        if (status != GSL_SUCCESS){
          printf ("error, return value=%d\n", status);
          exit(1);
        }
        this->spherical_collapse_evolution_of_delta[i][t] = y[0];
        this->spherical_collapse_evolution_of_delta_ddelta[i][t] = y[2];
        this->spherical_collapse_evolution_of_delta_ddelta2[i][t] = y[4];
      }
      else if(t>0){
        this->spherical_collapse_evolution_of_delta[i][t] = constants::max_contrast;
        this->spherical_collapse_evolution_of_delta_ddelta[i][t] = 0.0;
        this->spherical_collapse_evolution_of_delta_ddelta2[i][t] = 0.0;
        
      }
    }

    gsl_odeiv2_driver_free(d);
    
  }
  

  
}

/*******************************************************************************************************************************************************
 * 5.14 set_cylindrical_collapse_evolution_of_delta
 * Description:
 * - for a number of initial values for delta, this function computes their evolution in the
 *   cylindrical collapse model
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::set_cylindrical_collapse_evolution_of_delta(){

  //double delta_min = -10.0;
  //double delta_min = -3.0;
  double delta_min = -1.5;
  double delta_max = 1.455;
  delta_max = 2.0;
  //double delta_max = 1.4;
  //double ddelta = 0.001;
  double ddelta = 0.002;
  //double ddelta = 0.005;
  //double ddelta = 0.02;
  double delta = delta_min-ddelta;
  
  /* ISSUE: THIS HARDCODING IS REALLY BAD!! */
  /* CHANGE REQUIRED                        */
  double a_i = this->a[0];
  double D_i = this->Newtonian_growth_factor_of_delta[0];
  double e_i = this->eta[0];
  double D_prime_i = this->Newtonian_growth_factor_of_delta_prime[0];
  if(this->Newtonian_growth_factor_of_delta[0] > 0.001){
    a_i = a_i*0.001/D_i;
    D_i = 0.001;
    double t_i, H_i, H_prime_i;
    this->expansion_in_flat_matter_dominated_universe(a_i, this->return_Omega_m(), &t_i, &e_i, &H_i, &H_prime_i);
    D_prime_i = D_i*H_i;
  }
  
  int n_delta = 1+int((delta_max - delta_min)/ddelta);
  int number_of_time_steps = this->return_number_of_time_steps();
  this->delta_values_for_cylindrical_collapse.resize(n_delta, 0.0);
  this->cylindrical_collapse_evolution_of_delta.resize(n_delta, vector<double>(number_of_time_steps, 0.0));
  this->cylindrical_collapse_evolution_of_delta_ddelta.resize(n_delta, vector<double>(number_of_time_steps, 0.0));
  this->cylindrical_collapse_evolution_of_delta_ddelta2.resize(n_delta, vector<double>(number_of_time_steps, 0.0));
  
  integration_parameters params;
  params.Omega_m = this->return_Omega_m();
  params.pointer_to_Universe = this;
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {F_dF_ddF_cylindrical_wrt_delta_gsl, F_dF_ddF_cylindrical_wrt_delta_jac, 6, (void *) pointer_to_params};
  for(int d = 0; d < n_delta; d++){
    
    delta += ddelta;
    delta_values_for_cylindrical_collapse[d] = delta;
    
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    double e_start = e_i;
    double y[6] = { delta*D_i, delta*D_prime_i, D_i, D_prime_i, 0.0, 0.0};
    int status;
    for (int t = 0; t < number_of_time_steps; t++){
      double e_f = this->eta[t];
      if(t > 0 && this->cylindrical_collapse_evolution_of_delta[d][t-1] < constants::max_contrast){
        status = gsl_odeiv2_driver_apply(driver, &e_start, e_f, y);
        if (status != GSL_SUCCESS){
          printf ("error, return value=%d\n", status);
          exit(1);
        }
        this->cylindrical_collapse_evolution_of_delta[d][t] = y[0];
        this->cylindrical_collapse_evolution_of_delta_ddelta[d][t] = y[2];
        this->cylindrical_collapse_evolution_of_delta_ddelta2[d][t] = y[4];
      }
      else if(t > 0){
        this->cylindrical_collapse_evolution_of_delta[d][t] = constants::max_contrast;
        this->cylindrical_collapse_evolution_of_delta_ddelta[d][t] = 0.0;
        this->cylindrical_collapse_evolution_of_delta_ddelta2[d][t] = 0.0;
      }
    }
    gsl_odeiv2_driver_free(driver);
  }
  
}

vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_PDF_3D(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
  cout << "Computing CGF_data:\n";
  cout.flush();
  vector<vector<double> > CGF_data = this->compute_phi_of_lambda_3D(z, R_in_Mpc_over_h/constants::c_over_e5, f_NL, var_NL_rescale);
  cout << "Computing CGF_data done.\n";
  
  int n_lambda = 0;
  for(int i = 1; i < CGF_data[0].size(); i++){
    // ISSUE: this may fail for strong PNG, when there is two branch cuts in the CGF
    if(CGF_data[0][i-1] < CGF_data[0][i]){
      n_lambda = i+1;
    }
  }
  
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> d2phi_dlambda2(n_lambda, 0.0);
  for(int i = 0; i < n_lambda; i++){
    lambda_values[i] = CGF_data[0][i];
    d2phi_dlambda2[i] = CGF_data[3][i];
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = interpolate_neville_aitken(0.0, &lambda_values, &d2phi_dlambda2, constants::order_of_interpolation);
  double skew = interpolate_neville_aitken_derivative(0.0, &lambda_values, &d2phi_dlambda2, constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var, skew);
  double var_Gauss = log(1.0+var/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  double ddelta = (delta_max-delta_min)/double(n_delta);
  
  vector<double> delta_values(n_delta, 0.0);
  for(int i = 0; i < n_delta; i++){
    delta_values[i] = delta_min + double(i)*ddelta;
  }
  
  return this->compute_PDF_from_CGF(&delta_values, &CGF_data);
  
}


vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_PDF_2D(double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
  cout << "Computing CGF_data:\n";
  cout.flush();
  
  // var(L=L_input) \approx var(L=c/H_0)*(c/H_0/L_input)
  double L = L_in_Mpc_over_h/constants::c_over_e5;
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  double e = this->eta_at_a(1.0/(1.0+z));
  double D_growth = this->return_D_of_eta(e);
  
  this->current_P_NL = this->P_NL(e);
  double var_NL_Limber = variance_of_matter_within_R_NL_2D(R)/L;
  double var_NL_exact = variance_of_matter_within_R_NL_2D(R,L);
  
  vector<vector<double> > CGF_data = this->compute_phi_tilde_of_lambda_2D(e, R, f_NL, var_NL_rescale/L*var_NL_exact/var_NL_Limber);
  
  cout << "Computing CGF_data done.\n";
  
  int n_lambda = 0;
  for(int i = 1; i < CGF_data[0].size(); i++){
    // ISSUE: this may fail for strong PNG, when there is two branch cuts in the CGF
    if(CGF_data[0][i-1] < CGF_data[0][i]) n_lambda = i+1;
  }
  
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> d2phi_dlambda2(n_lambda, 0.0);
  for(int i = 0; i < n_lambda; i++){
    lambda_values[i] = CGF_data[0][i];
    d2phi_dlambda2[i] = CGF_data[3][i];
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = interpolate_neville_aitken(0.0, &lambda_values, &d2phi_dlambda2, constants::order_of_interpolation);
  double skew = interpolate_neville_aitken_derivative(0.0, &lambda_values, &d2phi_dlambda2, constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var, skew);
  double var_Gauss = log(1.0+var/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  double ddelta = (delta_max-delta_min)/double(n_delta);
  
  vector<double> delta_values(n_delta, 0.0);
  for(int i = 0; i < n_delta; i++){
    delta_values[i] = delta_min + double(i)*ddelta;
  }
  
  return this->compute_PDF_from_CGF(&delta_values, &CGF_data);
  
}





vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_PDF_at_chosen_deltas_2D(vector<double>* chosen_deltas, double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
  cout << "Computing phi_data:\n";
  cout.flush();
  
  // var(L=L_input) \approx var(L=c/H_0)*(c/H_0/L_input)
  double L = L_in_Mpc_over_h/constants::c_over_e5;
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  double e = this->eta_at_a(1.0/(1.0+z));
  double D_growth = this->return_D_of_eta(e);
  
  this->current_P_NL = this->P_NL(e);
  double var_NL_Limber = variance_of_matter_within_R_NL_2D(R)/L;
  double var_NL_exact = variance_of_matter_within_R_NL_2D(R,L);
  
  vector<vector<double> > CGF_data = this->compute_phi_tilde_of_lambda_2D(e, R, f_NL, var_NL_rescale/L*var_NL_exact/var_NL_Limber);
  cout << "Computing CGF_data done.\n";
  
  return this->compute_PDF_from_CGF(chosen_deltas, &CGF_data);
  
}



vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_PDF_from_CGF(vector<double>* delta_values, vector<vector<double> >* CGF_data){
  
  int n_delta = delta_values->size();
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  
  int n_lambda = 0;
  double lambda_c = (*CGF_data)[0][0];
  double delta_c = (*CGF_data)[2][0];
  double tau_c = 0.0;
  for(int i = 1; i < (*CGF_data)[0].size(); i++){
    // ISSUE: this may fail for strong PNG, when there are two branch cuts in the CGF
    if(((*CGF_data)[0][i-1] < (*CGF_data)[0][i])&&((*CGF_data)[0][i]*(*CGF_data)[2][i] > (*CGF_data)[1][i])){
      n_lambda = i+1;
      lambda_c = (*CGF_data)[0][i];
      delta_c = (*CGF_data)[2][i];
      tau_c = sqrt(2.0*((*CGF_data)[0][i]*(*CGF_data)[2][i] - (*CGF_data)[1][i]));
    }
    else if((*CGF_data)[0][i] > 0.0){
      i = 2*(*CGF_data)[0].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   * 0 = lambda
   * 1 = phi
   * 2 = phi'
   * 3 = phi''
   * 4 = phi'''
   * 5 = bias_term_1
   * 6 = bias_term_2
   * 7 = variance_term_1
   * 8 = variance_term_2
   * 9 = variance_term_3
   * 
   */
  vector<vector<double> > cut_CGF_data(10, vector<double>(n_lambda, 0.0));
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int c = 0; c < 10; c++){
    for(int i = 0; i < n_lambda; i++){
      cut_CGF_data[c][i] = (*CGF_data)[c][i];
    }
  }
  for(int i = 0; i < n_lambda; i++){
    tau_values[i] = sqrt(2.0*((*CGF_data)[0][i]*(*CGF_data)[2][i] - (*CGF_data)[1][i]));
    if(cut_CGF_data[0][i] < 0.0) tau_values[i] *= -1.0;
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  // The 0.9 is just to be save that we are sufficiently far away from the branch cuts.
  // ISSUE: This should also not be hardcoded!
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  tau_c = tau_max;
  
  delta_c = interpolate_neville_aitken(tau_max, &tau_values, &cut_CGF_data[2], constants::order_of_interpolation);
  lambda_c = interpolate_neville_aitken(tau_max, &tau_values, &cut_CGF_data[0], constants::order_of_interpolation);
  
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  // ISSUE: this shouldn't be hardcoded!
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  vector<double> bias_term_1_for_fit(n_tau,0.0);
  vector<double> bias_term_2_for_fit(n_tau,0.0);
  vector<double> variance_term_1_for_fit(n_tau,0.0);
  vector<double> variance_term_2_for_fit(n_tau,0.0);
  vector<double> variance_term_3_for_fit(n_tau,0.0);
  
  double dt = -tau_min/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i] = tau_min+double(i)*dt;
  }
  tau_for_fit[n_tau/2] = 0.0;
  dt = tau_max/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i+n_tau/2+1] = double(i+1)*dt;
  }
  
  for(int i = 0; i < n_tau; i++){
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &cut_CGF_data[2], constants::order_of_interpolation);
    bias_term_1_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &cut_CGF_data[5], constants::order_of_interpolation);
    bias_term_2_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &cut_CGF_data[6], constants::order_of_interpolation);
    variance_term_1_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &cut_CGF_data[7], constants::order_of_interpolation);
    variance_term_2_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &cut_CGF_data[8], constants::order_of_interpolation);
    variance_term_3_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &cut_CGF_data[9], constants::order_of_interpolation);
  }
  
  
  cout << "Done.\n";
  
  /*
   * Express functions as polynomials in tau.
   * 
   */
  cout << "Computing tau coefficients:\n";
  
  int n_coeff = constants::generating_function_coeff_order;
  
  vector<double> coefficients_phi_prime_of_tau = return_coefficients(&tau_for_fit, &phi_prime_for_fit, n_coeff);
  vector<double> coefficients_dphi_prime_dtau(coefficients_phi_prime_of_tau.size(), 0.0);
  vector<double> coefficients_d2phi_prime_dtau2(coefficients_phi_prime_of_tau.size(), 0.0);
  vector<double> coefficients_bias_term_1_of_tau = return_coefficients(&tau_for_fit, &bias_term_1_for_fit, n_coeff);
  vector<double> coefficients_bias_term_2_of_tau = return_coefficients(&tau_for_fit, &bias_term_2_for_fit, n_coeff);
  vector<double> coefficients_variance_term_1_of_tau = return_coefficients(&tau_for_fit, &variance_term_1_for_fit, n_coeff);
  vector<double> coefficients_variance_term_2_of_tau = return_coefficients(&tau_for_fit, &variance_term_2_for_fit, n_coeff);
  vector<double> coefficients_variance_term_3_of_tau = return_coefficients(&tau_for_fit, &variance_term_3_for_fit, n_coeff);
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  
  /*
   * ISSUE: one can enforce even more coefficients to their analytical value!
   */
  
  double var = interpolate_neville_aitken(0.0, &cut_CGF_data[0], &cut_CGF_data[3], constants::order_of_interpolation);
  coefficients_phi_prime_of_tau[0] = 0.0;
  coefficients_phi_prime_of_tau[1] = sqrt(var);
  coefficients_bias_term_1_of_tau[0] = 0.0;
  coefficients_bias_term_2_of_tau[0] = 0.0;
  coefficients_bias_term_2_of_tau[1] = 0.0;
  
  for(int i = 0; i < coefficients_phi_prime_of_tau.size()-1; i++){
    coefficients_dphi_prime_dtau[i] = coefficients_phi_prime_of_tau[i+1]*double(i+1);
  }
  for(int i = 0; i < coefficients_phi_prime_of_tau.size()-1; i++){
    coefficients_d2phi_prime_dtau2[i] = coefficients_dphi_prime_dtau[i+1]*double(i+1);
  }
  
  double delta;
  double tau_0, lambda_0;
  double dr;
  
  complex<double> lambda;
  complex<double> lambda_next;
  complex<double> tau, tau_next;
  complex<double> phi_prime, phi_prime_next;
  complex<double> b_term_1, b_term_1_next;
  complex<double> b_term_2, b_term_2_next;
  complex<double> v_term_1, v_term_1_next;
  complex<double> v_term_2, v_term_2_next;
  complex<double> v_term_3, v_term_3_next;
  complex<double> exponent, exponent_next;
  complex<double> dlambda;
  complex<double> step, first_step;
  
  vector<vector<double> > PDF_data(7, vector<double>(n_delta));
  
  cout << "Computing PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = (*delta_values)[i];
    
    PDF_data[0][i] = delta;
    PDF_data[1][i] = 0.0;
      
    if(delta < delta_c){
      tau_0 = interpolate_Newton(delta, &cut_CGF_data[2], &tau_values, constants::order_of_interpolation);
      lambda_0 = interpolate_Newton(delta, &cut_CGF_data[2], &cut_CGF_data[0], constants::order_of_interpolation);
    }
    else{
      tau_0 = tau_c;
      lambda_0 = lambda_c;
    }
    lambda = complex<double>(lambda_0, 0.0);
    tau = complex<double>(tau_0, 0.0);
    exponent = exp(-0.5*pow(tau, 2));
    phi_prime = return_polnomial_value(tau, &coefficients_phi_prime_of_tau);
    b_term_1 = return_polnomial_value(tau, &coefficients_bias_term_1_of_tau);
    b_term_2 = return_polnomial_value(tau, &coefficients_bias_term_2_of_tau);
    v_term_1 = return_polnomial_value(tau, &coefficients_variance_term_1_of_tau);
    v_term_2 = return_polnomial_value(tau, &coefficients_variance_term_2_of_tau);
    v_term_3 = return_polnomial_value(tau, &coefficients_variance_term_3_of_tau);
    
    // sigma_r^2 \approx 1/phi''(lambda_0)
    double sigma_frac = 0.001;
    dr = sigma_frac/sqrt(interpolate_neville_aitken(lambda_0, &cut_CGF_data[0], &cut_CGF_data[3], constants::order_of_interpolation));
    dlambda = complex<double>(0.0, dr);
    int j = 0;
    do{
      lambda_next = lambda + 0.5*dlambda;
      tau_next = get_tau_from_secant_method_complex_Bernardeau_notation_2D(lambda_next, tau, &coefficients_dphi_prime_dtau, &coefficients_d2phi_prime_dtau2);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_prime_of_tau);
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda_next = lambda + dlambda;
      tau_next = get_tau_from_secant_method_complex_Bernardeau_notation_2D(lambda_next, tau_next, &coefficients_dphi_prime_dtau, &coefficients_d2phi_prime_dtau2);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_prime_of_tau);
      exponent_next = exp(-lambda_next*(delta-phi_prime_next)-0.5*pow(tau_next, 2));
      b_term_1_next = return_polnomial_value(tau_next, &coefficients_bias_term_1_of_tau);
      b_term_2_next = return_polnomial_value(tau_next, &coefficients_bias_term_2_of_tau);
      v_term_1_next = return_polnomial_value(tau_next, &coefficients_variance_term_1_of_tau);
      v_term_2_next = return_polnomial_value(tau_next, &coefficients_variance_term_2_of_tau);
      v_term_3_next = return_polnomial_value(tau_next, &coefficients_variance_term_3_of_tau);
      
      step = 0.5*dlambda*(exponent_next*b_term_1_next+exponent*b_term_1);
      PDF_data[2][i] += step.imag();
      
      step = 0.5*dlambda*(exponent_next*b_term_2_next+exponent*b_term_2);
      PDF_data[3][i] += step.imag();
      
      step = 0.5*dlambda*(exponent_next*v_term_1_next+exponent*v_term_1);
      PDF_data[4][i] += step.imag();
      
      step = 0.5*dlambda*(exponent_next*v_term_2_next+exponent*v_term_2);
      PDF_data[5][i] += step.imag();
      
      step = 0.5*dlambda*(exponent_next*v_term_3_next+exponent*v_term_3);
      PDF_data[6][i] += step.imag();
      
      step = 0.5*dlambda*(exponent_next+exponent);
      PDF_data[1][i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      phi_prime = phi_prime_next;
      exponent = exponent_next;
      b_term_1 = b_term_1_next;
      b_term_2 = b_term_2_next;
      v_term_1 = v_term_1_next;
      v_term_2 = v_term_2_next;
      v_term_3 = v_term_3_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF_data[1][i] /= constants::pi;
    PDF_data[2][i] /= constants::pi*PDF_data[1][i];
    PDF_data[3][i] /= constants::pi*PDF_data[1][i];
    PDF_data[4][i] /= constants::pi*PDF_data[1][i];
    PDF_data[5][i] /= constants::pi*PDF_data[1][i];
    PDF_data[6][i] /= constants::pi*PDF_data[1][i];
    
    cout << PDF_data[0][i] << "   ";
    cout << PDF_data[1][i] << "\n";
    
  }
  cout << "Done.\n";
  
  return PDF_data;
  
}


vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_phi_of_lambda_3D(double z, double R, double f_NL, double var_NL_rescale){
  
  double e = this->eta_at_a(1.0/(1.0+z));
  double log_R = log(R);
  double RL, log_RL;
  
  this->current_P_NL = this->P_NL(e);
  this->current_P_L = this->P_L(e);
  double D_growth = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double var_NL_R = variance_of_matter_within_R_NL(R)*var_NL_rescale;
  double var_L_R = variance_of_matter_within_R(R);
  double var_L_RL, skew_L_RL, dvar_L_RL_dR, dskew_L_RL_dR;
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  vector<double> delta_NL_prime_prime_values;
  
  this->return_delta_NL_of_delta_L_and_dF_ddelta_3D(e, &delta_L_values, &delta_NL_values, &delta_NL_prime_values, &delta_NL_prime_prime_values);
  
  vector<vector<double> > data(10, vector<double>(delta_L_values.size(), 0.0));
  
  vector<double> j_values = delta_L_values;
  vector<double> lambda_values = delta_L_values;
  vector<double> dlambda_ddelta = delta_L_values;
  vector<double> phi_values = delta_L_values;
  
  int N_radii = this->log_top_hat_radii.size();
  double log_R_min = this->log_top_hat_radii[0];
  double log_R_max = this->log_top_hat_radii[N_radii-1];
  
  double var_L_R_min = this->top_hat_sphere_variances[0];
  double var_L_R_max = this->top_hat_sphere_variances[N_radii-1];
  
  double dvar_L_R_min_dR = this->dtop_hat_sphere_variances_dR[0];
  double dvar_L_R_max_dR = this->dtop_hat_sphere_variances_dR[N_radii-1];
  
  // auxilliary variables to store intermediate results:
  double dphi_dlnR;
  double d2phi_dlnR2, d2phi_dlnRdj, d2phi_dj2;
  double dvar_dlnR, d2var_dlnR2;
  double dskew_dlnR, d2skew_dlnR2; // skew is 3rd central moment
  double dF_dlambda;
  double one_over_one_plus_F;
  double delta_L_at_z;
  
  int n_delta = delta_L_values.size();
  int i_min = 0;
  int i_max = n_delta-1;
  double R_min = exp(this->log_top_hat_radii[0]);
  double R_max = exp(this->log_top_hat_radii[this->log_top_hat_radii.size()-1]);
  if(f_NL != 0.0){
    R_min = max(R_min, exp(this->log_top_hat_radii_for_skewnesses[0]));
    R_max = min(R_max, exp(this->log_top_hat_radii_for_skewnesses[this->log_top_hat_radii_for_skewnesses.size()-1]));
  }
  
  for(int d = 0; d < n_delta; d++){
    delta_L_at_z = delta_L_values[d]*D_growth;
    log_RL = log_R + log(1.0+delta_NL_values[d])/3.0;
    RL = exp(log_RL);
    if(RL < R_min) i_min++;
    if(RL > R_max) i_max--;
    
    var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
    dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
    skew_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->top_hat_sphere_skewnesses, constants::order_of_interpolation);
    dskew_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->dtop_hat_sphere_skewnesses_dR, constants::order_of_interpolation);
    
    skew_L_RL *= f_NL;
    dskew_L_RL_dR *= f_NL;
    
    if(f_NL != 0.0){
      
      if(this->skewness_initialisation == UNINITIALISED){
        cerr << "ERROR: skewness_initialisation == UNINITIALISED , i.e.\n";
        cerr << "no primordial skewnesses have been computed.\n";
        cerr << "either set f_NL = 0.0 or initialise skewnesses.\n";
        exit(1);
      }
      
      /*
       * phi(j) = var/2 * j**2 + skew/6*j**3
       * phi'(j) = var * j + skew/2*j**2
       * => 0 = j**2 + 2 var/skew * j - 2/skew * delta
       * => j = -var/skew \pm sqrt(var**2/skew**2 + 2/skew * delta)
       */
      
      j_values[d] = -var_L_RL/skew_L_RL;
      if(skew_L_RL > 0.0)
        j_values[d] += sqrt(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL);
      else
        j_values[d] -= sqrt(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL);
      
      lambda_values[d] = j_values[d]/delta_NL_prime_values[d];
      lambda_values[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*(dvar_L_RL_dR/2.0*pow(j_values[d],2) + dskew_L_RL_dR/6.0*pow(j_values[d],3));
      
      phi_values[d] = -lambda_values[d]*delta_NL_values[d]+delta_L_values[d]*j_values[d] - (var_L_RL/2.0*pow(j_values[d],2) + skew_L_RL/6.0*pow(j_values[d],3));
      phi_values[d] *= -1.0;
    }
    else{
      j_values[d] = delta_L_values[d]/var_L_RL;
      lambda_values[d] = j_values[d]/delta_NL_prime_values[d];
      lambda_values[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values[d],2);
      phi_values[d] = lambda_values[d]*delta_NL_values[d]-delta_L_values[d]*j_values[d] + var_L_RL/2.0*pow(j_values[d],2);
    }
    
    // now calculating dF/dlambda (needed in same applications for saddle point approxiamtion)
    dvar_dlnR = dvar_L_RL_dR*RL;
    d2var_dlnR2 = interpolate_neville_aitken_derivative(log_RL, &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
    d2var_dlnR2 = RL*d2var_dlnR2 + dvar_dlnR;
    dphi_dlnR = 0.5*pow(j_values[d],2)*dvar_dlnR;
    d2phi_dlnR2 = 0.5*pow(j_values[d],2)*d2var_dlnR2;
    d2phi_dlnRdj = j_values[d]*dvar_dlnR;
    d2phi_dj2 = var_L_RL;
    
    if(f_NL != 0.0){
      dskew_dlnR = dskew_L_RL_dR*RL;
      // this is d(dskew_dR)/dlnR:
      d2skew_dlnR2 = interpolate_neville_aitken_derivative(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->dtop_hat_sphere_skewnesses_dR, constants::order_of_interpolation);
      d2skew_dlnR2 = RL*d2skew_dlnR2 + dskew_dlnR;
      dphi_dlnR += pow(j_values[d],3)/6.0*dskew_dlnR;
      d2phi_dlnR2 += pow(j_values[d],3)/6.0*d2skew_dlnR2;
      d2phi_dlnRdj += pow(j_values[d],2)*0.5*dskew_dlnR;
      d2phi_dj2 += j_values[d]*skew_L_RL;
    }
    one_over_one_plus_F = 1.0/(1.0+delta_NL_values[d]);
    dlambda_ddelta[d] = pow(1.0 - one_over_one_plus_F/3.0*delta_NL_prime_values[d]*d2phi_dlnRdj, 2);
    dlambda_ddelta[d] /= delta_NL_prime_values[d]*d2phi_dj2;
    dlambda_ddelta[d] += pow(one_over_one_plus_F/3.0,2)*delta_NL_prime_values[d]*(3.0*dphi_dlnR - d2phi_dlnR2);
    dlambda_ddelta[d] += -delta_NL_prime_prime_values[d]*j_values[d]/pow(delta_NL_prime_values[d], 2);
    dF_dlambda = delta_NL_prime_values[d]/dlambda_ddelta[d];
    // finished calculating dF/dlambda
    
    phi_values[d] *= var_L_R/var_NL_R;
    lambda_values[d] *= var_L_R/var_NL_R;
    dlambda_ddelta[d] *= var_L_R/var_NL_R;
    dF_dlambda /= var_L_R/var_NL_R;
    
    data[0][d] = lambda_values[d];
    data[1][d] = phi_values[d];
    data[2][d] = delta_NL_values[d]; // dphi/dlambda
    data[3][d] = dF_dlambda;         // d^2phi/dlambda^2
    //                               // d^3phi/dlambda^3 see below
    /* 
     * bias terms, needed when calculating
     * <delta_g | delta_m> = delta_m + b_1^L bias_term_1 + b_2^L bias_term_2
     * 
     */
    data[5][d] = (1.0+data[2][d])*delta_L_at_z; // bias term 1
    data[6][d] = data[4][d]*delta_L_at_z;       // bias term 2
    /* 
     * variance terms, needed when calculating
     * <delta_g^2 | delta_m> = variance_term_1 + b_1^L variance_term_2 + (b_1^L)^2 variance_term_3
     * 
     */
    data[7][d] = dF_dlambda + pow(delta_NL_values[d],2);                                                                                                                              // variance term 1
    data[8][d] = 2.0*(dF_dlambda*(1.0+delta_NL_values[d]+delta_NL_prime_values[d]*delta_L_at_z)/delta_NL_prime_values[d] + delta_NL_values[d]*(1.0+delta_NL_values[d])*delta_L_at_z); // variance term 2
    data[9][d] = dF_dlambda*pow((1.0+delta_NL_values[d]+delta_NL_prime_values[d]*delta_L_at_z)/delta_NL_prime_values[d],2) + pow((1.0+delta_NL_values[d])*delta_L_at_z,2);            // variance term 3
    
  }
  
  //for(int d = 0; d < n_delta; d++){
  //  data[3][d] = interpolate_neville_aitken_derivative(data[0][d], &data[0], &data[2], constants::order_of_interpolation);
  //}
  for(int d = 0; d < n_delta; d++){
    data[4][d] = interpolate_neville_aitken_derivative(data[0][d], &data[0], &data[3], constants::order_of_interpolation);
  }
  
  // making sure the CGF data is returned over a range for R, for which the linear moments had been computed.
  int n_cut = i_max+1-i_min;
  
  if(n_cut == n_delta)
    return data;
  
  vector<vector<double> > data_cut(10, vector<double>(n_cut, 0.0));
  
  for(int c = 0; c < 10; c++){
    for(int d = 0; d < n_cut; d++){
      data_cut[c][d] = data[c][i_min+d];
    }
  }
  
  return data_cut;
  
}

vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_phi_tilde_of_lambda_2D(double e, double R, double f_NL, double var_NL_rescale){
  
  double log_R = log(R);
  double RL, log_RL;
  double D_growth = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
  this->current_P_L = this->P_L(e);
  this->current_P_NL = this->P_NL(e);
  
  double var_NL_R = variance_of_matter_within_R_NL_2D(R)*var_NL_rescale;
  double var_L_R = variance_of_matter_within_R_2D(R);
  double var_L_RL, skew_L_RL, dvar_L_RL_dR, dskew_L_RL_dR;
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  vector<double> delta_NL_prime_prime_values;
  
  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(e, &delta_L_values, &delta_NL_values, &delta_NL_prime_values, &delta_NL_prime_prime_values);
  
  vector<vector<double> > data(10, vector<double>(delta_L_values.size(), 0.0));
  vector<double> lambda_of_delta = delta_L_values;
  vector<double> dlambda_ddelta = delta_L_values;
  vector<double> phi_of_delta = delta_L_values;
  vector<double> j_values = delta_L_values;
  vector<double> j_values_Gauss = delta_L_values;
  vector<double> delta_g_prime_values = delta_L_values;
  
  // auxilliary variables to store intermediate results:
  double dphi_dlnR;
  double d2phi_dlnR2, d2phi_dlnRdj, d2phi_dj2;
  double dvar_dlnR, d2var_dlnR2;
  double dskew_dlnR, d2skew_dlnR2; // skew is 3rd central moment
  double dF_dlambda;
  double one_over_one_plus_F;
  double delta_L_at_z;
  
  int n_delta = delta_L_values.size();
  int i_min = 0;
  int i_max = n_delta-1;
  double R_min = exp(this->log_top_hat_cylinder_radii[0]);
  double R_max = exp(this->log_top_hat_cylinder_radii[this->log_top_hat_cylinder_radii.size()-1]);
  if(f_NL != 0.0){
    R_min = max(R_min, exp(this->log_top_hat_cylinder_radii_for_skewnesses[0]));
    R_max = min(R_max, exp(this->log_top_hat_cylinder_radii_for_skewnesses[this->log_top_hat_cylinder_radii_for_skewnesses.size()-1]));
  }
  
  for(int d = 0; d < delta_L_values.size(); d++){
    delta_L_at_z = delta_L_values[d]*D_growth;
    log_RL = log_R + log(1.0+delta_NL_values[d])/2.0;
    RL = exp(log_RL);
    if(RL < R_min) i_min++;
    if(RL > R_max) i_max--;
    
    var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii, &this->top_hat_cylinder_variances, constants::order_of_interpolation);
    dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii, &this->dtop_hat_cylinder_variances_dR, constants::order_of_interpolation);
    
    if(f_NL != 0.0){
      
      skew_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii_for_skewnesses, &this->top_hat_cylinder_skewnesses, constants::order_of_interpolation);
      dskew_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii_for_skewnesses, &this->dtop_hat_cylinder_skewnesses_dR, constants::order_of_interpolation);
    
      skew_L_RL *= f_NL;
      dskew_L_RL_dR *= f_NL;
      
      if(this->skewness_initialisation_2D == UNINITIALISED){
        cerr << "ERROR: skewness_initialisation == UNINITIALISED , i.e.\n";
        cerr << "no primordial skewnesses have been computed.\n";
        cerr << "either set f_NL = 0.0 or initialise skewnesses.\n";
        exit(1);
      }
      
      /*
       * phi(j) = var/2 * j**2 + skew/6*j**3
       * phi'(j) = var * j + skew/2*j**2
       * => 0 = j**2 + 2 var/skew * j - 2/skew * delta
       * => j = -var/skew \pm sqrt(var**2/skew**2 + 2/skew * delta)
       */
      
      if(skew_L_RL > 0.0){
        j_values[d] = -var_L_RL/skew_L_RL;
        j_values[d] += sqrt(max(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL, 0.0));
      }
      else if(skew_L_RL < 0.0){
        j_values[d] = -var_L_RL/skew_L_RL;
        j_values[d] -= sqrt(max(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL, 0.0));
      }
      else
        j_values[d] = delta_L_values[d]/var_L_RL;
      
      lambda_of_delta[d] = j_values[d]/delta_NL_prime_values[d];
      lambda_of_delta[d] -= RL/(2.0*(1.0+delta_NL_values[d]))*(dvar_L_RL_dR/2.0*pow(j_values[d],2) + dskew_L_RL_dR/6.0*pow(j_values[d],3));
      
      phi_of_delta[d] = -lambda_of_delta[d]*delta_NL_values[d]+delta_L_values[d]*j_values[d] - (var_L_RL/2.0*pow(j_values[d],2) + skew_L_RL/6.0*pow(j_values[d],3));
      phi_of_delta[d] *= -1.0;
    }
    else{
      j_values[d] = delta_L_values[d]/var_L_RL;
      lambda_of_delta[d] = j_values[d]/delta_NL_prime_values[d];
      lambda_of_delta[d] -= RL/(2.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values[d],2);
      phi_of_delta[d] = lambda_of_delta[d]*delta_NL_values[d]-delta_L_values[d]*j_values[d] + var_L_RL/2.0*pow(j_values[d],2);
    }
    
    
    // now calculating dF/dlambda (needed in same applications for saddle point approxiamtion)
    dvar_dlnR = dvar_L_RL_dR*RL;
    d2var_dlnR2 = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii, &this->d2top_hat_cylinder_variances_dR2, constants::order_of_interpolation);
    d2var_dlnR2 = pow(RL, 2)*d2var_dlnR2 + dvar_dlnR;
    dphi_dlnR = 0.5*pow(j_values[d],2)*dvar_dlnR;
    d2phi_dlnR2 = 0.5*pow(j_values[d],2)*d2var_dlnR2;
    d2phi_dlnRdj = j_values[d]*dvar_dlnR;
    d2phi_dj2 = var_L_RL;
    
    if(f_NL != 0.0){
      dskew_dlnR = dskew_L_RL_dR*RL;
      // this is d(dskew_dR)/dlnR:
      d2skew_dlnR2 = interpolate_neville_aitken_derivative(log_RL, &this->log_top_hat_cylinder_radii_for_skewnesses, &this->dtop_hat_cylinder_skewnesses_dR, constants::order_of_interpolation);
      d2skew_dlnR2 = RL*d2skew_dlnR2 + dskew_dlnR;
      dphi_dlnR += pow(j_values[d],3)/6.0*dskew_dlnR;
      d2phi_dlnR2 += pow(j_values[d],3)/6.0*d2skew_dlnR2;
      d2phi_dlnRdj += pow(j_values[d],2)*0.5*dskew_dlnR;
      d2phi_dj2 += j_values[d]*skew_L_RL;
    }
    one_over_one_plus_F = 1.0/(1.0+delta_NL_values[d]);
    dlambda_ddelta[d] = pow(1.0 - 0.5*one_over_one_plus_F*delta_NL_prime_values[d]*d2phi_dlnRdj, 2);
    dlambda_ddelta[d] /= delta_NL_prime_values[d]*d2phi_dj2;
    dlambda_ddelta[d] += 0.25*pow(one_over_one_plus_F,2)*delta_NL_prime_values[d]*(2.0*dphi_dlnR - d2phi_dlnR2);
    dlambda_ddelta[d] += -delta_NL_prime_prime_values[d]*j_values[d]/pow(delta_NL_prime_values[d], 2);
    dF_dlambda = delta_NL_prime_values[d]/dlambda_ddelta[d];
    // finished calculating dF/dlambda
    
    phi_of_delta[d] *= var_L_R/var_NL_R;
    lambda_of_delta[d] *= var_L_R/var_NL_R;
    dlambda_ddelta[d] *= var_L_R/var_NL_R;
    dF_dlambda /= var_L_R/var_NL_R;
    
    data[0][d] = lambda_of_delta[d];
    data[1][d] = phi_of_delta[d];
    data[2][d] = delta_NL_values[d]; // dphi/dlambda
    data[3][d] = dF_dlambda;         // d^2phi/dlambda^2
    //                               // d^3phi/dlambda^3 see below
    /* 
     * bias terms, needed when calculating
     * <delta_g | delta_m> = delta_m + b_1^L bias_term_1 + b_2^L bias_term_2
     * 
     */
    data[5][d] = (1.0+data[2][d])*delta_L_at_z; // bias term 1
    data[6][d] = data[4][d]*delta_L_at_z;       // bias term 2
    /* 
     * variance terms, needed when calculating
     * <delta_g^2 | delta_m> = variance_term_1 + b_1^L variance_term_2 + (b_1^L)^2 variance_term_3
     * 
     */
    data[7][d] = dF_dlambda + pow(delta_NL_values[d],2);                                                                                                                              // variance term 1
    data[8][d] = 2.0*(dF_dlambda*(1.0+delta_NL_values[d]+delta_NL_prime_values[d]*delta_L_at_z)/delta_NL_prime_values[d] + delta_NL_values[d]*(1.0+delta_NL_values[d])*delta_L_at_z); // variance term 2
    data[9][d] = dF_dlambda*pow((1.0+delta_NL_values[d]+delta_NL_prime_values[d]*delta_L_at_z)/delta_NL_prime_values[d],2) + pow((1.0+delta_NL_values[d])*delta_L_at_z,2);            // variance term 3
    
  }
  
  for(int d = 0; d < n_delta; d++){
    data[4][d] = interpolate_neville_aitken_derivative(data[0][d], &data[0], &data[3], constants::order_of_interpolation);
  }
  
  // making sure the CGF data is returned over a range for R, for which the linear moments had been computed.
  int n_cut = i_max+1-i_min;
  
  if(n_cut == n_delta)
    return data;
  
  vector<vector<double> > data_cut(10, vector<double>(n_cut, 0.0));
  
  
  for(int c = 0; c < 10; c++){
    for(int d = 0; d < n_cut; d++){
      data_cut[c][d] = data[c][i_min+d];
    }
  }
  
  return data_cut;
  
}









