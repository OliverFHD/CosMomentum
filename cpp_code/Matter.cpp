
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_legendre.h>

#include "delta0_and_kappa0.h"

#include "Matter.h"
#include "Eisenstein_and_Hu.h"
#include "generating_function_utils.h"
#include "Matter_integrands.cpp"
#include "Matter_configuration.cpp"
#include "Matter_initial_conditions.cpp"
#include "Matter_analytic_special_cases.cpp"
#include "Matter_halofit.cpp"
#include "Matter_variances.cpp"
#include "Matter_output.cpp"


/* _____ Content Matter.cpp
 * 
 * _____ 1. Initialization
 * 
 * ..... 1.1 Constructor
 * ..... 1.2 Destructor
 * ..... 1.3 set_wave_numbers
 * ..... 1.4 set_transfer_function_Eisenstein_and_Hu
 * ..... 1.5 change_cosmology (1)
 * ..... 1.6 change_cosmology (2)
 * 
 * _____ 3. Growth Factors and Power Spectra
 * 
 * ..... 3.2 initialize_linear_growth_factor_of_delta
 * ..... 3.5 Newtonian_linear_power_spectrum
 * ..... 3.8 variance_of_matter_within_R
 * 
 * _____ 5. Output and Checks
 * 
 * ..... 5.2 print_Newtonian_growth_factor
 * ..... 5.6 return_wave_numbers
 * ..... 5.7 transfer_function_at
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
 * ..... 1.3 set_wave_numbers              ****************
 * ..... 1.4 set_transfer_function_Eisenstein_and_Hu      *
 * ..... 1.5 set_transfer_function_Bond_and_Efstathiou    *
 * ..... 1.6 set_cylinder_variances        ****************
 *                                         *
 *******************************************
 *******************************************/

/*******************************************************************************************************************************************************
 * 1.1
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

Matter::Matter(Universe* uni){
  
  this->universe = uni;
  this->cosmology = this->universe->return_cosmology(); 
  
  this->a_initial = this->universe->return_a_initial();
  this->a_final = this->universe->return_a_final();
  this->eta_initial = this->universe->return_eta_initial();
  this->eta_final = this->universe->return_eta_final();
  
  this->set_initial_conditions();
  this->initialize_linear_growth_factor_of_delta();
  this->initialize_up_to_second_order_growth_factor_of_delta(2.0/7.0*this->D_initial*this->D_initial, 4.0/7.0*this->D_initial*this->D_prime_initial);
  
  double z_of_f_NL_rescale = 500.0;
  double a_of_f_NL_rescale = 1.0/(1.0+z_of_f_NL_rescale);
  double eta_of_f_NL_rescale = this->universe->eta_at_a(a_of_f_NL_rescale);
  this->f_NL_rescaling_factor = this->return_D_of_eta(eta_of_f_NL_rescale)/a_of_f_NL_rescale;
  cout << f_NL_rescaling_factor << "   HAHAHA\n";
  cout << this->D_initial/this->a_initial << "   HAHAHA\n";
  cout << this->a_initial << "   HAHAHA\n";
  
  this->set_wave_numbers();
  if(this->cosmology.Omega_b != 0.0){
    //this->set_transfer_function_Eisenstein_and_Hu();
    //this->set_transfer_function_from_file("Coras_PDFs/input_transfer_manera.dat");
    //this->set_transfer_function_from_file("Coras_PDFs/tkwmap5.dat");
    //this->set_transfer_function_from_file("Coras_PDFs/tk_Quijote.dat");
    //this->set_transfer_function_from_file("Coras_PDFs/tk_Quijote_z0.dat");
    //this->set_transfer_function_from_file("../Tobias_data/transfer_function_my_format.dat");
    
    // NOTE: THIS IS BETTER FOR FLEXIBILITY IN THE PYTHON INTERFACE!
    this->set_transfer_function_from_file("./transfer_function.dat");
  }
  else{
    this->set_transfer_function_Bond_and_Efstathiou();
  }
  this->norm = this->variance_of_matter_within_R_before_norm_was_determined(8.0);
  this->set_P_today();
  
  //cout << "Setting cylinder variances.\n";
  //this->set_cylinder_variances();
  cout << "Setting sphere variances.\n";
  this->set_sphere_variances();
  cout << "Setting sphere skewnesses.\n";
  //this->set_sphere_skewnesses();
  
  cout << "Setting sphere skewnesses from file.\n";
  //this->set_sphere_skewnesses_from_file("../data_for_PNG_paper/Final_plots/primordial_moments_Nishimishi_overwrite_save.dat");
  this->set_sphere_skewnesses_from_file("../data_for_PNG_paper/Final_plots/primordial_moments_Oriana_overwrite_save.dat");
  //this->set_sphere_skewnesses_from_file("../data_for_PNG_paper/Final_plots/primordial_moments_high_res.dat");
  
  cout << "Setting spherical collapse evolution.\n";
  this->set_spherical_collapse_evolution_of_delta(0.0, 3.0, 1000);
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

Matter::~Matter(){
}


/*******************************************************************************************************************************************************
 * 1.3 set_wave_numbers
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_wave_numbers(){
  
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

void Matter::set_transfer_function_Eisenstein_and_Hu(){
  
  tranfer_function_Eisenstein_and_Hu(&this->wave_numbers, &this->transfer_function, this->cosmology.Omega_m, this->cosmology.Omega_b, this->cosmology.h_100, this->cosmology.theta_27);
  
}

/*******************************************************************************************************************************************************
 * 1.5 set_transfer_function_Eisenstein_and_Hu
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_transfer_function_Bond_and_Efstathiou(){
  
  tranfer_function_Bond_and_Efstathiou(&this->wave_numbers, &this->transfer_function, this->cosmology.Omega_m, this->cosmology.h_100);
  
}

/*******************************************************************************************************************************************************
 * 1.6 set_transfer_function_from_file
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_transfer_function_from_file(string file){
  
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
      //cout << k/constants::c_over_e5;
      //cout << '\t';
      //cout << Tk/Tk_0 << '\n';
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
 * ..... 3.5 variance_of_matter_within_R_2D                  *
 * ..... 3.6 variance_of_matter_within_R_2D_NL               *
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
 * - In this method there are still NO radiation inhomogeneities includes. This would also require to take into account
 *   entropy perturbation!
 * 
*******************************************************************************************************************************************************/

void Matter::initialize_linear_growth_factor_of_delta(){

  integration_parameters params;
  params.pointer_to_Matter = this;
  params.pointer_to_Universe = this->universe;
  params.Omega_m = this->cosmology.Omega_m;
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {growth_factor_gsl, growth_factor_gsl_jac, 2, (void *) pointer_to_params};

  this->number_of_entries_Newton = 100000;
  
  this->eta_Newton.resize(this->number_of_entries_Newton);
  this->Newtonian_growth_factor_of_delta.resize(this->number_of_entries_Newton);
  this->Newtonian_growth_factor_of_delta_prime.resize(this->number_of_entries_Newton);
  this->Newtonian_growth_factor_of_delta_prime_prime.resize(this->number_of_entries_Newton);
  this->Newtonian_growth_factor_of_delta_prime_prime_prime.resize(this->number_of_entries_Newton);
  
  this->eta_Newton[0] = this->eta_initial;
  this->Newtonian_growth_factor_of_delta[0] = this->D_initial;
  this->Newtonian_growth_factor_of_delta_prime[0] = this->D_prime_initial;
  
  
  
  double e;
  double w;
  double de = (this->eta_final - this->eta_initial)/double(this->number_of_entries_Newton - 1);
  double e_i = this->eta_initial;
  double e_f = this->eta_initial + de;
  double D0;
  double y[2] = { this->Newtonian_growth_factor_of_delta[0], this->Newtonian_growth_factor_of_delta_prime[0]};
  
  // All quantities are of order 1.0 at their maximum.
  double hstart = 1.0*constants::gsl_hstart_relative;
  double eps_absolute = 1.0*constants::gsl_eps_relative;
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, eps_absolute, constants::gsl_eps_relative);

  for (int i = 1; i < this->number_of_entries_Newton; i++){
    
    
    this->eta_Newton[i] = e_f;
    int status = gsl_odeiv2_driver_apply(d, &e_i, e_f, y);
    
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    
    
    double scale = this->universe->a_at_eta(e_f);
    double H = this->universe->H_at_eta(e_f);
    double H_prime = this->universe->H_prime_at_eta(e_f);
    double Om_m = this->cosmology.Omega_m;
    
    this->Newtonian_growth_factor_of_delta[i] = y[0];
    this->Newtonian_growth_factor_of_delta_prime[i] = y[1];
    this->Newtonian_growth_factor_of_delta_prime_prime[i] = -H*y[1]+3.0/2.0*Om_m/scale*y[0];
    this->Newtonian_growth_factor_of_delta_prime_prime_prime[i] = -H_prime*y[1]-H*this->Newtonian_growth_factor_of_delta_prime_prime[i]-3.0/2.0*H*Om_m/scale*y[0]+3.0/2.0*Om_m/scale*y[1];
    e_f += de;
  }

  gsl_odeiv2_driver_free(d);
  
  e = this->universe->eta_at_a(1.0);
  
  D0 = interpolate_neville_aitken(e, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
  this->D_initial /= D0;
  this->D_prime_initial /= D0;
    
  for (int i = 0; i < this->number_of_entries_Newton; i++){
    this->Newtonian_growth_factor_of_delta[i] /= D0;
    this->Newtonian_growth_factor_of_delta_prime[i] /= D0;
    this->Newtonian_growth_factor_of_delta_prime_prime[i] /= D0;
    this->Newtonian_growth_factor_of_delta_prime_prime_prime[i] /= D0;
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

void Matter::initialize_up_to_second_order_growth_factor_of_delta(double D, double D_prime){

  integration_parameters params;
  params.pointer_to_Matter = this;
  params.pointer_to_Universe = this->universe;
  params.Omega_m = this->cosmology.Omega_m;
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {growth_factor_to_second_order_gsl, growth_factor_to_second_order_gsl_jac, 2, (void *) pointer_to_params};

  
  this->Newtonian_growth_factor_second_order.resize(this->number_of_entries_Newton);
  this->Newtonian_growth_factor_second_order_prime.resize(this->number_of_entries_Newton);
  
  this->Newtonian_growth_factor_second_order[0] = D;
  this->Newtonian_growth_factor_second_order_prime[0] = D_prime;
  
  
  
  double e;
  double de = (this->eta_final - this->eta_initial)/double(this->number_of_entries_Newton - 1);
  double e_i = this->eta_Newton[0];
  double e_f;
  double y[2] = { this->Newtonian_growth_factor_second_order[0], this->Newtonian_growth_factor_second_order_prime[0]};
    
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

  for (int i = 1; i < this->number_of_entries_Newton; i++){

    e_f = this->eta_Newton[i];
    int status = gsl_odeiv2_driver_apply(d, &e_i, e_f, y);
    
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    
    this->Newtonian_growth_factor_second_order[i] = y[0];
    this->Newtonian_growth_factor_second_order_prime[i] = y[1];
    
    e_f += de;
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

double Matter::Newtonian_linear_power_spectrum(double k, double e){

  
  double D = interpolate_neville_aitken(e, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double transfer_sq = this->transfer_function_at(k); transfer_sq *= transfer_sq;
  return D*D*pow(k, this->cosmology.n_s)*transfer_sq*pow(this->cosmology.sigma_8, 2)/this->norm;
  //double ln_k = log(k);
  //return D*D*interpolate_neville_aitken(ln_k, &this->CAMB_ln_k_values, &this->CAMB_power_spectrum, constants::order_of_interpolation);
    
}





/*******************************************************************************************************************************************************
 * 5.12 set_cylindrical_collapse_evolution_of_delta
 * Description:
 * - for a number of initial values for delta, this function computes their evolution in the
 *   cylindrical collapse model
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::set_cylindrical_collapse_evolution_of_delta(double z_min, double z_max, int n_time){

  double delta_min = -5.0;
  double delta_max = 1.4;
  double ddelta = 0.02;
  double delta = delta_min-ddelta;
  double eta_i = this->eta_Newton[0];
  double eta_min = this->universe->eta_at_a(1.0/(1+z_max));
  double eta_max = this->universe->eta_at_a(1.0/(1+z_min));
  double deta = (eta_max - eta_min)/double(n_time - 1);
  double D = this->Newtonian_growth_factor_of_delta[0];
  double D_prime_initial = this->Newtonian_growth_factor_of_delta_prime[0];
  double H_initial = this->universe->H_at_eta(eta_i);
  
  int n = 1+int((delta_max - delta_min)/ddelta);
  
  
  this->delta_values.resize(n, 0.0);
  this->eta_NL.resize(n_time, 0.0);
  this->cylindrical_collapse_evolution_of_delta.resize(n, vector<double>(n_time, 0.0));
  this->cylindrical_collapse_evolution_of_delta_ddelta.resize(n, vector<double>(n_time, 0.0));
  this->cylindrical_collapse_evolution_of_delta_ddelta2.resize(n, vector<double>(n_time, 0.0));
  this->F_prime_of_eta.resize(n_time, 0.0);
  this->F_prime_prime_of_eta.resize(n_time, 0.0);
  
  integration_parameters params;
  params.n_s = this->cosmology.n_s;
  params.Omega_m = this->cosmology.Omega_m;
  params.pointer_to_Universe = this->universe;
  params.pointer_to_Matter = this;
  params.top_hat_radius = 10.0;
  params.second_top_hat_radius = 10.0;
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {F_dF_ddF_cylindrical_wrt_delta_gsl, F_dF_ddF_cylindrical_wrt_delta_jac, 6, (void *) pointer_to_params};
  gsl_odeiv2_system sys_of_derivatives = {dF_cylindrical_ddelta_at_average_density_gsl, dF_cylindrical_ddelta_at_average_density_jac, 4, (void *) pointer_to_params};
    
  for(int i = 0; i < n; i++){
    
    delta += ddelta;
    delta_values[i] = delta;
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    eta_i = this->eta_Newton[0];
    double y[6] = { delta*D, delta*D_prime_initial, D, D_prime_initial, 0.0, 0.0};
    double eta = eta_min;
    for (int j = 0; j < n_time; j++){
      this->eta_NL[j] = eta;
      int status = gsl_odeiv2_driver_apply(d, &eta_i, eta, y);
      if (status != GSL_SUCCESS){
        printf ("error, return value=%d\n", status);
        break;
      }
      this->cylindrical_collapse_evolution_of_delta[i][j] = y[0];
      this->cylindrical_collapse_evolution_of_delta_ddelta[i][j] = y[2];
      this->cylindrical_collapse_evolution_of_delta_ddelta2[i][j] = y[4];
      eta += deta;
    }

    gsl_odeiv2_driver_free(d);
    
  }
  
  //cout << delta_max << setw(20) << this->cylindrical_collapse_evolution_of_delta[n-1][n_time-1] << endl;
  
  double eta_0 = this->universe->eta_at_a(1.0);
  double y[4] = {D, D*H_initial, 0.0, 0.0};
  double eta = eta_min;
  eta_i = this->eta_Newton[0];
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys_of_derivatives, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
  
  for (int j = 0; j < n_time; j++){
    int status = gsl_odeiv2_driver_apply(d, &eta_i, eta, y);
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    this->F_prime_of_eta[j] = y[0];
    this->F_prime_prime_of_eta[j] = y[2];
    eta += deta;
    
  }
  gsl_odeiv2_driver_free(d);
  
  
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

void Matter::set_spherical_collapse_evolution_of_delta(double z_min, double z_max, int n_time){

  double delta_min = -10.0;
  double delta_max = 1.674;
  double ddelta = 0.001;
  double delta = delta_min-ddelta;
  double eta_min = this->universe->eta_at_a(1.0/(1+z_max));
  double eta_max = this->universe->eta_at_a(1.0/(1+z_min));
  double deta = (eta_max - eta_min)/double(n_time - 1);
  
  //double z_i = 500.0;
  double z_i = 1.0/0.00025-1.0;
  double eta_i = this->universe->eta_at_a(1.0/(1+z_i));
  double D = interpolate_neville_aitken(eta_i, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_prime = interpolate_neville_aitken(eta_i, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta_prime, constants::order_of_interpolation);
  double H_initial = this->universe->H_at_eta(eta_i);
  
  int n = 1+int((delta_max - delta_min)/ddelta);
  
  
  this->delta_values_for_spherical_collapse.resize(n, 0.0);
  this->eta_NL_for_spherical_collapse.resize(n_time, 0.0);
  this->spherical_collapse_evolution_of_delta.resize(n, vector<double>(n_time, 0.0));
  this->spherical_collapse_evolution_of_delta_ddelta.resize(n, vector<double>(n_time, 0.0));
  this->spherical_collapse_evolution_of_delta_ddelta2.resize(n, vector<double>(n_time, 0.0));
  this->F_prime_of_eta_for_spherical_collapse.resize(n_time, 0.0);
  this->F_prime_prime_of_eta_for_spherical_collapse.resize(n_time, 0.0);
  
  integration_parameters params;
  params.n_s = this->cosmology.n_s;
  params.Omega_m = this->cosmology.Omega_m;
  params.pointer_to_Universe = this->universe;
  params.pointer_to_Matter = this;
  params.top_hat_radius = 10.0;
  params.second_top_hat_radius = 10.0;
  integration_parameters * pointer_to_params = &params;
  
  
  gsl_odeiv2_system sys = {F_dF_ddF_spherical_wrt_delta_gsl, F_dF_ddF_spherical_wrt_delta_jac, 6, (void *) pointer_to_params};
  gsl_odeiv2_system sys_of_derivatives = {dF_spherical_ddelta_at_average_density_gsl, dF_spherical_ddelta_at_average_density_jac, 4, (void *) pointer_to_params};
  for(int i = 0; i < n; i++){
    
    delta += ddelta;
    delta_values_for_spherical_collapse[i] = delta;
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    eta_i = this->universe->eta_at_a(1.0/(1+z_i));
    double y[6] = { delta*D, delta*D_prime, D, D_prime, 0.0, 0.0};
    double eta = eta_min;
    for (int j = 0; j < n_time; j++){
      this->eta_NL_for_spherical_collapse[j] = eta;
      int status = gsl_odeiv2_driver_apply(d, &eta_i, eta, y);
      if (status != GSL_SUCCESS){
        printf ("error, return value=%d\n", status);
        break;
      }
      this->spherical_collapse_evolution_of_delta[i][j] = y[0];
      //this->spherical_collapse_evolution_of_delta[i][j] = delta*interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
      this->spherical_collapse_evolution_of_delta_ddelta[i][j] = y[2];
      this->spherical_collapse_evolution_of_delta_ddelta2[i][j] = y[4];
      eta += deta;
    }

    gsl_odeiv2_driver_free(d);
    
  }
  
  double eta_0 = this->universe->eta_at_a(1.0);
  double y[4] = {D, D*H_initial, 0.0, 0.0};
  double eta = eta_min;
  eta_i = this->universe->eta_at_a(1.0/(1+z_i));
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys_of_derivatives, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
  
  for (int j = 0; j < n_time; j++){
    int status = gsl_odeiv2_driver_apply(d, &eta_i, eta, y);
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    this->F_prime_of_eta_for_spherical_collapse[j] = y[0];
    this->F_prime_prime_of_eta_for_spherical_collapse[j] = y[2];
    eta += deta;
    
  }
  gsl_odeiv2_driver_free(d);
  

  
}



vector<vector<double> > Matter::compute_phi_of_lambda_3D_with_rescaling(double z, double R, double f_NL, double var_NL_rescale, int modus){
  
  
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  double log_R = log(R);
  double RL, log_RL;
  
  this->current_P_NL = this->P_NL(eta);
  double D = interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_sq = pow(D,2);
  double var_NL_R = variance_of_matter_within_R_NL(R)*var_NL_rescale;
  double var_L_R = variance_of_matter_within_R(R);
  double var_L_RL, skew_L_RL, dvar_L_RL_dR, dskew_L_RL_dR;
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  
  
  vector<double> delta_L_values_EdS;
  vector<double> delta_NL_values_EdS;
  vector<double> delta_NL_prime_values_EdS;
  
  //return_EdS_spherical_collapse(D, &delta_L_values_EdS, &delta_NL_values_EdS, &delta_NL_prime_values_EdS);
  //return_EdS_spherical_collapse(D, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  this->return_delta_NL_of_delta_L_and_dF_ddelta_3D(eta, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  /*
  for(int i = 0; i < delta_L_values_EdS.size(); i++){
    cout << i << '\t';
    cout << delta_L_values_EdS[i] << '\t';
    cout << delta_L_values[i] << '\t';
    cout << delta_NL_values_EdS[i] << '\t';
    cout << delta_NL_values[i] << '\t';
    cout << delta_NL_prime_values_EdS[i] << '\t';
    cout << delta_NL_prime_values[i] << '\n';
  }
  */
  
  vector<vector<double> > data(9, vector<double>(delta_L_values.size(), 0.0));
  
  vector<double> j_values = delta_L_values;
  vector<double> lambda_values = delta_L_values;
  vector<double> phi_values = delta_L_values;
  
  vector<double> j_values_Gauss = delta_L_values;
  vector<double> lambda_values_Gauss = delta_L_values;
  vector<double> phi_values_Gauss = delta_L_values;
  
  int N_radii = this->log_top_hat_radii.size();
  double log_R_min = this->log_top_hat_radii[0];
  double log_R_max = this->log_top_hat_radii[N_radii-1];
  
  double var_L_R_min = this->top_hat_sphere_variances[0];
  double var_L_R_max = this->top_hat_sphere_variances[N_radii-1];
  
  double dvar_L_R_min_dR = this->dtop_hat_sphere_variances_dR[0];
  double dvar_L_R_max_dR = this->dtop_hat_sphere_variances_dR[N_radii-1];
  
  vector<double> *skewness_pointer;
  vector<double> *dskewness_dR_pointer;
  
  switch(modus){
    case 1:
      skewness_pointer = &this->top_hat_sphere_local_skewnesses;
      dskewness_dR_pointer = &this->dtop_hat_sphere_local_skewnesses_dR;
      break;
    case 2:
      skewness_pointer = &this->top_hat_sphere_equilateral_skewnesses;
      dskewness_dR_pointer = &this->dtop_hat_sphere_equilateral_skewnesses_dR;
      break;
    case 3:
      skewness_pointer = &this->top_hat_sphere_orthogonal_skewnesses;
      dskewness_dR_pointer = &this->dtop_hat_sphere_orthogonal_skewnesses_dR;
      break;
  }
  
  for(int d = 0; d < delta_L_values.size(); d++){
    log_RL = log_R + log(1.0+delta_NL_values[d])/3.0;
    RL = exp(log_RL);
    
    var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
    dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
    skew_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, skewness_pointer, constants::order_of_interpolation);
    dskew_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, dskewness_dR_pointer, constants::order_of_interpolation);
    
    
    j_values_Gauss[d] = delta_L_values[d]/var_L_RL;
    lambda_values_Gauss[d] = j_values_Gauss[d]/delta_NL_prime_values[d];
    lambda_values_Gauss[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values_Gauss[d],2);
    phi_values_Gauss[d] = lambda_values_Gauss[d]*delta_NL_values[d]-delta_L_values[d]*j_values_Gauss[d] + var_L_RL/2.0*pow(j_values_Gauss[d],2);
    
    skew_L_RL *= f_NL;
    dskew_L_RL_dR *= f_NL;
    
    j_values_Gauss[d] = delta_L_values[d]/var_L_RL;
    lambda_values_Gauss[d] = j_values_Gauss[d]/delta_NL_prime_values[d];
    lambda_values_Gauss[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values_Gauss[d],2);
    phi_values_Gauss[d] = lambda_values_Gauss[d]*delta_NL_values[d]-delta_L_values[d]*j_values_Gauss[d] + var_L_RL/2.0*pow(j_values_Gauss[d],2);
    
    if(f_NL != 0.0){
      
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
      j_values[d] = j_values_Gauss[d];
      lambda_values[d] = lambda_values_Gauss[d];
      phi_values[d] = phi_values_Gauss[d];
    }
    
    phi_values[d] *= D_sq*var_L_R/var_NL_R;
    phi_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
    
    //phi_values[d] += 0.5*lambda_values[d]*lambda_values[d]*(var_NL_R-D_sq*var_L_R);
    //phi_values_Gauss[d] += 0.5*lambda_values_Gauss[d]*lambda_values_Gauss[d]*(var_NL_R-D_sq*var_L_R);
    
    data[0][d] = delta_L_values[d];
    data[1][d] = delta_NL_values[d];
    data[2][d] = lambda_values[d];
    data[3][d] = phi_values[d];
    data[4][d] = lambda_values_Gauss[d];
    data[5][d] = phi_values_Gauss[d];
    data[6][d] = var_L_RL;
    data[7][d] = skew_L_RL;
    data[8][d] = RL*constants::c_over_e5;
    
  }
  
  return data;
  
}


vector<vector<double> > Matter::compute_phi_of_lambda_3D_with_rescaling_and_S3_powerlaw(double z, double R, double R_0, double A_S3, double n_S3, double var_NL_rescale){
  
  
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  double log_R = log(R);
  double RL, log_RL;
  
  this->current_P_NL = this->P_NL(eta);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  double D = interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_sq = pow(D,2);
  double var_NL_R = variance_of_matter_within_R_NL(R)*var_NL_rescale;
  double var_L_R = variance_of_matter_within_R(R);
  double var_L_RL, skew_L_RL, dvar_L_RL_dR, dskew_L_RL_dR;
  
  cout << R*constants::c_over_e5 << setw(20);
  cout << D << setw(20);
  cout << var_L_R << setw(20);
  cout << var_L_R*D_sq << setw(20);
  cout << var_NL_R << '\n';
  cout.flush();
  
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  this->return_delta_NL_of_delta_L_and_dF_ddelta_3D(eta, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  
  vector<vector<double> > data(9, vector<double>(delta_L_values.size(), 0.0));
  
  vector<double> j_values = delta_L_values;
  vector<double> lambda_values = delta_L_values;
  vector<double> phi_values = delta_L_values;
  
  vector<double> j_values_Gauss = delta_L_values;
  vector<double> lambda_values_Gauss = delta_L_values;
  vector<double> phi_values_Gauss = delta_L_values;
  
  cout << setprecision(10);
  
  int N_radii = this->log_top_hat_radii.size();
  double log_R_min = this->log_top_hat_radii[0];
  double log_R_max = this->log_top_hat_radii[N_radii-1];
  
  double var_L_R_min = this->top_hat_sphere_variances[0];
  double var_L_R_max = this->top_hat_sphere_variances[N_radii-1];
  
  double dvar_L_R_min_dR = this->dtop_hat_sphere_variances_dR[0];
  double dvar_L_R_max_dR = this->dtop_hat_sphere_variances_dR[N_radii-1];
  
  for(int d = 0; d < delta_L_values.size(); d++){
    log_RL = log_R + log(1.0+delta_NL_values[d])/3.0;
    RL = exp(log_RL);
    
    if(log_RL >= log_R_min && log_RL <= log_R_max){
      var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
      dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
    }
    else if(log_RL < log_R_min){
      var_L_RL = var_L_R_min*exp(2.0*(log_R_min-log_RL));
      dvar_L_RL_dR = -2.0*var_L_RL/RL;
    }
    else{
      var_L_RL = var_L_R_max*exp(2.0*(log_R_max-log_RL));
      dvar_L_RL_dR = -2.0*var_L_RL/RL;
    }
    
    double alpha = 1.5;
    skew_L_RL = A_S3*pow(RL/R_0, n_S3)*pow(var_L_RL, alpha);
    dskew_L_RL_dR = skew_L_RL*(n_S3/RL + alpha*dvar_L_RL_dR/var_L_RL);
    
    j_values_Gauss[d] = delta_L_values[d]/var_L_RL;
    lambda_values_Gauss[d] = j_values_Gauss[d]/delta_NL_prime_values[d];
    lambda_values_Gauss[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values_Gauss[d],2);
    phi_values_Gauss[d] = lambda_values_Gauss[d]*delta_NL_values[d]-delta_L_values[d]*j_values_Gauss[d] + var_L_RL/2.0*pow(j_values_Gauss[d],2);
    
    if(A_S3 != 0.0){
      
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
      j_values[d] = j_values_Gauss[d];
      lambda_values[d] = lambda_values_Gauss[d];
      phi_values[d] = phi_values_Gauss[d];
    }
    
    phi_values[d] *= D_sq*var_L_R/var_NL_R;
    phi_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
    
    //phi_values[d] += 0.5*lambda_values[d]*lambda_values[d]*(var_NL_R-D_sq*var_L_R);
    //phi_values_Gauss[d] += 0.5*lambda_values_Gauss[d]*lambda_values_Gauss[d]*(var_NL_R-D_sq*var_L_R);
    
    data[0][d] = delta_L_values[d];
    data[1][d] = delta_NL_values[d];
    data[2][d] = lambda_values[d];
    data[3][d] = phi_values[d];
    data[4][d] = lambda_values_Gauss[d];
    data[5][d] = phi_values_Gauss[d];
    data[6][d] = var_L_RL;
    data[7][d] = skew_L_RL;
    data[8][d] = RL*constants::c_over_e5;
    
  }
  
  return data;
  
}











