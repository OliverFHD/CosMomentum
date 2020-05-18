
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_legendre.h>

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
 *                                         ****************
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
  this->print_growth_history("linear_growth_history.dat");
  this->initialize_up_to_second_order_growth_factor_of_delta(2.0/7.0*this->D_initial*this->D_initial, 4.0/7.0*this->D_initial*this->D_prime_initial);
  
  double z_of_f_NL_rescale = 500.0;
  double a_of_f_NL_rescale = 1.0/(1.0+z_of_f_NL_rescale);
  double eta_of_f_NL_rescale = this->universe->eta_at_a(a_of_f_NL_rescale);
  this->f_NL_rescaling_factor = this->return_D_of_eta(eta_of_f_NL_rescale)/a_of_f_NL_rescale;
  
  this->set_wave_numbers();
  if(this->cosmology.Omega_b != 0.0){
    this->set_transfer_function_Eisenstein_and_Hu();
  }
  else{
    this->set_transfer_function_Bond_and_Efstathiou();
  }
  this->norm = this->variance_of_matter_within_R_before_norm_was_determined(8.0);
  this->set_P_today();
  
  cout << "Setting sphere variances.\n";
  this->set_sphere_variances();
  cout << "Setting cylinder variances.\n";
  this->set_cylinder_variances();
  this->skewness_initialisation = UNINITIALISED;
  this->skewness_initialisation_2D = UNINITIALISED;
  
  cout << "Setting spherical collapse evolution.\n";
  this->set_spherical_collapse_evolution_of_delta(0.0, 3.0, 1000);
  cout << "Setting cylindrical collapse evolution.\n";
  this->set_cylindrical_collapse_evolution_of_delta(0.0, 3.0, 1000);
  cout << "Done.\n";
  //this->set_phi_tilde_for_LOS_integration(20.0*constants::pi/180.0/60.0, 0.0);
  //cout << "Done.\n";
  
  /*cout << "loading pofz data\n";
  vector<vector<double> > nofz_data = read_table_double("Data/redshift_distributions/nofz_redmagic_0.2_0.45_true.dat");
  
  int Nz = nofz_data.size();
  vector<double> z_values(Nz, 0.0);
  vector<double> nofz_values(Nz, 0.0);
  for(int i = 0; i < Nz; i++){
    z_values[i] = nofz_data[i][0];
    nofz_values[i] = nofz_data[i][1];
    cout << i << "   " << z_values[i] << "   " << nofz_values[i] << "\n";
  }
  
  cout << "projecting CGF\n";
  
  this->return_LOS_integrated_phi_of_lambda(20.0*constants::arcmin, 0.0, z_values, nofz_values);
  */
  
  
}

Matter::Matter(Universe* uni, string file_for_transfer_function){
  
  this->universe = uni;
  this->cosmology = this->universe->return_cosmology(); 
  
  this->a_initial = this->universe->return_a_initial();
  this->a_final = this->universe->return_a_final();
  this->eta_initial = this->universe->return_eta_initial();
  this->eta_final = this->universe->return_eta_final();
  
  this->set_initial_conditions();
  this->initialize_linear_growth_factor_of_delta();
  this->print_growth_history("linear_growth_history.dat");
  this->initialize_up_to_second_order_growth_factor_of_delta(2.0/7.0*this->D_initial*this->D_initial, 4.0/7.0*this->D_initial*this->D_prime_initial);
  
  double z_of_f_NL_rescale = 500.0;
  double a_of_f_NL_rescale = 1.0/(1.0+z_of_f_NL_rescale);
  double eta_of_f_NL_rescale = this->universe->eta_at_a(a_of_f_NL_rescale);
  this->f_NL_rescaling_factor = this->return_D_of_eta(eta_of_f_NL_rescale)/a_of_f_NL_rescale;
  
  this->set_wave_numbers();
  this->set_transfer_function_from_file(file_for_transfer_function);
  this->norm = this->variance_of_matter_within_R_before_norm_was_determined(8.0);
  this->set_P_today();
  
  cout << "Setting sphere variances.\n";
  this->set_sphere_variances();
  cout << "Setting cylinder variances.\n";
  this->set_cylinder_variances();
  this->skewness_initialisation = UNINITIALISED;
  this->skewness_initialisation_2D = UNINITIALISED;
  
  cout << "Setting spherical collapse evolution.\n";
  this->set_spherical_collapse_evolution_of_delta(0.0, 3.0, 1000);
  cout << "Setting cylindrical collapse evolution.\n";
  this->set_cylindrical_collapse_evolution_of_delta(0.0, 3.0, 1000);
  cout << "Done.\n";
  //this->set_phi_tilde_for_LOS_integration(20.0*constants::pi/180.0/60.0, 0.0);
  //cout << "Done.\n";
  
  
  /*vector<double> delta_L_values0;
  vector<double> delta_NL_values0;
  vector<double> delta_NL_prime_values0;
  vector<double> delta_L_values1;
  vector<double> delta_NL_values1;
  vector<double> delta_NL_prime_values1;
  double a0 = 1.0/(1.0+0.0);
  double a1 = 1.0/(1.0+1.0);
  double eta0 = this->universe->eta_at_a(a0);
  double eta1 = this->universe->eta_at_a(a1);
  
  double D0 = this->return_D_of_eta(eta0);
  double D1 = this->return_D_of_eta(eta1);
  
  cout << D0 << "  " << D1 << '\n';
  
  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(eta0, &delta_L_values0, &delta_NL_values0, &delta_NL_prime_values0);
  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(eta1, &delta_L_values1, &delta_NL_values1, &delta_NL_prime_values1);
  
  cout << scientific << setprecision(10);
  cout << "# No      delta_L(z=0)      delta_NL(z=0)      delta_L(z=1)      delta_NL(z=1)\n";
  for(int i = 0; i < delta_L_values0.size(); i++){
    cout << setw(20) << i;
    cout << setw(20) << delta_L_values0[i];
    cout << setw(20) << delta_NL_values0[i];
    cout << setw(20) << delta_L_values1[i]*D1;
    cout << setw(20) << delta_NL_values1[i];
    cout << '\n';
  }
  */
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

void Matter::initialize_linear_growth_factor_of_delta(){

  integration_parameters params;
  params.pointer_to_Matter = this;
  params.pointer_to_Universe = this->universe;
  params.Omega_m = this->cosmology.Omega_m;
  integration_parameters * pointer_to_params = &params;
  
  gsl_odeiv2_system sys = {growth_factor_gsl, growth_factor_gsl_jac, 2, (void *) pointer_to_params};

  
  /* THIS HARDCODING IS REALLY BAD!! */
  /* CHANGE REQUIRED                 */
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
  double delta_max = 1.67;
  double ddelta = 0.01;
  double delta = delta_min-ddelta;
  double eta_min = this->universe->eta_at_a(1.0/(1+z_max));
  double eta_max = this->universe->eta_at_a(1.0/(1+z_min));
  double deta = (eta_max - eta_min)/double(n_time - 1);
  
  /* THIS HARDCODING IS REALLY BAD!! */
  /* CHANGE REQUIRED                 */
  double a_i = 0.00025;
  double eta_i = std::max(this->universe->eta_at_a(0.00025), this->eta_Newton[0]);
  double D = interpolate_neville_aitken(eta_i, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_prime = interpolate_neville_aitken(eta_i, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta_prime, constants::order_of_interpolation);
  double H_initial = this->universe->H_at_eta(eta_i);
  
  int n = 1+int((delta_max - delta_min)/ddelta);
  
  
  this->delta_values_for_spherical_collapse.resize(n, 0.0);
  this->eta_NL_for_spherical_collapse.resize(n_time, 0.0);
  this->spherical_collapse_evolution_of_delta.resize(n, vector<double>(n_time, 0.0));
  this->spherical_collapse_evolution_of_delta_ddelta.resize(n, vector<double>(n_time, 0.0));
  this->spherical_collapse_evolution_of_delta_ddelta2.resize(n, vector<double>(n_time, 0.0));
  
  integration_parameters params;
  params.n_s = this->cosmology.n_s;
  params.Omega_m = this->cosmology.Omega_m;
  params.pointer_to_Universe = this->universe;
  params.pointer_to_Matter = this;
  params.top_hat_radius = 10.0;
  params.second_top_hat_radius = 10.0;
  integration_parameters * pointer_to_params = &params;
  
  
  gsl_odeiv2_system sys = {F_dF_ddF_spherical_wrt_delta_gsl, F_dF_ddF_spherical_wrt_delta_jac, 6, (void *) pointer_to_params};
  for(int i = 0; i < n; i++){
    
    delta += ddelta;
    delta_values_for_spherical_collapse[i] = delta;
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    eta_i = std::max(this->universe->eta_at_a(0.00025), this->eta_Newton[0]);
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
      this->spherical_collapse_evolution_of_delta_ddelta[i][j] = y[2];
      this->spherical_collapse_evolution_of_delta_ddelta2[i][j] = y[4];
      eta += deta;
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

void Matter::set_cylindrical_collapse_evolution_of_delta(double z_min, double z_max, int n_time){

  double delta_min = -10.0;
  double delta_max = 1.45;
  double ddelta = 0.01;
  double delta = delta_min-ddelta;
  double eta_min = this->universe->eta_at_a(1.0/(1+z_max));
  double eta_max = this->universe->eta_at_a(1.0/(1+z_min));
  double deta = (eta_max - eta_min)/double(n_time - 1);
  
  /* THIS HARDCODING IS REALLY BAD!! */
  /* CHANGE REQUIRED                 */
  double a_i = 0.00025;
  double eta_i = std::max(this->universe->eta_at_a(0.00025), this->eta_Newton[0]);
  double D = interpolate_neville_aitken(eta_i, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_prime = interpolate_neville_aitken(eta_i, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta_prime, constants::order_of_interpolation);
  double H_initial = this->universe->H_at_eta(eta_i);
  
  int n = 1+int((delta_max - delta_min)/ddelta);
  
  
  this->delta_values_for_cylindrical_collapse.resize(n, 0.0);
  this->eta_NL_for_cylindrical_collapse.resize(n_time, 0.0);
  this->cylindrical_collapse_evolution_of_delta.resize(n, vector<double>(n_time, 0.0));
  this->cylindrical_collapse_evolution_of_delta_ddelta.resize(n, vector<double>(n_time, 0.0));
  this->cylindrical_collapse_evolution_of_delta_ddelta2.resize(n, vector<double>(n_time, 0.0));
  
  integration_parameters params;
  params.n_s = this->cosmology.n_s;
  params.Omega_m = this->cosmology.Omega_m;
  params.pointer_to_Universe = this->universe;
  params.pointer_to_Matter = this;
  params.top_hat_radius = 10.0;
  params.second_top_hat_radius = 10.0;
  integration_parameters * pointer_to_params = &params;
  
  
  gsl_odeiv2_system sys = {F_dF_ddF_cylindrical_wrt_delta_gsl, F_dF_ddF_cylindrical_wrt_delta_jac, 6, (void *) pointer_to_params};
  for(int d = 0; d < n; d++){
    
    delta += ddelta;
    delta_values_for_cylindrical_collapse[d] = delta;
    
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    eta_i = std::max(this->universe->eta_at_a(0.00025), this->eta_Newton[0]);
    double y[6] = { delta*D, delta*D_prime, D, D_prime, 0.0, 0.0};
    double eta = eta_min;
    for (int t = 0; t < n_time; t++){
      this->eta_NL_for_cylindrical_collapse[t] = eta;
      int status = gsl_odeiv2_driver_apply(driver, &eta_i, eta, y);
      if (status != GSL_SUCCESS){
        printf ("error, return value=%d\n", status);
        break;
      }
      this->cylindrical_collapse_evolution_of_delta[d][t] = y[0];
      this->cylindrical_collapse_evolution_of_delta_ddelta[d][t] = y[2];
      this->cylindrical_collapse_evolution_of_delta_ddelta2[d][t] = y[4];
      eta += deta;
    }

    gsl_odeiv2_driver_free(driver);
    
  }
  

  
}

vector<vector<double> > Matter::compute_PDF_3D(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
  cout << "Computing phi_data:\n";
  cout.flush();
    
  vector<vector<double> > phi_data = this->compute_phi_of_lambda_3D(z, R_in_Mpc_over_h/constants::c_over_e5, f_NL, var_NL_rescale);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[2][0];
  double delta_c = phi_data[1][0];
  double tau_c = 0.0;
  for(int i = 1; i < phi_data[4].size(); i++){
    if(phi_data[4][i-1] < phi_data[4][i]){
      n_lambda = i+1;
      lambda_c = phi_data[4][i];
      delta_c = phi_data[1][i];
      if(phi_data[1][i]*phi_data[4][i] > phi_data[3][i])
        tau_c = sqrt(2.0*(phi_data[1][i]*phi_data[4][i] - phi_data[3][i]));
    }
    else{
      i = 2*phi_data[4].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> delta_L_values(n_lambda, 0.0);
  vector<double> delta_NL_values(n_lambda, 0.0);
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> R_L_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    delta_L_values[i] = phi_data[0][i];
    delta_NL_values[i] = phi_data[1][i];
    lambda_values[i] = phi_data[2][i];
    phi_values[i] = phi_data[3][i];
    R_L_values[i] = phi_data[8][i];
    if(phi_data[1][i]*phi_data[2][i] > phi_data[3][i]){
      tau_values[i] = sqrt(2.0*(phi_data[1][i]*phi_data[2][i] - phi_data[3][i]));
      if(lambda_values[i] < 0.0) tau_values[i] *= -1.0;
    }
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  
  double R_max = interpolate_neville_aitken(tau_max, &tau_values, &R_L_values, constants::order_of_interpolation);
  
  double R_min = interpolate_neville_aitken(tau_min, &tau_values, &R_L_values, constants::order_of_interpolation);
  
  
  if(R_max > exp(this->log_top_hat_radii_for_skewnesses[this->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5){
    tau_max = min(tau_max, interpolate_neville_aitken(exp(this->log_top_hat_radii_for_skewnesses[this->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  if(R_min < exp(this->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5){
    tau_min = max(tau_min, interpolate_neville_aitken(exp(this->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  tau_c = tau_max;
  
  delta_c = interpolate_neville_aitken(tau_max, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  lambda_c = interpolate_neville_aitken(tau_max, &tau_values, &lambda_values, constants::order_of_interpolation);
  
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
  vector<double> phi_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  
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
    lambda_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &lambda_values, constants::order_of_interpolation);
    phi_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &phi_values, constants::order_of_interpolation);
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &delta_NL_values, constants::order_of_interpolation);
  }
  
  
  cout << "Done.\n";
  
  /*
   * Express functions as polynomials in tau.
   * 
   */
  cout << "Computing tau coefficients:\n";
  
  int n_coeff = constants::generating_function_coeff_order;
  
  vector<double> coefficients_lambda_of_tau = return_coefficients(&tau_for_fit, &lambda_for_fit, n_coeff);
  vector<double> coefficients_lambda_of_tau_prime(coefficients_lambda_of_tau.size(), 0.0);
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_lambda_of_tau_prime[i] = coefficients_lambda_of_tau[i+1]*double(i+1);
  
  vector<double> coefficients_phi_of_tau = return_coefficients(&tau_for_fit, &phi_for_fit, n_coeff);
  vector<double> coefficients_phi_of_tau_prime = return_coefficients(&tau_for_fit, &phi_prime_for_fit, n_coeff);
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = this->return_non_linear_variance(z, R_in_Mpc_over_h);
  double delta_min = interpolate_neville_aitken(tau_min, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  double delta_max = max(delta_c, 6.0*sqrt(var));
  // --> ISSUE: choosing delta_max to be 6*sigma may not be 100% reasonable for very skewed PDFs
  //            Maybe choose it by some quantile in a log-normal PDF that approximates the real PDF?
  
  double ddelta = (delta_max-delta_min)/double(n_delta-1);
  double delta;
  double tau_0, lambda_0;
  double dr;
  
  complex<double> lambda;
  complex<double> lambda_next;
  complex<double> tau, tau_next;
  complex<double> phi_prime, phi_prime_next;
  complex<double> exponent, exponent_next;
  complex<double> dlambda;
  complex<double> step, first_step;
  
  vector<vector<double> > PDF_data(2, vector<double>(n_delta));
  
  cout << "Computing PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    
    PDF_data[0][i] = delta;
    PDF_data[1][i] = 0.0;
      
    if(delta < delta_c){
      tau_0 = interpolate_Newton(delta, &delta_NL_values, &tau_values, constants::order_of_interpolation);
      lambda_0 = interpolate_Newton(delta, &delta_NL_values, &lambda_values, constants::order_of_interpolation);
    }
    else{
      tau_0 = tau_c;
      lambda_0 = lambda_c;
    }
    lambda = complex<double>(lambda_0, 0.0);
    tau = complex<double>(tau_0, 0.0);
    exponent = exp(-lambda*delta + return_polnomial_value(tau, &coefficients_phi_of_tau));
    
    // sigma_r^2 \approx 1/phi''(lambda_0)
    double sigma_frac = 0.001;
    dr = sigma_frac/sqrt(interpolate_neville_aitken_derivative(lambda_0, &lambda_values, &delta_NL_values, constants::order_of_interpolation));
    dlambda = complex<double>(0.0, dr);
    int j = 0;
    do{
      lambda_next = lambda + 0.5*dlambda;
      tau_next = Newtons_method_complex(tau, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda_next = lambda + dlambda;
      tau_next = Newtons_method_complex(tau_next, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      exponent_next = exp(-lambda_next*delta + return_polnomial_value(tau_next, &coefficients_phi_of_tau));
      
      step = 0.5*dlambda*(exponent_next+exponent);
      PDF_data[1][i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      exponent = exponent_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF_data[1][i] /= constants::pi;
    
  }
  cout << "Done.\n";
  
  return PDF_data;
  
}


vector<vector<double> > Matter::compute_phi_of_lambda_3D(double z, double R, double f_NL, double var_NL_rescale){
  
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
    
    var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
    dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
    skew_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->top_hat_sphere_skewnesses, constants::order_of_interpolation);
    dskew_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->dtop_hat_sphere_skewnesses_dR, constants::order_of_interpolation);
    
    j_values_Gauss[d] = delta_L_values[d]/var_L_RL;
    lambda_values_Gauss[d] = j_values_Gauss[d]/delta_NL_prime_values[d];
    lambda_values_Gauss[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values_Gauss[d],2);
    phi_values_Gauss[d] = lambda_values_Gauss[d]*delta_NL_values[d]-delta_L_values[d]*j_values_Gauss[d] + var_L_RL/2.0*pow(j_values_Gauss[d],2);
    
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
      j_values[d] = j_values_Gauss[d];
      lambda_values[d] = lambda_values_Gauss[d];
      phi_values[d] = phi_values_Gauss[d];
    }
    
    phi_values[d] *= D_sq*var_L_R/var_NL_R;
    phi_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
        
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



vector<vector<double> > Matter::compute_phi_of_lambda_3D_EdS(double z, double R, double f_NL, double var_NL_rescale){
  
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
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  
  return_EdS_spherical_collapse(D, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  
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
  
  for(int d = 0; d < delta_L_values.size(); d++){
    log_RL = log_R + log(1.0+delta_NL_values[d])/3.0;
    RL = exp(log_RL);
    
    var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
    dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
    skew_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->top_hat_sphere_skewnesses, constants::order_of_interpolation);
    dskew_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_radii_for_skewnesses, &this->dtop_hat_sphere_skewnesses_dR, constants::order_of_interpolation);
    
    j_values_Gauss[d] = delta_L_values[d]/var_L_RL;
    lambda_values_Gauss[d] = j_values_Gauss[d]/delta_NL_prime_values[d];
    lambda_values_Gauss[d] -= RL/(3.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values_Gauss[d],2);
    phi_values_Gauss[d] = lambda_values_Gauss[d]*delta_NL_values[d]-delta_L_values[d]*j_values_Gauss[d] + var_L_RL/2.0*pow(j_values_Gauss[d],2);
    
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
      j_values[d] = j_values_Gauss[d];
      lambda_values[d] = lambda_values_Gauss[d];
      phi_values[d] = phi_values_Gauss[d];
    }
    
    phi_values[d] *= D_sq*var_L_R/var_NL_R;
    phi_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values[d] *= D_sq*var_L_R/var_NL_R;
    lambda_values_Gauss[d] *= D_sq*var_L_R/var_NL_R;
        
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



vector<vector<double> > Matter::compute_phi_of_lambda_2D(double z, double R, double L, double f_NL, double var_NL_rescale, int LOS_modus){
  
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  double log_R = log(R);
  double RL, log_RL;
  
  this->current_P_NL = this->P_NL(eta);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  double D = interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_sq = pow(D,2);
  cout << "Test 1\n";
  double var_NL_R;
  // LOS_modus: 1 == delta function, 2 == Gauss, 3 == tophat
  if(LOS_modus == 1)
    var_NL_R = variance_of_matter_within_R_NL_2D(R)/L*var_NL_rescale;
  else if(LOS_modus == 2)
    var_NL_R = variance_of_matter_within_R_NL_2D_GaussianLOS(R, L)*var_NL_rescale;
  else if(LOS_modus == 3)
    var_NL_R = variance_of_matter_within_R_NL_2D(R, L)*var_NL_rescale;
  else{
    cerr << "ERROR: LOS_modus not valid in Matter::compute_phi_of_lambda_2D.\n";
    cerr << "LOS_modus = " << LOS_modus << "\n";
    exit(1);
  }
  double var_L_R = variance_of_matter_within_R_2D(R)/L;
  double var_L_RL, skew_L_RL, dvar_L_RL_dR, dskew_L_RL_dR;
  cout << "Test 2\n";
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  
  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(eta, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  
  cout << "Test 3\n";
  vector<vector<double> > data(9, vector<double>(delta_L_values.size(), 0.0));
  
  vector<double> j_values = delta_L_values;
  vector<double> lambda_values = delta_L_values;
  vector<double> phi_values = delta_L_values;
  
  vector<double> j_values_Gauss = delta_L_values;
  vector<double> lambda_values_Gauss = delta_L_values;
  vector<double> phi_values_Gauss = delta_L_values;
  
  int N_radii = this->log_top_hat_radii.size();
  double log_R_min = this->log_top_hat_cylinder_radii[0];
  double log_R_max = this->log_top_hat_cylinder_radii[N_radii-1];
  
  double var_L_R_min = this->top_hat_cylinder_variances[0];
  double var_L_R_max = this->top_hat_cylinder_variances[N_radii-1];
  
  double dvar_L_R_min_dR = this->dtop_hat_cylinder_variances_dR[0];
  double dvar_L_R_max_dR = this->dtop_hat_cylinder_variances_dR[N_radii-1];
  
  cout << "Test 4\n";
  for(int d = 0; d < delta_L_values.size(); d++){
    log_RL = log_R + log(1.0+delta_NL_values[d])/2.0;
    RL = exp(log_RL);
    
    var_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii, &this->top_hat_cylinder_variances, constants::order_of_interpolation);
    dvar_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii, &this->dtop_hat_cylinder_variances_dR, constants::order_of_interpolation);
    
    var_L_RL /= L;
    dvar_L_RL_dR /= L;
    
    j_values_Gauss[d] = delta_L_values[d]/var_L_RL;
    lambda_values_Gauss[d] = j_values_Gauss[d]/delta_NL_prime_values[d];
    lambda_values_Gauss[d] -= RL/(2.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values_Gauss[d],2);
    phi_values_Gauss[d] = lambda_values_Gauss[d]*delta_NL_values[d]-delta_L_values[d]*j_values_Gauss[d] + var_L_RL/2.0*pow(j_values_Gauss[d],2);
    
    
    if(f_NL != 0.0){
      
      skew_L_RL = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii_for_skewnesses, &this->top_hat_cylinder_skewnesses, constants::order_of_interpolation);
      dskew_L_RL_dR = interpolate_neville_aitken(log_RL, &this->log_top_hat_cylinder_radii_for_skewnesses, &this->dtop_hat_cylinder_skewnesses_dR, constants::order_of_interpolation);
    
      skew_L_RL *= f_NL/(L*L);
      dskew_L_RL_dR *= f_NL/(L*L);
      
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
      
      j_values[d] = -var_L_RL/skew_L_RL;
      if(skew_L_RL > 0.0)
        j_values[d] += sqrt(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL);
      else
        j_values[d] -= sqrt(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL);
      
      lambda_values[d] = j_values[d]/delta_NL_prime_values[d];
      lambda_values[d] -= RL/(2.0*(1.0+delta_NL_values[d]))*(dvar_L_RL_dR/2.0*pow(j_values[d],2) + dskew_L_RL_dR/6.0*pow(j_values[d],3));
      
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
  cout << "Test 5\n";
  
  return data;
  
}



void Matter::compute_phi_tilde_of_lambda_2D(double eta, double R, double f_NL, vector<double> * lambda_of_delta, vector<double> * phi_of_delta, vector<double> * phi_prime_of_delta){
  
  double log_R = log(R);
  double RL, log_RL;
  
  this->current_P_NL = this->P_NL(eta);
  this->current_P_L = this->P_L(this->universe->eta_at_a(1.0));
  double D = interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_sq = pow(D,2);
  
  double var_NL_R = variance_of_matter_within_R_NL_2D(R);
  double var_L_R = variance_of_matter_within_R_2D(R);
  double var_L_RL, skew_L_RL, dvar_L_RL_dR, dskew_L_RL_dR;
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  
  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(eta, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  
  if(lambda_of_delta->size()!=delta_L_values.size()) (*lambda_of_delta) = delta_L_values;
  if(phi_of_delta->size()!=delta_L_values.size()) (*phi_of_delta) = delta_L_values;
  (*phi_prime_of_delta) = delta_NL_values;
  
  vector<double> j_values = delta_L_values;
  vector<double> j_values_Gauss = delta_L_values;
  
  int N_radii = this->log_top_hat_radii.size();
  double log_R_min = this->log_top_hat_cylinder_radii[0];
  double log_R_max = this->log_top_hat_cylinder_radii[N_radii-1];
  
  double var_L_R_min = this->top_hat_cylinder_variances[0];
  double var_L_R_max = this->top_hat_cylinder_variances[N_radii-1];
  
  double dvar_L_R_min_dR = this->dtop_hat_cylinder_variances_dR[0];
  double dvar_L_R_max_dR = this->dtop_hat_cylinder_variances_dR[N_radii-1];
  
  for(int d = 0; d < delta_L_values.size(); d++){
    log_RL = log_R + log(1.0+delta_NL_values[d])/2.0;
    RL = exp(log_RL);
    
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
      
      j_values[d] = -var_L_RL/skew_L_RL;
      if(skew_L_RL > 0.0)
        j_values[d] += sqrt(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL);
      else
        j_values[d] -= sqrt(pow(var_L_RL/skew_L_RL,2) + 2.0*delta_L_values[d]/skew_L_RL);
      
      (*lambda_of_delta)[d] = j_values[d]/delta_NL_prime_values[d];
      (*lambda_of_delta)[d] -= RL/(2.0*(1.0+delta_NL_values[d]))*(dvar_L_RL_dR/2.0*pow(j_values[d],2) + dskew_L_RL_dR/6.0*pow(j_values[d],3));
      
      (*phi_of_delta)[d] = -(*lambda_of_delta)[d]*delta_NL_values[d]+delta_L_values[d]*j_values[d] - (var_L_RL/2.0*pow(j_values[d],2) + skew_L_RL/6.0*pow(j_values[d],3));
      (*phi_of_delta)[d] *= -1.0;
    }
    else{
      j_values[d] = delta_L_values[d]/var_L_RL;
      (*lambda_of_delta)[d] = j_values[d]/delta_NL_prime_values[d];
      (*lambda_of_delta)[d] -= RL/(2.0*(1.0+delta_NL_values[d]))*dvar_L_RL_dR/2.0*pow(j_values[d],2);
      (*phi_of_delta)[d] = (*lambda_of_delta)[d]*delta_NL_values[d]-delta_L_values[d]*j_values[d] + var_L_RL/2.0*pow(j_values[d],2);
    }
    
    (*phi_of_delta)[d] *= D_sq*var_L_R/var_NL_R;
    (*lambda_of_delta)[d] *= D_sq*var_L_R/var_NL_R;
    
  }  
  
}








