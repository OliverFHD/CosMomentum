


/*******************************************************************************************************************************************************
* These functions are for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************
* These functions are for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************
* These functions are for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************
* These functions are for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/


/*******************************************************
 *******************************************************
 **__________ 4. FITTING FORMULA OF SMITH + __________**
 *******************************************************
 ******************************************************* 
 *                                                     *
 * ..... 4.3  sig_sq                                   *
 * ..... 4.4  c_and_n_NL                               *
 * ..... 4.5  k_NL                                     *
 * ..... 4.6  Delta_Q_sq                               *
 * ..... 4.7  Delta_H_sq                               *
 * ..... 4.8  Delta_H_prime_sq                         *
 * ..... 4.9  P_NL_at                                  *
 * ..... 4.10 P_NL                                     *
 *                                                     *
 *******************************************************
 *******************************************************/

/*******************************************************************************************************************************************************
 * 4.3  sig_sq
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::sig_sq(double R, double e){

	double s_sq = 0;
  double prefactors;
  double D_sq = interpolate_neville_aitken(e, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
	integration_parameters params;
  integration_parameters * pointer_to_params = &params;
	
  D_sq *= D_sq;
  prefactors = D_sq/(2.0*constants::pi_sq)*pow(this->cosmology.sigma_8, 2)/this->norm;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
  
  s_sq = int_gsl_integrate_medium_precision(halofit_sig_sq_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  return prefactors*s_sq;
  
}

/*******************************************************************************************************************************************************
 * 4.4  c_and_n_NL
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

vector<double> Matter::c_and_n_NL(double R, double e){
  
  double lnk_max = log(maximal_wave_number_in_H0_units);
  double prefactors;
  double D_sq = interpolate_neville_aitken(e, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);  

	integration_parameters params;
  integration_parameters * pointer_to_params = &params;

  vector<double> integral(2, 0.0);
  
	
  D_sq *= D_sq;
  prefactors = D_sq/(2*pi_sq)*pow(this->cosmology.sigma_8, 2)/this->norm;
  params.top_hat_radius = R;
  params.n_s = this->cosmology.n_s;
  params.pointer_to_Matter = this;
	
	integral[0] = 2.0*prefactors*int_gsl_integrate_medium_precision(halofit_C_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  integral[1] = 4.0*prefactors*int_gsl_integrate_medium_precision(halofit_n_gsl,(void*)pointer_to_params,log(minimal_wave_number_in_H0_units),log(maximal_wave_number_in_H0_units),NULL,1000);
  
  return integral;
  
}


/*******************************************************************************************************************************************************
 * 4.5 k_NL
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::k_NL(double k_min, double k_max, double e){
  
  double k = (k_max + k_min)/2.0;
  double s = this->sig_sq(1.0/k, e);  
  

  if(abs(s-1.0) <= sig_sq_precision){
    return k;
  }
  else if(s-1.0 > sig_sq_precision){
    return this->k_NL(k_min, k, e);
  }
  else{
    s = this->sig_sq(1.0/k_max, e);
    if(s-1.0 > sig_sq_precision)
      return this->k_NL(k, k_max, e);
  }

  
  return this->k_NL(0.9*k_max, 2*k_max, e);
  
}


/*******************************************************************************************************************************************************
 * 4.6 Delta_Q_sq
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::Delta_Q_sq(double k, double e){
  
  double h = this->cosmology.h_100;

  double n = this->cosmology.n_s;
  double k_s = this->current_k_non_linear;
  double y = k/k_s;
  
  double Delta_L_sq = this->Newtonian_linear_power_spectrum(k, e)*pow(k, 3.0)/(2*pi_sq);
  
  double result = Delta_L_sq;

  result *= pow(1+Delta_L_sq, betan);
  result /= 1 + alphan*Delta_L_sq;
  result /= exp(f_no_index(y));

  
  return result;
  
}


/*******************************************************************************************************************************************************
 * 4.7 Delta_H_sq
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::Delta_H_sq(double k){
  
  double n = this->current_n_eff;
  double k_s = this->current_k_non_linear;
  double y = k/k_s;
  
  double result = Delta_H_prime_sq(k);
  
  result /= 1 + mun/y + nun/(y*y);

  
  return result;
  
}


/*******************************************************************************************************************************************************
 * 4.8 Delta_H_prime_sq
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::Delta_H_prime_sq(double k){

  double n = this->current_n_eff;
  double k_s = this->current_k_non_linear;
  double y = k/k_s;
  double C_sm = this->current_C_sm;
  double Om_l = this->universe->rho_L_of_a(this->current_scale);
  double w = this->universe->w_L_of_a(this->current_scale);
  
  double result = 1.0;
  
  result *= pow(y, 3*f_1);
  result *= a_n(n, C_sm, Om_l, w);
  result /= 1 + b_n(n, C_sm, Om_l, w)*pow(y, f_2) + pow(c_n(n, C_sm)*f_3*y , 3.0 - gamma_n(n, C_sm));

  
  return result;
  
}


/*******************************************************************************************************************************************************
 * 4.9 P_NL_at
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::P_NL_at(double k, double e){

  double q = this->Delta_Q_sq(k, e);
  double h = this->Delta_H_sq(k);
  
  return (q+h)/pow(k, 3)*(2*pi*pi);
  
}


/*******************************************************************************************************************************************************
 * 4.10 current_P_NL_at
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::current_P_NL_at(double ln_k){
  
  return interpolate_neville_aitken(ln_k, &this->log_wave_numbers, &this->current_P_NL, constants::order_of_interpolation);
  
}


/*******************************************************************************************************************************************************
 * 4.11 current_P_L_at
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::current_P_L_at(double ln_k){

  return interpolate_neville_aitken(ln_k, &this->log_wave_numbers, &this->current_P_L, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 4.12 P_NL
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/


vector<double> Matter::P_NL(double e){
  
  int n = this->wave_numbers.size();
  vector<double> P(n, 0.0);
  double scale = this->universe->a_at_eta(e);
  double D = return_D_of_eta(e);

  double k;
  double k_s = this->k_NL(0., 3., e);
  double hubble;
  double q;
  double h;
  double n_e;
  double C_sm;
  double Om_m, Om_l;
  double interpol;
  
  vector<double> c_and_n_eff = this->c_and_n_NL(1.0/k_s, e);
  
  this->current_scale = scale;
  hubble = this->universe->H_at_eta(e);
  Om_m = this->universe->rho_m_of_a(scale)*scale*scale/(hubble*hubble);
  Om_l = this->universe->rho_L_of_a(scale)*scale*scale/(hubble*hubble);
  interpol = Om_l/(1.0-Om_m);
  double w = this->universe->w_L_of_a(scale);
  w = interpol*w + (1.0 - interpol)*(-1.0/3.0);
  interpol = -(3.0*w+1.0)/2.0;
  
  this->f_1 = interpol*f_1b(Om_m) + (1.0-interpol)*f_1a(Om_m);
  this->f_2 = interpol*f_2b(Om_m) + (1.0-interpol)*f_2a(Om_m);
  this->f_3 = interpol*f_3b(Om_m) + (1.0-interpol)*f_3a(Om_m);
  
  
  
  this->current_k_non_linear = k_s;
  this->current_n_eff = c_and_n_eff[0]-3.0;
  this->current_C_sm = c_and_n_eff[0]*c_and_n_eff[0]+c_and_n_eff[1];
  n_e = this->current_n_eff;
  C_sm = this->current_C_sm;

  this->mun = mu_n(n_e);
  this->nun = nu_n(n_e);
  this->alphan = alpha_n(n_e, C_sm);
  this->betan = beta_n(n_e, C_sm);

  double z = 1.0/scale - 1.0;
  for(int i = 0; i < n; i++){
    k = this->wave_numbers[i];
    P[i] = this->P_NL_at(k, e);
    if(P[i] < 0.0) P[i] = 0.0;
  }
  
  return P;
  
}





/*******************************************************************************************************************************************************
 * 4.10 P_NL
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

vector<double> Matter::P_L(double e){

  int n = this->wave_numbers.size();

  vector<double> P(n, 0.0);

  for(int i = 0; i < n; i++){
    P[i] = this->Newtonian_linear_power_spectrum(this->wave_numbers[i], e);
  }
  
  return P;
}


/*******************************************************************************************************************************************************
* These functions were for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************
* These functions were for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/
/*******************************************************************************************************************************************************
* These functions were for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/
 
