
using namespace std;

void FlatInhomogeneousUniverseLCDM::set_initial_conditions_for_growth(){
  
  error_handling::general_warning("WARNING in FlatInhomogeneousUniverseLCDM::set_initial_conditions(): growth factor is computed for matter density fluctuations and all other components are assumed to be homogeneous!");
  
  STATE_OF_EARLY_UNIVERSE early_state = UNDEFINED;
  
  double Omega_m = this->return_Omega_m();
  double Omega_r = this->return_Omega_r();
  double Omega_L = this->return_Omega_L();
  
  
  if(Omega_r != 0.0){
    early_state = RADIATION_DOMINATED;
  }
  else if(Omega_m != 0.0){
    early_state = MATTER_DOMINATED;
  }
  else{
    error_handling::general_error_message("ERROR in FlatInhomogeneousUniverseLCDM::set_initial_conditions(): universe has no energy density content besides Lambda!");
  }
  
  double a_i = this->return_a_initial(); // early time scale factor
  double t_i; // early time physical time
  double e_i; // early time conformal time
  double H_i; // early time conformal expansion rate
  double H_prime_i; // early time derivative of conformal expansion rate
  double D_i; // early time conformal expansion rate
  double D_prime_i; // early time derivative of conformal expansion rate
  
  switch(early_state){
    case RADIATION_DOMINATED:
      if(Omega_m != 0.0 || Omega_L != 0.0){
        a_i = min(0.0005*Omega_r/Omega_m, pow(0.0005*Omega_r/Omega_L, 0.25) ); // go to a time when other components where at most 1/1000 of total energy density
        a_i = min(a_i, this->return_a_initial());
      }
      this->expansion_in_flat_radiation_dominated_universe(a_i, &t_i, &e_i, &H_i, &H_prime_i);
      D_i = a_i;
      D_prime_i = 0.0;
      break;
    case MATTER_DOMINATED:
      if(Omega_L != 0.0){
        a_i = 0.0001*pow(Omega_m/Omega_L, 1.0/3.0); // go to a time when other components where at most 1/1000 of total energy density
        a_i = min(a_i, this->return_a_initial());
      }
      this->expansion_in_flat_matter_dominated_universe(a_i, Omega_m, &t_i, &e_i, &H_i, &H_prime_i);
      D_i = a_i;
      D_prime_i = a_i*H_i;
      break;
    case LAMBDA_DOMINATED:
      error_handling::general_error_message("ERROR in FlatInhomogeneousUniverseLCDM::set_initial_conditions(): universe has no energy density content besides Lambda!");
      break;
    case UNDEFINED:
      error_handling::general_error_message("ERROR in FlatInhomogeneousUniverseLCDM::set_initial_conditions(): universe has no energy density content besides Lambda!");
      break;
  }
  
  if(a_i < this->return_a_initial()){
    double e_f = this->return_eta_initial();
    double y[2] = { D_i, D_prime_i};
    
    integration_parameters params;
    params.pointer_to_Universe = this;
    params.Omega_m = Omega_m;
    integration_parameters * pointer_to_params = &params;
    
    gsl_odeiv2_system sys = {growth_factor_gsl, growth_factor_gsl_jac, 2, (void *) pointer_to_params};
    
    double hstart = 1.0*constants::gsl_hstart_relative;
    double eps_absolute = 1.0*constants::gsl_eps_relative;
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, eps_absolute, constants::gsl_eps_relative);
    
    int status = gsl_odeiv2_driver_apply(d, &e_i, e_f, y);
    
    D_i = y[0];
    D_prime_i = y[1];
    
    gsl_odeiv2_driver_free(d);
  }
  
  this->D_initial = D_i;
  this->D_prime_initial = D_prime_i;
  
}
