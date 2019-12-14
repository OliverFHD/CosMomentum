
using namespace std;


enum STATE_OF_EARLY_UNIVERSE {UNDEFINED, MATTER_DOMINATED, RADIATION_DOMINATED, LAMBDA_DOMINATED};


void Universe::set_initial_conditions(){
  
  error_handling::general_warning("WARNING in Universe::set_initial_conditions(): initial conditions module still assumes flat universe!");
  
  STATE_OF_EARLY_UNIVERSE early_state = UNDEFINED;
  
  if(this->cosmology.Omega_r != 0.0){
    early_state = RADIATION_DOMINATED;
  }
  else if(this->cosmology.Omega_m != 0.0){
    early_state = MATTER_DOMINATED;
  }
  else if(this->cosmology.Omega_L != 0.0){
    early_state = LAMBDA_DOMINATED;
  }
  else{
    error_handling::general_error_message("Error: universe has no matter content!");
  }
  
  double a_i = this->a_initial; // early time scale factor
  double t_i; // early time physical time
  double e_i; // early time conformal time
  double H_i; // early time conformal expansion rate
  double H_prime_i; // early time derivative of conformal expansion rate
  
  switch(early_state){
    case RADIATION_DOMINATED:
      if(this->cosmology.Omega_m != 0.0 || this->cosmology.Omega_L != 0.0)
        a_i = min( min(0.0005*this->cosmology.Omega_r/this->cosmology.Omega_m, pow(0.0005*this->cosmology.Omega_r/this->cosmology.Omega_L, 0.25) ), this->a_initial); // go to a time when other components where at most 1/1000 of total energy density
      expansion_in_flat_radiation_dominated_universe(a_i, &t_i, &e_i, &H_i, &H_prime_i);
      break;
    case MATTER_DOMINATED:
      if(this->cosmology.Omega_L != 0.0)
        a_i = min( 0.00001*pow(this->cosmology.Omega_m/this->cosmology.Omega_L, 1.0/3.0), this->a_initial); // go to a time when other components where at most 1/1000 of total energy density
      expansion_in_flat_matter_dominated_universe(a_i, this->cosmology.Omega_m, &t_i, &e_i, &H_i, &H_prime_i);
      break;
    case LAMBDA_DOMINATED:
      expansion_in_flat_Lambda_dominated_universe(a_i, &t_i, &e_i, &H_i, &H_prime_i);
      break;
    case UNDEFINED:
      error_handling::general_error_message("Universe should have non-zero energy density content in order to compute initial conditions!");
      break;
  }
  
  if(a_i >= this->a_initial){
    this->t_initial = t_i;
    this->eta_initial = e_i;
    this->H_initial = H_i;
  }
  else{
    double a_f = this->a_initial;
    double y[3] = { e_i, H_i, t_i};
    
    integration_parameters_Universe params;
    params.pointer_to_Universe = this;
    integration_parameters_Universe * pointer_to_params = &params;
    gsl_odeiv2_system sys = {scale_factor_gsl, scale_factor_gsl_jac, 3, (void *) pointer_to_params};
    
    double hstart = (a_f-a_i)*constants::gsl_hstart_relative;
    double eps_absolute = constants::gsl_eps_relative*std::max(e_i, std::max(H_i, t_i));
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, eps_absolute, constants::gsl_eps_relative);
    
    int status = gsl_odeiv2_driver_apply(d, &a_i, a_f, y);
    
    this->eta_initial = y[0];
    this->H_initial = y[1];
    this->t_initial = y[2];
    
    gsl_odeiv2_driver_free(d);
  }
  
}
