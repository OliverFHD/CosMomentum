

void Matter::set_P_today(){
  
  double eta = this->universe->eta_at_a(1.0);
  this->P_L_today = this->P_L(eta);
  this->P_NL_today = this->P_NL(eta);
  
}



double Matter::P_L_today_at(double ln_k){

  return interpolate_neville_aitken(ln_k, &this->log_wave_numbers, &this->P_L_today, constants::order_of_interpolation);
  
}


