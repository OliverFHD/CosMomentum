#include <gsl/gsl_multifit.h>



double return_polnomial_value_real(double x, vector<double> *coefficients){
  int order = coefficients->size()-1;
  double polynomial_value = (*coefficients)[order];
  for(int i = order-1; i > -1; i--){
    polynomial_value = x*polynomial_value + (*coefficients)[i];
  }
  
  return polynomial_value;
  
}

complex<double> return_polnomial_value(complex<double> x, vector<double> *coefficients){
  int order = coefficients->size()-1;
  complex<double> polynomial_value = (*coefficients)[order];
  for(int i = order-1; i > -1; i--){
    polynomial_value = x*polynomial_value + (*coefficients)[i];
  }
  
  return polynomial_value;
  
}


vector<double> sub_array(vector<double> * x_sub, vector<double> * x, vector<double>* f){
  
  int n = x_sub->size();
  vector<double> sub_f(n, 0.0);
  for(int i = 0; i < n; i++){
    sub_f[i] = interpolate_neville_aitken((*x_sub)[i] , x, f, 3);
  }
  
  return sub_f;
  
}





complex<double> Newtons_method_complex(complex<double> x_start, complex<double> f_0, vector<double> *f_coefficients, vector<double> *f_prime_coefficients){ 

  double re, im;

  complex<double> x_now = x_start;
  complex<double> x_previous = x_now - 0.1;
  complex<double> x_previous_previous;
  complex<double> f;
  complex<double> f_previous;
  complex<double> f_prime;
  complex<double> dx;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(x_now - x_previous)  > constants::eps_Newton){
    
    f = return_polnomial_value(x_now, f_coefficients)-f_0;
    f_previous = return_polnomial_value(x_previous, f_coefficients)-f_0;
    f_prime = return_polnomial_value(x_now, f_prime_coefficients);

    x_previous_previous = x_previous;
    x_previous = x_now;
 
    dx = f/f_prime;
    if(abs(dx) > 1.0) dx /= abs(dx); 
    if(abs(f_prime) != 0.0)
      x_now = x_previous - dx;
    else
      x_now += 0.1;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(x_now - x_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      x_now += complex<double>(re, im);
    }
  }

  return x_now;

}




complex<double> get_tau_from_secant_method_complex(complex<double> y, complex<double> tau_c, vector<double> *tau_values, vector<double> *G_prime_values){ 

  int n = tau_values->size();
  
  double re, im;

  complex<double> tau_now = tau_c;
  complex<double> tau_previous = tau_now - 0.1;
  complex<double> tau_previous_previous;
  complex<double> f;
  complex<double> f_previous;
  complex<double> f_prime;
  complex<double> Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = y*Newton_interpolation_complex(tau_now, tau_values, G_prime_values)-tau_now;
    f_previous = y*Newton_interpolation_complex(tau_previous, tau_values, G_prime_values)-tau_previous;
    f_prime = (f - f_previous)/(tau_now - tau_previous);

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 1.0;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      tau_now += complex<double>(re, im);
    }
  }

  return tau_now;

}



complex<double> get_tau_from_secant_method_complex(complex<double> y, complex<double> tau_c, vector<double> *G_prime_coefficients){ 
  
  double re, im;

  complex<double> tau_now = tau_c;
  complex<double> tau_previous = tau_now - 0.1;
  complex<double> tau_previous_previous;
  complex<double> f;
  complex<double> f_previous;
  complex<double> f_prime;
  complex<double> Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = y*return_polnomial_value(tau_now, G_prime_coefficients)-tau_now;
    f_previous = y*return_polnomial_value(tau_previous, G_prime_coefficients)-tau_previous;
    f_prime = (f - f_previous)/(tau_now - tau_previous);

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 1.0;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      tau_now += complex<double>(re, im);
    }
  }

  return tau_now;

}



complex<double> get_tau_2D_from_secant_method_complex(double variance, complex<double> lambda, complex<double> tau_c, vector<double> *tau_values, vector<double> *G_prime_values){ 

  int n = tau_values->size();
  
  double re, im;

  complex<double> tau_now = tau_c;
  complex<double> tau_previous = tau_now - 0.1;
  complex<double> tau_previous_previous;
  complex<double> f;
  complex<double> f_previous;
  complex<double> f_prime;
  complex<double> Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = lambda*Newton_interpolation_complex(tau_now, tau_values, G_prime_values)-tau_now/variance;
    f_previous = lambda*Newton_interpolation_complex(tau_previous, tau_values, G_prime_values)-tau_previous/variance;
    f_prime = (f - f_previous)/(tau_now - tau_previous);

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 1.0;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      tau_now += complex<double>(re, im);
    }
  }

  return tau_now;

}



complex<double> get_tau_from_secant_method_complex_Bernardeau_notation(complex<double> lambda, complex<double> tau_c, vector<double> *tau_values, vector<double> *G_prime_values){ 

  int n = tau_values->size();
  
  double re, im;

  complex<double> tau_now = tau_c;
  complex<double> tau_previous = tau_now - 0.1;
  complex<double> tau_previous_previous;
  complex<double> f;
  complex<double> f_previous;
  complex<double> f_prime;
  complex<double> Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = lambda*Newton_interpolation_complex(tau_now, tau_values, G_prime_values)-tau_now;
    f_previous = lambda*Newton_interpolation_complex(tau_previous, tau_values, G_prime_values)-tau_previous;
    f_prime = (f - f_previous)/(tau_now - tau_previous);

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 1.0;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      tau_now += complex<double>(re, im);
    }
  }

  return tau_now;

}


vector<double> return_coefficients(vector<double> *x_in, vector<double> *y_in, int order){
  
  int n = x_in->size();
  double chisq = 0.0;
  vector<double> coefficients(order+1, 0.0);
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;  
  
  X = gsl_matrix_alloc (n, order+1);
  y = gsl_vector_alloc (n);
  c = gsl_vector_alloc (order+1);
  cov = gsl_matrix_alloc (order+1, order+1);
  
  for (int i = 0; i < n; i++){
    for(int j = 0; j <= order; j++)
      gsl_matrix_set (X, i, j, pow((*x_in)[i], j));
    
    gsl_vector_set (y, i, (*y_in)[i]);
  }
  
  {
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, order+1);
    gsl_multifit_linear (X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free (work);
  }
  
  for(int i = 0; i <= order; i++)
    coefficients[i] = gsl_vector_get(c, i);
  
  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (c);
  gsl_matrix_free (cov);
  
  return coefficients;
  
}








complex<double> get_tau_from_secant_method_complex_Bernardeau_notation_2D(complex<double> lambda, complex<double> tau_c, vector<double> *G_prime_coefficients){ 

  double re, im;

  complex<double> tau_now = tau_c;
  complex<double> tau_previous = tau_now - 0.1;
  complex<double> tau_previous_previous;
  complex<double> f;
  complex<double> f_previous;
  complex<double> f_prime;
  complex<double> Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = lambda*return_polnomial_value(tau_now, G_prime_coefficients)-tau_now;
    f_previous = lambda*return_polnomial_value(tau_previous, G_prime_coefficients)-tau_previous;
    f_prime = (f - f_previous)/(tau_now - tau_previous);

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 0.1;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      tau_now += complex<double>(re, im);
    }
  }

  return tau_now;

}

complex<double> get_tau_from_secant_method_complex_Bernardeau_notation_2D(complex<double> lambda, complex<double> tau_c, vector<double> *G_prime_coefficients, vector<double> *G_prime_prime_coefficients){ 

  double re, im;

  complex<double> tau_now = tau_c;
  complex<double> tau_previous = tau_now - 0.1;
  complex<double> tau_previous_previous;
  complex<double> f;
  complex<double> f_prime;
  complex<double> Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = lambda*return_polnomial_value(tau_now, G_prime_coefficients)-tau_now;
    f_prime = lambda*return_polnomial_value(tau_now, G_prime_prime_coefficients) - 1.0;

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 0.1;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      re = double(gen())/double(gen.max());
      im = double(gen())/double(gen.max());
      tau_now += complex<double>(re, im);
    }
  }

  return tau_now;

}


double get_delta_2_from_delta_1_nested(double w, double dmin, double dmax, double fmin, double fmax, double delta_L_1, double V_1, double R_L_1, double R_NL_2, vector<double> *delta_values, vector<double> *F_values, Matter* pointer_to_Matter){
  
  double delta_min = dmin;
  double delta_max = dmax;
  double delta_pivot = 0.5*(dmax+dmin);
  double f_min = fmin;
  double f_max = fmax;
  double f;
  double alpha = delta_L_1/V_1;
  double F;
  double R_L_2;
  double V_12;
  
  int i = 0;
  
  while(abs(delta_max-delta_min) > constants::eps_nested){
    
    F = interpolate_neville_aitken(delta_pivot, delta_values, F_values, constants::order_of_interpolation);
    if(F < -1.0){
      R_L_2 = 1.0/maximal_wave_number;
    }
    else{
      R_L_2 = R_NL_2*sqrt(1.0+F);
    }
    pointer_to_Matter->set_linear_Legendres(2, R_L_2/w);
    //V_12 = pointer_to_Matter->covariance_of_matter_within_R_2D(R_L_1, R_L_2);
    V_12 = pointer_to_Matter->covariance_of_matter_within_R_2D();
    f = delta_pivot - alpha*V_12;
    
    if(f*f_max < 0.0){
      f_min = f;
      delta_min = delta_pivot;
    }
    else{
      f_max = f;
      delta_max = delta_pivot;
    }
    
    delta_pivot = 0.5*(delta_max+delta_min);
    
  }
  
  return delta_pivot;
  
}


double get_delta_2_from_delta_1_final(double w, double d_start, double delta_L_1, double V_1, double R_L_1, double R_NL_2, vector<double> *delta_values, vector<double> *F_values, Matter* pointer_to_Matter){ 
  
  int i = 0;
  
  double d = d_start;
  double d_step;
  double alpha = delta_L_1/V_1;
  double F;
  double F_prime;
  double R_L_2;
  double dR_L_2;
  double V_12;
  double dV_12;
  double f;
  double f_prime;
  
  vector<double> d_values(0, 0.0);
  vector<double> f_values(0, 0.0);
  
  
  F = interpolate_neville_aitken(d, delta_values, F_values, constants::order_of_interpolation);
  F_prime = interpolate_neville_aitken_derivative(d, delta_values, F_values, constants::order_of_interpolation);
  if(F < -1.0){
    R_L_2 = 1.0/maximal_wave_number;
    dR_L_2 = 0.5*F_prime/sqrt(1.0+(*F_values)[0])*R_NL_2;
  }
  else{
    R_L_2 = R_NL_2*sqrt(1.0+F);
    dR_L_2 = 0.5*F_prime/sqrt(1.0+F)*R_NL_2;
  }
  pointer_to_Matter->set_linear_Legendres(3, R_L_2/w);
  
  V_12 = pointer_to_Matter->covariance_of_matter_within_R_2D();
  dV_12 = pointer_to_Matter->dcov2_at()/w;
  f = d - alpha*V_12;
  f_prime =  1.0 - alpha*dV_12*dR_L_2;
  
  d_step = -f/f_prime;
  d_values.push_back(d);
  f_values.push_back(f);
  i = f_values.size();
  
  
  d += d_step;

  
  while(abs(d_step) > constants::eps_Newton){

    
    F = interpolate_neville_aitken(d, delta_values, F_values, 2);
    F_prime = interpolate_neville_aitken_derivative(d, delta_values, F_values, 2);
    if(F < -1.0){
      R_L_2 = 1.0/maximal_wave_number;
      dR_L_2 = 0.0;
    }
    else{
      R_L_2 = R_NL_2*sqrt(1.0+F);
      dR_L_2 = 0.5*F_prime/sqrt(1.0+F)*R_NL_2;
    }
    pointer_to_Matter->set_linear_Legendres(3, R_L_2/w);
    V_12 = pointer_to_Matter->covariance_of_matter_within_R_2D();
    dV_12 = pointer_to_Matter->dcov2_at()/w;
    f = d - alpha*V_12;
    f_prime =  1.0 - alpha*dV_12*dR_L_2;
    d_step = -f/f_prime;
    
    d_values.push_back(d);
    f_values.push_back(f);
    i = f_values.size();
    
    
    if(i > 5){
    /*
    cout << d << setw(20);
    cout << f << setw(20);
    cout << alpha << setw(20);
    cout << dV_12 << setw(20);
    cout << dR_L_2 << setw(20);
    cout << f_prime << '\n';
    */
    
      if(f_values[i-1]*f_values[i-2] < 0.0){
        if(d_values[i-1] < d_values[i-2]){
          return get_delta_2_from_delta_1_nested(w, d_values[i-1], d_values[i-2], f_values[i-1], f_values[i-2], delta_L_1, V_1, R_L_1, R_NL_2, delta_values, F_values, pointer_to_Matter);
        }
        else{
          return get_delta_2_from_delta_1_nested(w, d_values[i-2], d_values[i-1], f_values[i-2], f_values[i-1], delta_L_1, V_1, R_L_1, R_NL_2, delta_values, F_values, pointer_to_Matter);
        }
      }
    }

    
    d += d_step;
    
    
    
        
    if(i > 100){
      cerr << "Error: Newton method stuck (get_delta_2_from_delta_1)\n";
      exit(1);
    }
    
  }

  return d;

}





















double get_lambda_final(double delta, double lambda, vector<double> *G_coefficients, vector<double> *G_prime_coefficients){ 

  double random_step;

  double lambda_now = lambda;
  double lambda_previous = lambda_now - 0.1;
  double lambda_previous_previous;
  double f;
  double f_prime;
  double Dlambda;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(lambda_now - lambda_previous)  > constants::eps_Newton){

    f = return_polnomial_value_real(lambda_now, G_coefficients)-delta;
    f_prime = return_polnomial_value_real(lambda_now, G_prime_coefficients);

    lambda_previous_previous = lambda_previous;
    lambda_previous = lambda_now;
 
    Dlambda = f/f_prime;
    //if(abs(Dlambda) > 1.0) Dlambda /= abs(Dlambda); 
    //if(abs(f_prime) != 0.0)
      lambda_now = lambda_previous - Dlambda;
    //else
    //  lambda_now += 0.1;
    
    // shaking the algorithm if it starts jumping between two points.
    //if(abs(lambda_now - lambda_previous_previous)  < constants::eps_Newton){
    //  random_step = double(gen())/double(gen.max());
    //  lambda_now += random_step;
    //}
  }

  return lambda_now;

}




double get_lambda_from_tau(double tau, double lambda, vector<double> *tau_coefficients, vector<double> *tau_prime_coefficients){ 

  double random_step;

  double lambda_now = lambda;
  double lambda_previous = lambda_now - 0.1;
  double lambda_previous_previous;
  double f;
  double f_prime;
  double Dlambda;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(lambda_now - lambda_previous)  > constants::eps_Newton){

    f = return_polnomial_value_real(lambda_now, tau_coefficients)-tau;
    f_prime = return_polnomial_value_real(lambda_now, tau_prime_coefficients);

    lambda_previous_previous = lambda_previous;
    lambda_previous = lambda_now;
 
    Dlambda = f/f_prime;
    if(abs(Dlambda) > 1.0) Dlambda /= abs(Dlambda); 
    if(abs(f_prime) != 0.0)
      lambda_now = lambda_previous - Dlambda;
    else
      lambda_now += 0.1;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(lambda_now - lambda_previous_previous)  < constants::eps_Newton){
      random_step = double(gen())/double(gen.max());
      lambda_now += random_step;
    }
  }

  return lambda_now;

}


double get_tau_final(double delta, double tau, vector<double> *G_coefficients, vector<double> *G_prime_coefficients){ 

  double random_step;

  double tau_now = tau;
  double tau_previous = tau_now - 0.1;
  double tau_previous_previous;
  double f;
  double f_prime;
  double Dtau;
  
  /*to shake the algorithm when it becomes stiff.*/
  mt19937 gen(2);

  while(abs(tau_now - tau_previous)  > constants::eps_Newton){

    f = return_polnomial_value_real(tau_now, G_coefficients)-delta;
    f_prime = return_polnomial_value_real(tau_now, G_prime_coefficients);

    tau_previous_previous = tau_previous;
    tau_previous = tau_now;
 
    Dtau = f/f_prime;
    if(abs(Dtau) > 1.0) Dtau /= abs(Dtau); 
    if(abs(f_prime) != 0.0)
      tau_now = tau_previous - Dtau;
    else
      tau_now += 0.1;
    
    // shaking the algorithm if it starts jumping between two points.
    if(abs(tau_now - tau_previous_previous)  < constants::eps_Newton){
      random_step = double(gen())/double(gen.max());
      tau_now += random_step;
    }
  }

  return tau_now;
  
}



vector<double> get_tau_coefficients(vector<double> *phi_coefficients){
  
  int n = phi_coefficients->size();
  
  vector<double> tau_coefficients(n-1, 0.0);
  vector<double> tau_squared_coefficients(n, 0.0);
  
  for(int i = 2; i < n; i++){
    tau_squared_coefficients[i] = 2.0*(*phi_coefficients)[i]*double(i-1);
  }
  
  //positive sign here, because tau monotonically increases with y:
  tau_coefficients[1] = sqrt(tau_squared_coefficients[2]);
  cout << 1 << setw(20) << tau_coefficients[1] << '\n';
  for(int i = 2; i < n-1; i++){
    tau_coefficients[i] = tau_squared_coefficients[i+1];
    for(int j = 2; j < i; j++){
      tau_coefficients[i] -= tau_coefficients[j]*tau_coefficients[i+1-j];
    }
    tau_coefficients[i] /= 2.0*tau_coefficients[1];
    cout << i << setw(20) << tau_coefficients[i] << '\n';
  }
  
  return tau_coefficients;
  
}

