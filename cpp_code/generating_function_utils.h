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
  for(int i = 2; i < n-1; i++){
    tau_coefficients[i] = tau_squared_coefficients[i+1];
    for(int j = 2; j < i; j++){
      tau_coefficients[i] -= tau_coefficients[j]*tau_coefficients[i+1-j];
    }
    tau_coefficients[i] /= 2.0*tau_coefficients[1];
  }
  
  return tau_coefficients;
  
}

