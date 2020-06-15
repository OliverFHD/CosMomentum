
using namespace std;

/*******************************************************************************************************************************************************
 * Vector Utils
 * 
*******************************************************************************************************************************************************/



void log_binning(double x_min, double x_max, int n, vector<double> *x){
  
  (*x) = vector<double>(n+1, 0.0);
  
  for(int i = 0; i < n+1; i++){
    (*x)[i] = x_min*pow(x_max/x_min, double(i)/double(n));
  }
  
}

int find_index(double x0, vector<double> *x){

  int i1 = 0;
  int i2 = x->size();
  int pivot;

  while(i2-i1 > 1){
    pivot = (i2+i1)/2;
    if((*x)[pivot]<=x0)
      i1 = pivot;
    else
      i2 = pivot;
  }
  
  if(i1 == (*x).size()-1)
    i1--;

  return i1;

}



/*******************************************************************************************************************************************************
 * Neville-Aitken Interpolation
 * 
*******************************************************************************************************************************************************/

double neville_aitken(double x0, vector<double> *x, vector<double> *f){
  
  int n = (*x).size();
  
  if(n == 2){
    
    double x1 = (*x)[0];
    double f1 = (*f)[0];
    
    return f1 + (x0-x1)*((*f)[1]-f1)/((*x)[1]-x1);
  }
  else{
    int n_m1 = n-1;
  
    double x1 = (*x)[0];
    double x2 = (*x)[n_m1];
  
    vector<double> x_new_1(n_m1);
    vector<double> x_new_2(n_m1);
    vector<double> f_new_1(n_m1);
    vector<double> f_new_2(n_m1);
  
    for(int i = 0; i < n_m1; i++){
       x_new_1[i] = (*x)[i];
       x_new_2[i] = (*x)[i+1];
     
       f_new_1[i] = (*f)[i];
       f_new_2[i] = (*f)[i+1];
    }
  
    return ((x0-x1)*neville_aitken(x0, &x_new_2, &f_new_2)-(x0-x2)*neville_aitken(x0, &x_new_1, &f_new_1))/(x2-x1);
  }
  
  return -1;
  
}


double neville_aitken_more_clever(int index, int length, double x0, vector<double> *x, vector<double> *f){
  
  if(length == 2){
    
    double x1 = (*x)[index];
    double f1 = (*f)[index];
    
    return f1 + (x0-x1)*((*f)[index+1]-f1)/((*x)[index+1]-x1);
  }
  else{
    int n_m1 = length-1;
  
    double x1 = (*x)[index];
    double x2 = (*x)[index+n_m1];
  
    /*vector<double> x_new_1(n_m1);
    vector<double> x_new_2(n_m1);
    vector<double> f_new_1(n_m1);
    vector<double> f_new_2(n_m1);
  
    for(int i = 0; i < n_m1; i++){
       x_new_1[i] = (*x)[i];
       x_new_2[i] = (*x)[i+1];
     
       f_new_1[i] = (*f)[i];
       f_new_2[i] = (*f)[i+1];
    }*/
  
    return ((x0-x1)*neville_aitken_more_clever(index+1, n_m1, x0, x, f)-(x0-x2)*neville_aitken_more_clever(index, n_m1, x0, x, f))/(x2-x1);
  }
  
  return -1;
  
}

double neville_aitken_derivative(double x0, vector<double> *x, vector<double> *f){
  
  int n = (*x).size();
  
  if(n == 2){
    
    double x1 = (*x)[0];
    double f1 = (*f)[0];
    
    return ((*f)[1]-f1)/((*x)[1]-x1);
  }
  else{
    int n_m1 = n-1;
  
    double x1 = (*x)[0];
    double x2 = (*x)[n_m1];
  
    vector<double> x_new_1(n_m1);
    vector<double> x_new_2(n_m1);
    vector<double> f_new_1(n_m1);
    vector<double> f_new_2(n_m1);
  
    for(int i = 0; i < n_m1; i++){
       x_new_1[i] = (*x)[i];
       x_new_2[i] = (*x)[i+1];
     
       f_new_1[i] = (*f)[i];
       f_new_2[i] = (*f)[i+1];
    }
  
    return (neville_aitken(x0, &x_new_2, &f_new_2) + (x0-x1)*neville_aitken_derivative(x0, &x_new_2, &f_new_2)-neville_aitken(x0, &x_new_1, &f_new_1)-(x0-x2)*neville_aitken_derivative(x0, &x_new_1, &f_new_1))/(x2-x1);
  }
  
  return -1;
  
}

double neville_aitken_derivative_more_clever(int index, int length, double x0, vector<double> *x, vector<double> *f){
    
  if(length == 2){
    
    double x1 = (*x)[index];
    double f1 = (*f)[index];
    
    return ((*f)[index+1]-f1)/((*x)[index+1]-x1);
  }
  else{
    int n_m1 = length-1;
  
    double x1 = (*x)[index];
    double x2 = (*x)[index+n_m1];
  
    return (neville_aitken_more_clever(index+1, n_m1, x0, x, f) + (x0-x1)*neville_aitken_derivative_more_clever(index+1, n_m1, x0, x, f)-neville_aitken_more_clever(index, n_m1, x0, x, f)-(x0-x2)*neville_aitken_derivative_more_clever(index, n_m1, x0, x, f))/(x2-x1);
  }
  
  return -1;
  
}

double neville_aitken_2nd_derivative_more_clever(int index, int length, double x0, vector<double> *x, vector<double> *f){
  
  if(length < 3){
    return 0.0;
  }
  else if(length == 3){
    
    double x1 = (*x)[index];
    double x2 = (*x)[index+1];
    double x3 = (*x)[index+2];
    double f1 = (*f)[index];
    double f2 = (*f)[index+1];
    double f3 = (*f)[index+2];
    double h1 = x2-x1;
    double h2 = x3-x2;
    
    return 2.0*((f1*h2+f3*h1)/(h1+h2)-f2)/(h1*h2);
  }
  else{
    int n_m1 = length-1;
  
    double x1 = (*x)[index];
    double x2 = (*x)[index+n_m1];
      
    double value = 2.0*neville_aitken_derivative_more_clever(index+1, n_m1, x0, x, f) + (x0-x1)*neville_aitken_2nd_derivative_more_clever(index+1, n_m1, x0, x, f);
    value -= 2.0*neville_aitken_derivative_more_clever(index, n_m1, x0, x, f) + (x0-x2)*neville_aitken_2nd_derivative_more_clever(index, n_m1, x0, x, f);
    value /= (x2-x1);
    return value;
    
  }
  
  return -1;
  
}


double interpolate_neville_aitken(double x0, vector<double> *x, vector<double> *f, int order){
  
  int index = find_index(x0, x);
  int n = (*x).size();
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order = 1;
  
  vector<double> x_new(order+1);
  vector<double> f_new(order+1);
  
  
  index -= (order+1)/2;
  if(index < 0)
    index = 0;
  if(index + order > n - 1)
    index = n - order - 1;
  
  for(int i = 0; i <= order; i++){
    x_new[i] = (*x)[index+i];
    f_new[i] = (*f)[index+i];
  }
  
  return neville_aitken_more_clever(0, order+1, x0, &x_new, &f_new);
  //return neville_aitken(x0, &x_new, &f_new);
  
}

double interpolate_neville_aitken_derivative(double x0, vector<double> *x, vector<double> *f, int order){
  
  int index = find_index(x0, x);
  int n = (*x).size();
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order = 1;
		
  
  vector<double> x_new(order+1);
  vector<double> f_new(order+1);
  
  
  index -= (order+1)/2;
  if(index < 0)
    index = 0;
  if(index + order > n - 1)
    index = n - order - 1;
  
  for(int i = 0; i <= order; i++){
    x_new[i] = (*x)[index+i];
    f_new[i] = (*f)[index+i];
  }
  
  return neville_aitken_derivative_more_clever(0, order+1, x0, &x_new, &f_new);
  //return neville_aitken_derivative(x0, &x_new, &f_new);
  
}

double interpolate_neville_aitken_2nd_derivative(double x0, vector<double> *x, vector<double> *f, int order){
  
  int index = find_index(x0, x);
  int n = (*x).size();
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order = 1;
		
  
  vector<double> x_new(order+1);
  vector<double> f_new(order+1);
  
  
  index -= (order+1)/2;
  if(index < 0)
    index = 0;
  if(index + order > n - 1)
    index = n - order - 1;
  
  for(int i = 0; i <= order; i++){
    x_new[i] = (*x)[index+i];
    f_new[i] = (*f)[index+i];
  }
  
  return neville_aitken_2nd_derivative_more_clever(0, order+1, x0, &x_new, &f_new);
  
}

double interpolate_neville_aitken_grid(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x, int order_t, int order_x){
  
  int index = find_index(t0, t);
  int n = (*t).size();
	
	
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order_x = 1;
	if(t0 > (*t)[n-1] || t0 < (*t)[0])
		order_t = 1;
  
  vector<double> t_new(order_t+1, 0.0);
  vector<double> P_new(order_t+1, 0.0);
  
  index -= (order_t+1)/2;
  if(index < 0)
    index = 0;
  if(index + order_t > n - 1)
    index = n - order_t - 1;
  
  for(int i = 0; i <= order_t; i++){
    t_new[i] = (*t)[i+index];
    P_new[i] = interpolate_neville_aitken(x0, x, &(*f_of_t_x)[i+index], order_x);
  }
  
  //return neville_aitken(t0, &t_new, &P_new);
  return neville_aitken_more_clever(0, order_t+1, t0, &t_new, &P_new);
  
}

double interpolate_neville_aitken_dgrid_dt(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x, int order_t, int order_x){
  
  int index = find_index(t0, t);
  int n = (*t).size();
	
	
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order_x = 1;
	if(t0 > (*t)[n-1] || t0 < (*t)[0])
		order_t = 1;
  
  vector<double> t_new(order_t+1, 0.0);
  vector<double> P_new(order_t+1, 0.0);
  
  index -= (order_t+1)/2;
  if(index < 0)
    index = 0;
  if(index + order_t > n - 1)
    index = n - order_t - 1;
  
  for(int i = 0; i <= order_t; i++){
    t_new[i] = (*t)[i+index];
    P_new[i] = interpolate_neville_aitken(x0, x, &(*f_of_t_x)[i+index], order_x);
  }
  
  //return neville_aitken(t0, &t_new, &P_new);
  return neville_aitken_derivative_more_clever(0, order_t+1, t0, &t_new, &P_new);
  
}

double interpolate_neville_aitken_dgrid_dx(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x, int order_t, int order_x){
  
  int index = find_index(t0, t);
  int n = (*t).size();
	
	
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order_x = 1;
	if(t0 > (*t)[n-1] || t0 < (*t)[0])
		order_t = 1;
  
  vector<double> t_new(order_t+1, 0.0);
  vector<double> P_new(order_t+1, 0.0);
  
  index -= (order_t+1)/2;
  if(index < 0)
    index = 0;
  if(index + order_t > n - 1)
    index = n - order_t - 1;
  
  for(int i = 0; i <= order_t; i++){
    t_new[i] = (*t)[i+index];
    P_new[i] = interpolate_neville_aitken_derivative(x0, x, &(*f_of_t_x)[i+index], order_x);
  }
  
  //return neville_aitken(t0, &t_new, &P_new);
  return neville_aitken_more_clever(0, order_t+1, t0, &t_new, &P_new);
  
}







/*
complex<double> neville_aitken_complex(complex<double> x0, vector<double> *x, vector<double> *f){
  
  int n = (*x).size();
  if(n == 2){
    
    complex<double> x1 = complex<double>((*x)[0], 0.0);
    complex<double> f1 = complex<double>((*f)[0], 0.0);
    
    return f1 + (x0-x1)*(complex<double>((*f)[1], 0.0)-f1)/(complex<double>((*x)[1], 0.0)-x1);
  }
  else{
    int n_m1 = n-1;
  
    complex<double> x1 = complex<double>((*x)[0], 0.0);
    complex<double> x2 = complex<double>((*x)[n_m1], 0.0);
  
    vector<double> x_new_1(n_m1, 0.0);
    vector<double> x_new_2(n_m1, 0.0);
    vector<double> f_new_1(n_m1, 0.0);
    vector<double> f_new_2(n_m1, 0.0);
  
    for(int i = 0; i < n_m1; i++){
       x_new_1[i] = (*x)[i];
       x_new_2[i] = (*x)[i+1];
     
       f_new_1[i] = (*f)[i];
       f_new_2[i] = (*f)[i+1];
    }
  
    return ((x0-x1)*neville_aitken_complex(x0, &x_new_2, &f_new_2)-(x0-x2)*neville_aitken_complex(x0, &x_new_1, &f_new_1))/(x2-x1);
  }
  
  return complex<double>(-1.0, 0.0);
  
}
*/



complex<double> Newton_interpolation_complex(complex<double> x0, vector<double> *x, vector<double> *f){
  
  int n = (*x).size();
  complex<double> p_of_x(0.0, 0.0);
  complex<double> product;
  
  for(int i = 0; i < n; i++){
    product = complex<double>(1.0, 0.0);
    for(int j = 0; j < n; j++){
      if(j !=i){
	product *= (x0-(*x)[j])/((*x)[i]-(*x)[j]);
      }
    }
    p_of_x += product*(*f)[i];
  }
  
  return p_of_x;
  
}


double Newton_interpolation(double x0, vector<double> *x, vector<double> *f){
  
  int n = (*x).size();
  double p_of_x = 0.0;
  double product;
  
  for(int i = 0; i < n; i++){
    product = 1.0;
    for(int j = 0; j < n; j++){
      if(j !=i){
				product *= (x0-(*x)[j])/((*x)[i]-(*x)[j]);
      }
    }
    p_of_x += product*(*f)[i];
  }
  
  return p_of_x;
  
}

double interpolate_Newton(double x0, vector<double> *x, vector<double> *f, int order){
  
  int index = find_index(x0, x);
  int n = (*x).size();
	if(x0 > (*x)[n-1] || x0 < (*x)[0])
		order = 1;
  
  vector<double> x_new(order+1);
  vector<double> f_new(order+1);
  
  
  index -= (order+1)/2;
  if(index < 0)
    index = 0;
  if(index + order > n - 1)
    index = n - order - 1;
  
  for(int i = 0; i <= order; i++){
    x_new[i] = (*x)[index+i];
    f_new[i] = (*f)[index+i];
  }
  
  return Newton_interpolation(x0, &x_new, &f_new);
  
}

