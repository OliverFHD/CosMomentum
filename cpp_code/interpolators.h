
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cblas.h>

using namespace std;

/*******************************************************************************************************************************************************
 * Vector Utils
 * 
*******************************************************************************************************************************************************/




void invert_matrix(vector<vector<double> > *matrix, vector<vector<double> > *inverse){
	
	int n_data = matrix->size();
	
	(*inverse) = vector<vector<double> >(n_data, vector<double>(n_data, 0.0));
	
	gsl_matrix * m = gsl_matrix_alloc (n_data, n_data);
	gsl_matrix * mAux = gsl_matrix_alloc (n_data, n_data);
	gsl_matrix * mInverse = gsl_matrix_alloc (n_data, n_data);
	gsl_permutation * perm1 = gsl_permutation_alloc (n_data);
	int s1;
	
	for(int i = 0; i < n_data; i++){
		for(int j = 0; j < n_data; j++){
			gsl_matrix_set(m, i, j, (*matrix)[i][j]);
		}
	}
	
	gsl_matrix_memcpy (mAux, m);
	gsl_linalg_LU_decomp (mAux, perm1, &s1);
	gsl_linalg_LU_invert (mAux, perm1, mInverse);
	
	
	for(int i = 0; i < n_data; i++){
		for(int j = 0; j < n_data; j++){
			(*inverse)[i][j] = gsl_matrix_get(mInverse, i, j);
		}
	}
	
	
	gsl_matrix_free(m);
	gsl_matrix_free(mAux);
	gsl_matrix_free(mInverse);
	gsl_permutation_free(perm1);
	
}


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
  if(order == 0){
    return (*f)[index];
  }
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

double interpolate_neville_aitken_grid(double t0, double x0, vector<double> * t, vector<double> * x, vector<vector<double> > * f_of_t_x, int order_t, int order_x){
  
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





double interpolate_grid_bilinear(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x){
  
  int n_t = (*t).size();
  int n_x = (*x).size();
  int index_t = find_index(t0, t);
  int index_x = find_index(x0, x);
  if(index_t >= n_t - 1) index_t = n_t - 2;
  if(index_x >= n_x - 1) index_x = n_x - 2;
  
  /*
   * (index_t, index_x)---> p00  p01
   *                        p10  
   * ISSUE: this treatment is not isotropic.
   */
  
  double dt = t0 - (*t)[index_t];
  double dx = x0 - (*x)[index_x];
  double dt1 = (*t)[index_t+1] - (*t)[index_t];
  double dx1 = (*x)[index_x+1] - (*x)[index_x];
  
  /*
   * f - f00 = df/dt*dt + df/dx*dx
   */
  
  double f00 = (*f_of_t_x)[index_t][index_x];
  double f10 = (*f_of_t_x)[index_t+1][index_x];
  double f01 = (*f_of_t_x)[index_t][index_x+1];
  
  return f00 + dt*(f10-f00)/dt1 + dx*(f01-f00)/dx1;
  
}


double interpolate_grid_biquadratic(double t0, double x0, vector<double> * t, vector<double> * x, vector<vector<double> > * f_of_t_x){
  
  int n_t = (*t).size();
  int n_x = (*x).size();
  int index_t = find_index(t0, t);
  int index_x = find_index(x0, x);
  if(index_t >= n_t - 2) index_t = n_t - 3;
  if(index_x >= n_x - 2) index_x = n_x - 3;
  
  /*
   * (index_t, index_x)---> p00  p01  p02
   *                        p10  p11
   *                        p20
   * ISSUE: this treatment is not isotropic.
   * However, OF didn't want to go for bi-cubic interpolation,
   * because it is so slow.
   */
  
  double dt = t0 - (*t)[index_t];
  double dx = x0 - (*x)[index_x];
  double dt1 = (*t)[index_t+1] - (*t)[index_t];
  double dt2 = (*t)[index_t+2] - (*t)[index_t];
  double dx1 = (*x)[index_x+1] - (*x)[index_x];
  double dx2 = (*x)[index_x+2] - (*x)[index_x];
  
  /*
   * f - f00 = gt*dt + gx*dx + Htt*dt^2 + 2*Htx*dt*dx + Hxx*dx^2
   * --> need 5D matrix
   */
  
  double f00 = (*f_of_t_x)[index_t][index_x];
  vector<double> f_reduced(5, 0.0);
  vector<double> coeffs(5, 0.0);
  /* f10_reduced: */ f_reduced[0] = (*f_of_t_x)[index_t+1][index_x]-f00;
  /* f20_reduced: */ f_reduced[1] = (*f_of_t_x)[index_t+2][index_x]-f00;
  /* f01_reduced: */ f_reduced[2] = (*f_of_t_x)[index_t][index_x+1]-f00;
  /* f11_reduced: */ f_reduced[3] = (*f_of_t_x)[index_t+1][index_x+1]-f00;
  /* f02_reduced: */ f_reduced[4] = (*f_of_t_x)[index_t][index_x+2]-f00;
  vector<vector<double> > A(5, vector<double>(5, 0.0));
  vector<vector<double> > A_inverse;
  A[0][0] = dt1;
  //A[0][1] = 0.0;
  A[0][2] = dt1*dt1;
  //A[0][3] = 0.0;
  //A[0][4] = 0.0;
  
  A[1][0] = dt2;
  //A[1][1] = 0.0;
  A[1][2] = dt2*dt2;
  //A[1][3] = 0.0;
  //A[1][4] = 0.0;
  
  //A[2][0] = 0.0;
  A[2][1] = dx1;
  //A[2][2] = 0.0;
  //A[2][3] = 0.0;
  A[2][4] = dx1*dx1;
  
  A[3][0] = dt1;
  A[3][1] = dx1;
  A[3][2] = dt1*dt1;
  A[3][3] = dt1*dx1;
  A[3][4] = dx1*dx1;
  
  //A[4][0] = 0.0;
  A[4][1] = dx2;
  //A[4][2] = 0.0;
  //A[4][3] = 0.0;
  A[4][4] = dx2*dx2;
  
  
  invert_matrix(&A, &A_inverse);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      coeffs[i] += A_inverse[i][j]*f_reduced[j];
    }
  }
  
  // coeffs[3] = 2.0*Hxt
  return f00 + coeffs[0]*dt + coeffs[1]*dx + coeffs[2]*dt*dt + coeffs[3]*dt*dx + coeffs[4]*dx*dx;
  
}


double interpolate_grid_biquadratic_with_mask(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x, vector<vector<int> > * mask){
  
  int n_t = (*t).size();
  int n_x = (*x).size();
  int index_t = find_index(t0, t);
  int index_x = find_index(x0, x);
  if(index_t >= n_t - 2) index_t = n_t - 3;
  if(index_x >= n_x - 2) index_x = n_x - 3;
  
  /*
   * (index_t, index_x)---> p00  p01  p02
   *                        p10  p11
   *                        p20
   * ISSUE: this treatment is not isotropic.
   * However, OF didn't want to go for bi-cubic interpolation,
   * because it is so slow.
   */
  
  double dt = t0 - (*t)[index_t];
  double dx = x0 - (*x)[index_x];
  double dt1 = (*t)[index_t+1] - (*t)[index_t];
  double dt2 = (*t)[index_t+2] - (*t)[index_t];
  double dx1 = (*x)[index_x+1] - (*x)[index_x];
  double dx2 = (*x)[index_x+2] - (*x)[index_x];
  
  int masked = (*mask)[index_t][index_x];
  masked *= (*mask)[index_t+1][index_x];
  masked *= (*mask)[index_t+2][index_x];
  masked *= (*mask)[index_t][index_x+1];
  masked *= (*mask)[index_t+1][index_x+1];
  masked *= (*mask)[index_t][index_x+2];
  
  if(masked == 0)
    return 0.0;
  
  /*
   * f - f00 = gt*dt + gx*dx + Htt*dt^2 + 2*Htx*dt*dx + Hxx*dx^2
   * --> need 5D matrix
   */
  
  double f00 = (*f_of_t_x)[index_t][index_x];
  vector<double> f_reduced(5, 0.0);
  vector<double> coeffs(5, 0.0);
  /* f10_reduced: */ f_reduced[0] = (*f_of_t_x)[index_t+1][index_x]-f00;
  /* f20_reduced: */ f_reduced[1] = (*f_of_t_x)[index_t+2][index_x]-f00;
  /* f01_reduced: */ f_reduced[2] = (*f_of_t_x)[index_t][index_x+1]-f00;
  /* f11_reduced: */ f_reduced[3] = (*f_of_t_x)[index_t+1][index_x+1]-f00;
  /* f02_reduced: */ f_reduced[4] = (*f_of_t_x)[index_t][index_x+2]-f00;
  vector<vector<double> > A(5, vector<double>(5, 0.0));
  vector<vector<double> > A_inverse;
  A[0][0] = dt1;
  //A[0][1] = 0.0;
  A[0][2] = dt1*dt1;
  //A[0][3] = 0.0;
  //A[0][4] = 0.0;
  
  A[1][0] = dt2;
  //A[1][1] = 0.0;
  A[1][2] = dt2*dt2;
  //A[1][3] = 0.0;
  //A[1][4] = 0.0;
  
  //A[2][0] = 0.0;
  A[2][1] = dx1;
  //A[2][2] = 0.0;
  //A[2][3] = 0.0;
  A[2][4] = dx1*dx1;
  
  A[3][0] = dt1;
  A[3][1] = dx1;
  A[3][2] = dt1*dt1;
  A[3][3] = dt1*dx1;
  A[3][4] = dx1*dx1;
  
  //A[4][0] = 0.0;
  A[4][1] = dx2;
  //A[4][2] = 0.0;
  //A[4][3] = 0.0;
  A[4][4] = dx2*dx2;
  
  
  invert_matrix(&A, &A_inverse);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      coeffs[i] += A_inverse[i][j]*f_reduced[j];
    }
  }
  
  // coeffs[3] = 2.0*Hxt
  return f00 + coeffs[0]*dt + coeffs[1]*dx + coeffs[2]*dt*dt + coeffs[3]*dt*dx + coeffs[4]*dx*dx;
  
}


double interpolate_grid_biquadratic_dt(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x){
  
  int n_t = (*t).size();
  int n_x = (*x).size();
  int index_t = find_index(t0, t);
  int index_x = find_index(x0, x);
  if(index_t >= n_t - 2) index_t = n_t - 3;
  if(index_x >= n_x - 2) index_x = n_x - 3;
  
  /*
   * (index_t, index_x)---> p00  p01  p02
   *                        p10  p11
   *                        p20
   * ISSUE: this treatment is not isotropic.
   * However, OF didn't want to go for bi-cubic interpolation,
   * because it is so slow.
   */
  
  double dt = t0 - (*t)[index_t];
  double dx = x0 - (*x)[index_x];
  double dt1 = (*t)[index_t+1] - (*t)[index_t];
  double dt2 = (*t)[index_t+2] - (*t)[index_t];
  double dx1 = (*x)[index_x+1] - (*x)[index_x];
  double dx2 = (*x)[index_x+2] - (*x)[index_x];
  
  /*
   * f - f00 = gt*dt + gx*dx + Htt*dt^2 + 2*Htx*dt*dx + Hxx*dx^2
   * --> need 5D matrix
   */
  
  double f00 = (*f_of_t_x)[index_t][index_x];
  vector<double> f_reduced(5, 0.0);
  vector<double> coeffs(5, 0.0);
  /* f10_reduced: */ f_reduced[0] = (*f_of_t_x)[index_t+1][index_x]-f00;
  /* f20_reduced: */ f_reduced[1] = (*f_of_t_x)[index_t+2][index_x]-f00;
  /* f01_reduced: */ f_reduced[2] = (*f_of_t_x)[index_t][index_x+1]-f00;
  /* f11_reduced: */ f_reduced[3] = (*f_of_t_x)[index_t+1][index_x+1]-f00;
  /* f02_reduced: */ f_reduced[4] = (*f_of_t_x)[index_t][index_x+2]-f00;
  vector<vector<double> > A(5, vector<double>(5, 0.0));
  vector<vector<double> > A_inverse;
  A[0][0] = dt1;
  //A[0][1] = 0.0;
  A[0][2] = dt1*dt1;
  //A[0][3] = 0.0;
  //A[0][4] = 0.0;
  
  A[1][0] = dt2;
  //A[1][1] = 0.0;
  A[1][2] = dt2*dt2;
  //A[1][3] = 0.0;
  //A[1][4] = 0.0;
  
  //A[2][0] = 0.0;
  A[2][1] = dx1;
  //A[2][2] = 0.0;
  //A[2][3] = 0.0;
  A[2][4] = dx1*dx1;
  
  A[3][0] = dt1;
  A[3][1] = dx1;
  A[3][2] = dt1*dt1;
  A[3][3] = dt1*dx1;
  A[3][4] = dx1*dx1;
  
  //A[4][0] = 0.0;
  A[4][1] = dx2;
  //A[4][2] = 0.0;
  //A[4][3] = 0.0;
  A[4][4] = dx2*dx2;
  
  
  invert_matrix(&A, &A_inverse);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      coeffs[i] += A_inverse[i][j]*f_reduced[j];
    }
  }
  
  // coeffs[3] = 2.0*Hxt
  return coeffs[0] + 2.0*coeffs[2]*dt + coeffs[3]*dx;
  
}

double interpolate_grid_biquadratic_dx(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x){
  
  int n_t = (*t).size();
  int n_x = (*x).size();
  int index_t = find_index(t0, t);
  int index_x = find_index(x0, x);
  if(index_t >= n_t - 2) index_t = n_t - 3;
  if(index_x >= n_x - 2) index_x = n_x - 3;
  
  /*
   * (index_t, index_x)---> p00  p01  p02
   *                        p10  p11
   *                        p20
   * ISSUE: this treatment is not isotropic.
   * However, OF didn't want to go for bi-cubic interpolation,
   * because it is so slow.
   */
  
  double dt = t0 - (*t)[index_t];
  double dx = x0 - (*x)[index_x];
  double dt1 = (*t)[index_t+1] - (*t)[index_t];
  double dt2 = (*t)[index_t+2] - (*t)[index_t];
  double dx1 = (*x)[index_x+1] - (*x)[index_x];
  double dx2 = (*x)[index_x+2] - (*x)[index_x];
  
  /*
   * f - f00 = gt*dt + gx*dx + Htt*dt^2 + 2*Htx*dt*dx + Hxx*dx^2
   * --> need 5D matrix
   */
  
  double f00 = (*f_of_t_x)[index_t][index_x];
  vector<double> f_reduced(5, 0.0);
  vector<double> coeffs(5, 0.0);
  /* f10_reduced: */ f_reduced[0] = (*f_of_t_x)[index_t+1][index_x]-f00;
  /* f20_reduced: */ f_reduced[1] = (*f_of_t_x)[index_t+2][index_x]-f00;
  /* f01_reduced: */ f_reduced[2] = (*f_of_t_x)[index_t][index_x+1]-f00;
  /* f11_reduced: */ f_reduced[3] = (*f_of_t_x)[index_t+1][index_x+1]-f00;
  /* f02_reduced: */ f_reduced[4] = (*f_of_t_x)[index_t][index_x+2]-f00;
  vector<vector<double> > A(5, vector<double>(5, 0.0));
  vector<vector<double> > A_inverse;
  A[0][0] = dt1;
  //A[0][1] = 0.0;
  A[0][2] = dt1*dt1;
  //A[0][3] = 0.0;
  //A[0][4] = 0.0;
  
  A[1][0] = dt2;
  //A[1][1] = 0.0;
  A[1][2] = dt2*dt2;
  //A[1][3] = 0.0;
  //A[1][4] = 0.0;
  
  //A[2][0] = 0.0;
  A[2][1] = dx1;
  //A[2][2] = 0.0;
  //A[2][3] = 0.0;
  A[2][4] = dx1*dx1;
  
  A[3][0] = dt1;
  A[3][1] = dx1;
  A[3][2] = dt1*dt1;
  A[3][3] = dt1*dx1;
  A[3][4] = dx1*dx1;
  
  //A[4][0] = 0.0;
  A[4][1] = dx2;
  //A[4][2] = 0.0;
  //A[4][3] = 0.0;
  A[4][4] = dx2*dx2;
  
  
  invert_matrix(&A, &A_inverse);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      coeffs[i] += A_inverse[i][j]*f_reduced[j];
    }
  }
  
  // coeffs[3] = 2.0*Hxt
  return coeffs[1] + coeffs[3]*dt + 2.0*coeffs[4]*dx;
  
}




vector<double> interpolate_grid_biquadratic_derivs(double t0, double x0, vector<double> * t, vector<double> * x,
				       vector<vector<double> > * f_of_t_x){
  
  int n_t = (*t).size();
  int n_x = (*x).size();
  int index_t = find_index(t0, t);
  int index_x = find_index(x0, x);
  if(index_t >= n_t - 2) index_t = n_t - 3;
  if(index_x >= n_x - 2) index_x = n_x - 3;
  
  /*
   * (index_t, index_x)---> p00  p01  p02
   *                        p10  p11
   *                        p20
   * ISSUE: this treatment is not isotropic.
   * However, OF didn't want to go for bi-cubic interpolation,
   * because it is so slow.
   */
  
  double dt = t0 - (*t)[index_t];
  double dx = x0 - (*x)[index_x];
  double dt1 = (*t)[index_t+1] - (*t)[index_t];
  double dt2 = (*t)[index_t+2] - (*t)[index_t];
  double dx1 = (*x)[index_x+1] - (*x)[index_x];
  double dx2 = (*x)[index_x+2] - (*x)[index_x];
  
  /*
   * f - f00 = gt*dt + gx*dx + Htt*dt^2 + 2*Htx*dt*dx + Hxx*dx^2
   * --> need 5D matrix
   */
  
  double f00 = (*f_of_t_x)[index_t][index_x];
  vector<double> f_reduced(5, 0.0);
  vector<double> coeffs(5, 0.0);
  /* f10_reduced: */ f_reduced[0] = (*f_of_t_x)[index_t+1][index_x]-f00;
  /* f20_reduced: */ f_reduced[1] = (*f_of_t_x)[index_t+2][index_x]-f00;
  /* f01_reduced: */ f_reduced[2] = (*f_of_t_x)[index_t][index_x+1]-f00;
  /* f11_reduced: */ f_reduced[3] = (*f_of_t_x)[index_t+1][index_x+1]-f00;
  /* f02_reduced: */ f_reduced[4] = (*f_of_t_x)[index_t][index_x+2]-f00;
  vector<vector<double> > A(5, vector<double>(5, 0.0));
  vector<vector<double> > A_inverse;
  A[0][0] = dt1;
  //A[0][1] = 0.0;
  A[0][2] = dt1*dt1;
  //A[0][3] = 0.0;
  //A[0][4] = 0.0;
  
  A[1][0] = dt2;
  //A[1][1] = 0.0;
  A[1][2] = dt2*dt2;
  //A[1][3] = 0.0;
  //A[1][4] = 0.0;
  
  //A[2][0] = 0.0;
  A[2][1] = dx1;
  //A[2][2] = 0.0;
  //A[2][3] = 0.0;
  A[2][4] = dx1*dx1;
  
  A[3][0] = dt1;
  A[3][1] = dx1;
  A[3][2] = dt1*dt1;
  A[3][3] = dt1*dx1;
  A[3][4] = dx1*dx1;
  
  //A[4][0] = 0.0;
  A[4][1] = dx2;
  //A[4][2] = 0.0;
  //A[4][3] = 0.0;
  A[4][4] = dx2*dx2;
  
  
  invert_matrix(&A, &A_inverse);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      coeffs[i] += A_inverse[i][j]*f_reduced[j];
    }
  }
  
  // coeffs[3] = 2.0*Hxt
  vector<double> output(3,0.0);
  output[0] = f00 + coeffs[0]*dt + coeffs[1]*dx + coeffs[2]*dt*dt + coeffs[3]*dt*dx + coeffs[4]*dx*dx;
  output[1] = coeffs[0] + 2.0*coeffs[2]*dt + coeffs[3]*dx;
  output[2] = coeffs[1] + coeffs[3]*dt + 2.0*coeffs[4]*dx;
  return output;
  
}
