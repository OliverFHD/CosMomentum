#include <math.h>

using namespace std;
using namespace constants;

double z_equality(double Omega_0, double h, double theta_27){
  return 2.5*Omega_0*h*h*pow(10.0/theta_27, 4.0);	
}

double z_drag(double Omega_0, double Omega_b, double h, double theta_27){
  double Om_times_h_sq = Omega_0*h*h;
  double Ob_times_h_sq = Omega_b*h*h;
  double b1 = 0.313*pow(Om_times_h_sq, -0.419)*(1.0+0.607*pow(Om_times_h_sq, 0.674));
  double b2 = 0.238*pow(Om_times_h_sq, 0.223);	
	
  return 1291.0*pow(Om_times_h_sq, 0.251)/(1.0+0.659*pow(Om_times_h_sq, 0.828))*(1.0+b1*pow(Ob_times_h_sq, b2));	
}

double k_equality(double Omega_0, double h, double theta_27){
  double z_eq = z_equality(Omega_0, h, theta_27);
  //return 7.46e-2*Omega_0*h*h/theta_27/theta_27/h*c_over_e5;	
  return sqrt(2.0*Omega_0*z_eq);	
}

double k_equality(double Omega_0, double h, double theta_27, double z_eq){
  //return 7.46e-2*Omega_0*h*h/theta_27/theta_27/h*c_over_e5;
  return sqrt(2.0*Omega_0*z_eq);	
}

double k_silk(double Omega_0, double Omega_b, double h, double theta_27){
  return 1.6*pow(Omega_b*h*h, 0.52)*pow(Omega_0*h*h, 0.73)*(1.0+pow(10.4*Omega_0*h*h, -0.95))*c_over_e5/h;
}

double baryon_photon_momentum_ratio(double Omega_b, double h, double theta_27, double z){
  //return 0.75*Omega_b/(2.469e-5/h/h)/(1+z);
  return 31.5*Omega_b*h*h/pow(theta_27, 4.0)*1000.0/z;
}

// alpha_c
double CDM_suppression(double Omega_0, double Omega_b, double h){
  double Om_times_h_sq = Omega_0*h*h;
  double Omega_ratio = Omega_b/Omega_0;
  double a1 = pow(46.9*Om_times_h_sq, 0.670)*(1.0+pow(32.1*Om_times_h_sq, -0.532));
  double a2 = pow(12.0*Om_times_h_sq, 0.424)*(1.0+pow(45.0*Om_times_h_sq, -0.582));
  return 1.0/(pow(a1, Omega_ratio)*pow(a2,pow(Omega_ratio, 3.0)));
}

// beta_c
double CDM_shift(double Omega_0, double Omega_c, double h){
  double Om_times_h_sq = Omega_0*h*h;
  double Omega_ratio = Omega_c/Omega_0;
  double b1 = 0.944/(1.0+pow(458.0*Om_times_h_sq, -0.708));
  double b2 = pow(0.395*Om_times_h_sq, -0.0266);
  return 1.0/(1.0+b1*(pow(Omega_ratio, b2)-1.0));
}

double sound_horizon(double Omega_0, double Omega_b, double h, double theta_27){
  double z_eq = z_equality(Omega_0, h, theta_27);
  double z_d = z_drag(Omega_0, Omega_b, h, theta_27);
  double k_eq = k_equality(Omega_0, h, theta_27, z_eq);
  double R_eq = baryon_photon_momentum_ratio(Omega_b, h, theta_27, z_eq);
  double R_d = baryon_photon_momentum_ratio(Omega_b, h, theta_27, z_d);
  
  return 2.0/(3.0*k_eq)*sqrt(6.0/R_eq)*log((sqrt(1.0+R_d)+sqrt(R_d+R_eq))/(1.0+sqrt(R_eq)));
}

double f_factor(double k, double s){
  return 1.0/(1.0+pow(k*s/5.4,4.0));
}

double T_tilde_0(double k, double alpha, double beta, double k_eq){
  double q = k/(13.41*k_eq);
  //q = k*1.009259259*1.009259259/(0.3*0.7)/c_over_e5;
  double C = 14.2/alpha+386.0/(1.0+69.9*pow(q, 1.08));
  double log_term = log(eulers_constant+1.8*beta*q);
  return log_term/(log_term+C*q*q);
}
double transfer_c(double k, double s, double alpha_c, double beta_c, double k_eq){
  double f = f_factor(k, s);
  double T0_1 = T_tilde_0(k, 1.0, beta_c, k_eq);
  double T0_2 = T_tilde_0(k, alpha_c, beta_c, k_eq);
  return f*T0_1+(1.0-f)*T0_2;
}

double s_tilde(double k, double s, double Om_times_h_sq){
  double beta_node = 8.41*pow(Om_times_h_sq, 0.435);
  return s/(pow(1.0+pow(beta_node/(k*s), 3.0),1.0/3.0));
}
double G_function(double y){
  double sqrt_one_plus_y = sqrt(1.0+y);
  return y*(-6.0*sqrt_one_plus_y+(2.0+3.0*y)*log((sqrt_one_plus_y+1.0)/(sqrt_one_plus_y-1.0)));
}
double baryon_suppression(double k_eq, double s, double z_eq, double z_d, double R_d){
  double ratio = (1.0+z_eq)/(1.0+z_d);
  return 2.07*k_eq*s*pow(1.0+R_d, -0.75)*G_function(ratio);
}
double baryon_shift(double Omega_0, double Omega_b, double h){
  double ratio = Omega_b/Omega_0;
  return 0.5 + ratio + (3.0-2.0*ratio)*sqrt(pow(17.2*Omega_0*h*h, 2.0)+1.0);
}
double transfer_b(double k, double s, double alpha_b, double beta_b, double k_eq, double k_s, double Omega_0, double h){
  double T0 = T_tilde_0(k, 1.0, 1.0, k_eq);
  double ss = s_tilde(k, s, Omega_0*h*h);
  double j0 = sin(k*ss)/(k*ss);
  return j0*(T0/(1.0+pow(k*s/5.2,2.0)) + alpha_b/(1.0+pow(beta_b/(k*s),3.0))*exp(-pow(k/k_s,1.4)));
}

void tranfer_function_Eisenstein_and_Hu(vector<double> *wave_numbers, vector<double> *transfer_function, double Omega_0, double Omega_b, double h, double theta_27){
  
  double Omega_c = Omega_0-Omega_b;/**/
  double ratio_b0 = Omega_b/Omega_0;/**/
  double ratio_c0 = Omega_c/Omega_0;/**/
  double z_eq = z_equality(Omega_0, h, theta_27);/**/
  double z_d = z_drag(Omega_0, Omega_b, h, theta_27);/**/
  double R_d = baryon_photon_momentum_ratio(Omega_b, h, theta_27, z_d);/*naja*/
  double k_eq = k_equality(Omega_0, h, theta_27, z_eq);/**/
  double k_s = k_silk(Omega_0, Omega_b, h, theta_27);/**/
  double s = sound_horizon(Omega_0, Omega_b, h, theta_27);/**/
  double alpha_c = CDM_suppression(Omega_0, Omega_b, h);/**/
  double alpha_b = baryon_suppression(k_eq, s, z_eq, z_d, R_d);/**/
  double beta_b = baryon_shift(Omega_0, Omega_b, h);/**/
  double beta_c = CDM_shift(Omega_0, Omega_c, h);/**/
  double k, Tc, Tb;
  int n = (*wave_numbers).size();
  
  transfer_function->resize(n);
  for(int i = 0; i < n; i++){
    k = (*wave_numbers)[i];
    Tb = transfer_b(k, s, alpha_b, beta_b, k_eq, k_s, Omega_0, h);
    Tc = transfer_c(k, s, alpha_c, beta_c, k_eq);
    (*transfer_function)[i] = ratio_b0*Tb+ratio_c0*Tc;
  }
}




// NOTE you have made significant changes to this in March 2019, due to a change of units
void tranfer_function_Bond_and_Efstathiou(vector<double> *wave_numbers, vector<double> *transfer_function, double Omega_0, double h){
  

  int n = (*wave_numbers).size();
	double k;
	double q;
	double a = 6.4*h/c_over_e5;
	double b = 3.0*h/c_over_e5;
	double c = 1.7*h/c_over_e5;
	double nu = 1.13;
	double Gamma = Omega_0;
  
  transfer_function->resize(n);
  for(int i = 0; i < n; i++){
    k = (*wave_numbers)[i];
		q = k/Gamma;
    (*transfer_function)[i] = 1.0/pow(1 + pow(a*q + pow(b*q, 1.5) + pow(c*q, 2), nu), 1.0/nu);
  }
}

