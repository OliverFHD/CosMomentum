#ifndef _GalaxySample
#define _GalaxySample

#include <gsl/gsl_sf_gamma.h>

class GalaxySample {

 public:

   GalaxySample(FlatInhomogeneousUniverseLCDM* universe, double b1, double b2, double a0, double a1, BIAS_MODEL b_model);
   ~GalaxySample();
   
   // setting class attributes
   void set_universe(FlatInhomogeneousUniverseLCDM* universe);
   void set_parameters(double b1, double b2, double a0, double a1, BIAS_MODEL b_model);
   void set_error_flag_negative_density(int flag){this->error_flag_negative_density = flag;};
   double set_b2_to_minimise_negative_densities(double variance); // also returns new b2
   
   // returning (pointers to) class attributes
   int get_error_flag_negative_density(){return this->error_flag_negative_density;};
   
   // class output
   BIAS_MODEL return_bias_model();
   int return_N_max(double N_bar, double variance);
   double delta_g_Eulerian(double delta, double variance);
   vector<vector<double> > delta_g_Lagrangian(vector<vector<double> > *PDF_data);
   double return_P_of_N_given_delta_g(int N, double N_bar, double delta_m, double delta_g);
   vector<double> return_CIC_from_matter_density_PDF(double N_bar, vector<vector<double> > PDF_data);
   
   FlatInhomogeneousUniverseLCDM* pointer_to_universe(){return this->universe;};
  
   
 private:
   
   BIAS_MODEL bias_model = EULERIAN;
   double linear_bias;
   double quadratic_bias;
   double linear_Lagrangian_bias;
   double quadratic_Lagrangian_bias;
   double alpha_0; // shot-noise parameters a la https://arxiv.org/abs/1710.05162
   double alpha_1; // shot-noise parameters a la https://arxiv.org/abs/1710.05162
   // --> ISSUE: this pretends that quadratic_bias, alpha_0 and alpha_1 can be the same on all scales. That is 
   //            however not possible. Scale dependence should be made explicit in the future!
   
   // These are only needed if the Lagrangian bias model is used:
   // <delta_g | delta_m > = bias_term_1 + b_1_Lagrange*bias_term_2 + b_2_Lagrange*bias_term_3
   vector<double> bias_term_1;
   vector<double> bias_term_2;
   vector<double> bias_term_3;
   
   FlatInhomogeneousUniverseLCDM* universe;
   cosmological_model cosmology;
   
   int error_flag_negative_density = 0;
   
};

#include "GalaxySample.cpp"
#endif
