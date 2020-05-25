#ifndef _GalaxySample
#define _GalaxySample

#include <gsl/gsl_sf_gamma.h>

class GalaxySample {

 public:

   GalaxySample(FlatInhomogeneousUniverseLCDM* universe, double b1, double b2, double a0, double a1);
   ~GalaxySample();
   
   // setting class attributes
   void set_universe(FlatInhomogeneousUniverseLCDM* universe);
   void set_parameters(double b1, double b2, double a0, double a1);
   void set_bias_model_from_br_parametrisation(double b_tilde, double r, double N_bar, double variance, double skewness);
   void set_error_flag_negative_density(int flag){this->error_flag_negative_density = flag;};
   double set_b2_to_minimise_negative_densities(double variance); // also returns new b2
   
   // returning (pointers to) class attributes
   int get_error_flag_negative_density(){return this->error_flag_negative_density;};
   FlatInhomogeneousUniverseLCDM* pointer_to_universe(){return this->universe;};
   
   // class output
   int return_N_max(double N_bar, double variance);
   double return_P_of_N_given_delta(int N, double N_bar, double delta, double variance); // variance only used when quadratic_bias != 0
   vector<double> return_CIC_from_matter_density_PDF(double N_bar, vector<vector<double> > PDF_data);
   
   
 private:
   
   double linear_bias;
   double quadratic_bias;
   double alpha_0; // shot-noise parameters a la https://arxiv.org/abs/1710.05162
   double alpha_1; // shot-noise parameters a la https://arxiv.org/abs/1710.05162
   // --> ISSUE: this pretends that quadratic_bias, alpha_0 and alpha_1 can be the same on all scales. That is 
   //            however not possible. Scale dependence should be made explicit in the future!
   
   FlatInhomogeneousUniverseLCDM* universe;
   cosmological_model cosmology;
   
   int error_flag_negative_density = 0;
   
};

#include "GalaxySample.cpp"
#endif
