#ifndef _ProjectedGalaxySample
#define _ProjectedGalaxySample
#include "GalaxySample.h"


class ProjectedGalaxySample : public GalaxySample {

 public:

   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, BIAS_MODEL b_model, string input_file);
   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, double a_IA, double alph_IA, BIAS_MODEL b_model, string input_file);
   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, double a_IA, double alph_IA, double a_lens, BIAS_MODEL b_model, string input_file);
   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, string input_file, string bias_and_shotnoise);
   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double a_IA, double alph_IA, string input_file, string bias_and_shotnoise);
   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double a_IA, double alph_IA, double a_lens, string input_file, string bias_and_shotnoise);
   ~ProjectedGalaxySample();
   
   void set_n_of_w_data(string n_of_z_input_file);
   void set_parameters_projected(double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, BIAS_MODEL b_model, string input_file);
   double set_b2_to_minimise_negative_densities_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   double compute_variance_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   double compute_skewness_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   
   int return_N_max_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   vector<double> return_CiC_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<double> return_CiC_PDF_in_angular_tophat_Gauss(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<double> return_CiC_saddle_point_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, int next_to_leading);
   vector<vector<double> > return_kappa_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<vector<double> > return_matter_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<vector<double> > return_matter_saddle_point_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, int next_to_leading);
   void return_joint_saddle_point_PDF_Ng_kappaCMB_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
   void return_joint_saddle_point_PDF_Ng_kappaCMB_noisy_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_CMB_noise_variance, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
   void return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, vector<double> w_values_lensing_kernel, vector<double> lensing_kernel, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
   void return_LOS_data(vector<double> *w_vals, vector<double> *n_of_w_vals, vector<double> *lensing_kernel_vals);
   
 private:
   
   vector<double> w_values;
   vector<double> n_of_w_values;
   vector<double> lensing_kernel_values;
   
   int grid_initialised = 0;
   double variance = 0;
   vector<vector<double> > d_grid;
   vector<vector<double> > k_grid;
   vector<vector<double> > p_grid;
   
   double density; // in radians
   double A_lens   = 1.0;
   double A_IA     = 0.0;
   double alpha_IA = 0.0;
   double z0       = 0.62;
   
};


#include "ProjectedGalaxySample.cpp"
#endif
