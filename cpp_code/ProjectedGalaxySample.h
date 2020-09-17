#ifndef _ProjectedGalaxySample
#define _ProjectedGalaxySample
#include "GalaxySample.h"


class ProjectedGalaxySample : public GalaxySample {

 public:

   ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, string input_file);
   ~ProjectedGalaxySample();
   
   void set_n_of_w_data(string n_of_z_input_file);
   void set_parameters_projected(double density_in_arcmin_squared, double b1, double b2, double a0, double a1, string input_file);
   double set_b2_to_minimise_negative_densities_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   void set_projected_bias_model_from_br_parametrisation(double b_tilde, double r, double theta_in_arcmin, double f_NL, double var_NL_rescale);
   double compute_variance_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   double compute_skewness_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   
   int return_N_max_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   vector<double> return_CiC_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<double> return_CiC_saddle_point_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<vector<double> > return_kappa_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<vector<double> > return_matter_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   vector<vector<double> > return_matter_saddle_point_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   void return_joint_saddle_point_PDF_Ng_kappaCMB_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
   void return_joint_saddle_point_PDF_Ng_kappaCMB_noisy_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_CMB_noise_variance, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
   void return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, vector<double> w_values_lensing_kernel, vector<double> lensing_kernel, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
   void return_LOS_data(vector<double> *w_vals, vector<double> *n_of_w_vals, vector<double> *lensing_kernel_vals);
   
 private:
   
   vector<double> w_values;
   vector<double> n_of_w_values;
   vector<double> lensing_kernel_values;
   
   double density; // in radians
   
};


#include "ProjectedGalaxySample.cpp"
#endif
