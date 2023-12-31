#include <Rcpp.h>

using namespace Rcpp;

void ks_matrix_R(Rcpp::NumericVector t_X, Rcpp::NumericVector t_R,
            Rcpp::IntegerVector t_sidxs, Rcpp::IntegerVector t_n_genes,
            Rcpp::IntegerVector t_geneset_idxs, Rcpp::IntegerVector t_n_geneset,
            Rcpp::NumericVector t_tau, Rcpp::IntegerVector t_n_samples,
            Rcpp::IntegerVector t_mx_diff, Rcpp::IntegerVector t_abs_rnk);

double
ks_sample(double* x, int* x_sort_indxs, int n_genes, int* geneset_mask,
          int* geneset_idxs, int n_geneset, double tau, int mx_diff, int abs_rnk){


	double dec = 1.0 / (n_genes - n_geneset);
	double sum_gset = 0.0;
	for(int i = 0; i < n_geneset; ++i){
		sum_gset += pow(x[geneset_idxs[i]-1], tau);
	}

	//double mx_value = 0.0;
	double mx_value_sign = 0.0;
	double cum_sum = 0.0;

	double mx_pos = 0.0;
	double mx_neg = 0.0;

	int idx;
	for(int i = 0; i < n_genes; ++i){
		idx = x_sort_indxs[i]-1;

		if(geneset_mask[idx] == 1){
			cum_sum += pow(x[idx], tau) / sum_gset;
		}else{
			cum_sum -= dec;
		}

		if(cum_sum > mx_pos){ mx_pos = cum_sum; }
		if(cum_sum < mx_neg){ mx_neg = cum_sum; }
	}

	if (mx_diff != 0) {
		mx_value_sign = mx_pos + mx_neg;
    if (abs_rnk != 0)
      mx_value_sign = mx_pos - mx_neg;
	} else {
		mx_value_sign = (mx_pos > fabs(mx_neg)) ? mx_pos : mx_neg;
	}
	return mx_value_sign;
}


/**
 * X <- gene density scores
 * R <- result
 * sidxs <- sorted gene densities idxs
 */
void ks_matrix(double* X, double* R, int* sidxs, int n_genes, int* geneset_idxs,
               int n_geneset, double tau, int n_samples, int mx_diff, int abs_rnk){
	int geneset_mask[n_genes];
	for(int i = 0; i < n_genes; ++i){
		geneset_mask[i] = 0;
	}

	for(int i = 0; i < n_geneset; ++i){
		geneset_mask[geneset_idxs[i]-1] = 1;
	}

	for(int j = 0; j < n_samples; ++j){
		int offset = j * n_genes;
    R[j] = ks_sample(&X[offset], &sidxs[offset], n_genes, &geneset_mask[0],
                     geneset_idxs, n_geneset, tau, mx_diff, abs_rnk);
	}
}

// [[Rcpp::export]]
void ks_matrix_R(Rcpp::NumericVector t_X, Rcpp::NumericVector t_R,
                 Rcpp::IntegerVector t_sidxs, Rcpp::IntegerVector t_n_genes,
                 Rcpp::IntegerVector t_geneset_idxs, Rcpp::IntegerVector t_n_geneset,
                 Rcpp::NumericVector t_tau, Rcpp::IntegerVector t_n_samples,
                 Rcpp::IntegerVector t_mx_diff, Rcpp::IntegerVector t_abs_rnk) {
  double* X = t_X.begin();
  double* R = t_R.begin();
  int* sidxs = t_sidxs.begin();
  int* n_genes = t_n_genes.begin();
  int* geneset_idxs = t_geneset_idxs.begin();
  int* n_geneset = t_n_geneset.begin();
  double* tau = t_tau.begin();
  int* n_samples = t_n_samples.begin();
  int* mx_diff = t_mx_diff.begin();
  int* abs_rnk = t_abs_rnk.begin();

	ks_matrix(X, R, sidxs, *n_genes, geneset_idxs, *n_geneset, *tau, *n_samples, *mx_diff, *abs_rnk);
}
