#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
NumericVector h (
    NumericVector& y,
    NumericVector& cvals
){
  int   k = y.size();
  NumericVector ans(k);
  ans   = ans + cvals(0);
  for (int j=0; j<4; j++) {
    ans += pow(y,j+1)*cvals(j+1);
  }
  return ans;
} // END h



// [[Rcpp::export]]
NumericVector hp (
    NumericVector&  y,
    double&         pt
){
  int   k = y.size();
  NumericVector ans(k);
  for (int j=0; j<4; j++) {
    ans += pow(-y/pt, j);
    ans = ans * pow(-pt, -1);
  }
  return ans;
} // END hp



// [[Rcpp::export]]
NumericVector hpp (
    NumericVector&  y,
    double&         pt
){
  int   k = y.size();
  NumericVector ans(k);
  for (int j=0; j<3; j++) {
    ans += (j + 1) * pow(-y/pt, j);
    ans = ans * pow(-pt, -2);
  }
  return ans;
} // END hpp  



// [[Rcpp::export]]
mat mllog (
  vec&          x,        // kx1 vector
  double&       eps,      // default value 1/k
  double        M,        // default value Inf
  const int     der       // default value 0
) {
  int   k = x.n_elem;
  NumericVector   xx = wrap(x);
  
  if (eps>M) stop("Thresholds out of order");
  
  LogicalVector   lo = xx < eps;
  LogicalVector   hi = xx > M;
  LogicalVector   md = (!lo) & (!hi);
  
  // Coefficients for 4th order Taylor approx below eps
  NumericVector coefs(5);
  coefs(0)      = -log(eps);
  for (int i=1; i<5; i++) {
    coefs(i)    = pow(-eps, -i)/i;
  }
  
  // Coefficients for 4th order Taylor approx above M
  NumericVector Coefs(5);
  Coefs(0)      = -log(M);
  for (int i=1; i<5; i++) {
    Coefs(i)    = pow(-M, -i)/i;
  }
  
  // function value
  NumericVector f(k);
  
  NumericVector f_tmp_lo = xx[lo] ;
  f_tmp_lo      = f_tmp_lo - eps;
  f[lo]         = h(f_tmp_lo, coefs);
  
  NumericVector f_tmp_hi = xx[hi] ;
  f_tmp_hi      = f_tmp_hi - M;
  f[hi]         = h(f_tmp_hi, Coefs);
  
  NumericVector xxmd      = xx[md];
  NumericVector f_tmp_md  = -log(xxmd);
  f[md]         = f_tmp_md;
  
  if (der<1) return as<colvec>(f);
  
  // first derivative
  NumericVector fp(k);
  fp[lo]        = hp(f_tmp_lo, eps);
  fp[hi]        = hp(f_tmp_hi, M);
  
  NumericVector fp_tmp_md  = -pow(xxmd, -1);
  fp[md]        = fp_tmp_md;
  
  if (der <2) return join_rows(as<colvec>(f), as<colvec>(fp));
  
  // second derivative
  NumericVector fpp(k);
  fpp[lo]       = hpp(f_tmp_lo, eps);
  fpp[hi]       = hpp(f_tmp_hi, M);
  NumericVector fpp_tmp_md  = pow(xxmd, -2);
  fpp[md]       = fpp_tmp_md;
  
  return join_rows(as<colvec>(f), as<colvec>(fp), as<colvec>(fpp));
} // END mllog



// [[Rcpp::export]]
vec svdlm (
  mat& X,
  mat& y
) {
  // Tolerances for generalized inverse via SVD
  const double RELTOL = 1e-9;
  const double ABSTOL = 1e-100;
  
  mat U;
  vec d;
  mat V;
  svd_econ(U,d,V,X);

  int D = d.n_elem;
  vec dinv(D);
  const double  threshold = RELTOL * max(d) + ABSTOL;
  for (int i=0; i<D; i++) {
    if (d(i) < threshold) {
      dinv(i) = 0;
    } else {
      dinv(i) = pow(d(i), -1);
    }
  }
  mat Xplus   = V * diagmat(dinv) * U.head_cols(D).t();
  
  return Xplus * y;
} // END svdlm

  
  
// [[Rcpp::export]]
List emplik (
  mat     z,
  Nullable<NumericVector> mu      = R_NilValue,      // These values are filled by default with the same values as in scel.R
  Nullable<NumericVector> lam     = R_NilValue,
  Nullable<double>        eps     = R_NilValue,
  Nullable<double>        M       = R_NilValue,
  Nullable<double>        thresh  = R_NilValue,
  Nullable<int>           itermax = R_NilValue,
  Nullable<bool>          verbose = R_NilValue
) {
  const int n = z.n_rows;
  const int d = z.n_cols;
  double    nnn = n;

  // Backtracking line search parameters [Tweak only with extreme caution.]
  // See Boyd and Vandenberghe, pp464-466.
  const double ALPHA = 0.3;               // seems better than 0.01 on some 2d test data (sometimes fewer iters)
  const double BETA  = 0.8;
  // We need  0 < ALPHA < 1/2   and 0 < BETA < 1

  // Backtrack threshold: you can miss by this much.
  const double BACKEPS = 0;
  // Consider replacing 0 by 1e-10 if backtracking seems to be
  // failing due to round off.

  vec MU;
  if (mu.isNotNull()) {
    MU = as<vec>(mu);
  } else {
    MU = zeros<vec>(d);
  }
  z.each_row() -= MU.t();

  double EPS;
  if (eps.isNotNull()) {
    EPS = as<double>(eps);
  } else {
    EPS = 1/nnn;
  }

  double MM;
  if (M.isNotNull()) {
    MM = as<double>(M);
  } else {
    const double MM_tmp = datum::inf;
    MM = MM_tmp;
  }

  double THRESH;
  if (thresh.isNotNull()) {
    THRESH = as<double>(thresh);
  } else {
    THRESH = 1e-30;
  }

  int ITERMAX;
  if (itermax.isNotNull()) {
    ITERMAX = as<int>(itermax);
  } else {
    ITERMAX = 100;
  }

  bool VERBOSE;
  if (verbose.isNotNull()) {
    VERBOSE = as<bool>(verbose);
  } else {
    VERBOSE = false;
  }

  // Use lam = 0 or initial lam, whichever is best
  vec iota    = ones<vec>(n);
  mat init0   = mllog(iota, EPS, MM, 2);            // i.e. lam = 0

  vec LAM;
  mat init;
  if (lam.isNotNull()) {
    LAM     = as<vec>(lam);
    vec init_tmp  = iota + z*LAM;
    init    = mllog(init_tmp, EPS, MM, 2);
    if ( accu(init0.head_cols(1)) < accu(init.head_cols(1)) ) {
      LAM     = zeros<vec>(d);
      init    = init0;
    }
  } else {
    LAM     = zeros<vec>(d);
    init    = init0;
  }

  // Initial f, g
  double fold = accu(init.head_cols(1));
  mat ZZ_INIT = z;
  for (int i=0; i<d; i++) {
    ZZ_INIT.col(i) = ZZ_INIT.col(i) % init.col(1);
  }
  rowvec gold = sum(ZZ_INIT, 0);

  bool  converged = false;
  int   iter      = 0;
  mat   oldvals   = init;

  double ndec;
  double gradnorm;

  while ( !converged ) {
    iter++;
  
    // Get Newton Step
    vec rootllpp  = pow(oldvals.col(2), 0.5);         // sqrt 2nd deriv of -llog lik
    mat zt        = z;
    for (int i=0; i<d; i++) zt.col(i) = zt.col(i) % rootllpp;
    vec yt        = oldvals.col(1) % pow(rootllpp, -1);
    vec step      = -svdlm(zt, yt);                   // more reliable than step = -lm( yt~zt-1 )$coef
    
    bool backtrack = false;
    int s = 1;                                        // usually called t
    while( !backtrack ) {
      vec nv_tmp  = iota + z * (LAM + s * step);
      mat newvals = mllog(nv_tmp, EPS, MM, 2);
      double fnew = accu(newvals.col(0));
      double targ = fold + ALPHA * s * accu(gold.t() % step) + BACKEPS;      //  (BACKEPS for roundoff, should not be needed)
      bool   fff  = fnew <= targ;
      
      if ( fff ) {
        // backtracking has converged
        backtrack = true;
        oldvals   = newvals;
        fold      = fnew;
        ZZ_INIT   = z;
        for (int i=0; i<d; i++) {
          ZZ_INIT.col(i) = ZZ_INIT.col(i) % oldvals.col(1);
        }
        gold      = sum(ZZ_INIT, 0);
        // take the step
        LAM      += s * step;
      } else {
        s         = s * BETA;
      }
    } // END while( !backtrack )
    
    // Newton decrement and gradient norm
    ndec     = pow(accu(pow(gold.t() % step, 2)), 0.5);
    gradnorm = pow( accu(pow(gold, 2)), 0.5);

    NumericVector lll = wrap(LAM);
    // if ( verbose ) Rcout << NumericVector::create(fold, gradnorm, ndec) << " " << lll << endl;
    converged       = pow(ndec, 2) <= THRESH;
    if ( iter>= ITERMAX ) break;
  } // END while (!converged)
  
  vec     wts_tmp   = iota + z * LAM;
  vec     wts       = (1/nnn) * pow(wts_tmp, -1);
  double  logelr    = accu( mllog(wts_tmp, EPS, MM, 0) );

  return List::create(
    _["logelr"]     = logelr,
    _["lam"]        = LAM,
    _["wts"]        = wts,
    _["converged"]  = converged,
    _["iter"]       = iter,
    _["ndec"]       = ndec,
    _["gradnorm"]   = gradnorm
  );
} // END emplik

