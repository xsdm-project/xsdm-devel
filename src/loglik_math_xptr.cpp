// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(ucminfcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <memory>
#include <vector>
#include "ucminf_core.hpp"

namespace {

double eval_loglik_math(
    const Rcpp::Function& loglik_math_fn,
    const Rcpp::NumericVector& par,
    const Rcpp::RObject& env_dat,
    const Rcpp::RObject& occ,
    const Rcpp::RObject& mask,
    const int num_threads
) {
  return Rcpp::as<double>(loglik_math_fn(
    Rcpp::_["param_vector"] = par,
    Rcpp::_["env_dat"] = env_dat,
    Rcpp::_["occ"] = occ,
    Rcpp::_["mask"] = mask,
    Rcpp::_["num_threads"] = num_threads,
    Rcpp::_["negative"] = true
  ));
}

} // namespace

// [[Rcpp::export]]
SEXP make_loglik_math_xptr(
    SEXP env_dat,
    SEXP occ,
    SEXP mask = R_NilValue,
    int num_threads = 1,
    std::string grad = "central",
    Rcpp::NumericVector gradstep = Rcpp::NumericVector::create(1e-6, 1e-8)
) {
  if (gradstep.size() != 2) {
    Rcpp::stop("`gradstep` must have length 2.");
  }
  if (grad != "forward" && grad != "central") {
    Rcpp::stop("`grad` must be either 'forward' or 'central'.");
  }

  const double gradstep_rel = gradstep[0];
  const double gradstep_abs = gradstep[1];
  const bool use_central = (grad == "central");

  Rcpp::Environment ns = Rcpp::Environment::namespace_env("xsdm");
  Rcpp::Function loglik_math_fn = ns["loglik_math"];
  Rcpp::RObject env_dat_obj(env_dat);
  Rcpp::RObject occ_obj(occ);
  Rcpp::RObject mask_obj(mask);

  auto fn = std::make_unique<ucminf::ObjFun>(
    [loglik_math_fn, env_dat_obj, occ_obj, mask_obj, num_threads,
     gradstep_rel, gradstep_abs, use_central]
    (const std::vector<double>& x, std::vector<double>& g, double& f) {
      const int n = static_cast<int>(x.size());
      Rcpp::NumericVector x_r(n);
      for (int i = 0; i < n; ++i) {
        x_r[i] = x[i];
      }

      f = eval_loglik_math(
        loglik_math_fn, x_r, env_dat_obj, occ_obj, mask_obj, num_threads
      );

      for (int i = 0; i < n; ++i) {
        const double xi = x_r[i];
        const double dx = std::abs(xi) * gradstep_rel + gradstep_abs;

        x_r[i] = xi + dx;
        const double f_plus = eval_loglik_math(
          loglik_math_fn, x_r, env_dat_obj, occ_obj, mask_obj, num_threads
        );

        if (use_central) {
          x_r[i] = xi - dx;
          const double f_minus = eval_loglik_math(
            loglik_math_fn, x_r, env_dat_obj, occ_obj, mask_obj, num_threads
          );
          g[i] = (f_plus - f_minus) / (2.0 * dx);
        } else {
          g[i] = (f_plus - f) / dx;
        }

        x_r[i] = xi;
      }
    }
  );

  return Rcpp::XPtr<ucminf::ObjFun>(fn.release(), true);
}
