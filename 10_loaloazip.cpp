#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_VECTOR(N);
  // Params
  PARAMETER_VECTOR(eta1);
  PARAMETER_VECTOR(eta2);

  // Constants
  int n = y.size();

  // Transformations
  // Prob of nonstructural zero
  vector<Type> zeroprob(n);
  // Prob of infection
  vector<Type> infectionprob(n);
  for (int i = 0;i<n;i++) {
    zeroprob(i) = 1.0 / (1.0 + exp(-1.0 * eta1(i)));
    infectionprob(i) = 1.0 / (1.0 + exp(-1.0 * eta2(i)));
  }

  // Log likelihood
  Type ll = 0;
  Type L = 0;
  for (int i = 0;i < n;i++) {
    // Zero prob
    if (y(i) == 0) {
      L += 1.0 - zeroprob(i);
    }
    // Binomial variability
    L += zeroprob(i) * exp(lfactorial(N(i)) - lfactorial(y(i)) - lfactorial(N(i) - y(i))) *
      pow(infectionprob(i),y(i)) *
      pow(1.0 - infectionprob(i),N(i) - y(i));

    ll += log(L);
    L = 0;
  }
  return ll;
}
