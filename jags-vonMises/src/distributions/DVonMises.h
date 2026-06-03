#ifndef DVONMISES_H_
#define DVONMISES_H_

#include <distribution/ScalarDist.h>

namespace jags {
namespace vonmises {

/**
 * von Mises distribution on [0, 2*pi).
 *
 * Parameters:
 *   mu        - mean direction in [0, 2*pi)
 *   precision - concentration kappa >= 0
 *
 * Density: exp(kappa * cos(x - mu)) / (2*pi * I_0(kappa))
 *
 * References:
 *   Forbes et al. (2011). Statistical Distributions, 4th ed.
 *   Best & Fisher (1979). J. R. Statist. Soc. Ser. C 28, 152-157.
 */
class DVonMises : public ScalarDist {
 public:
  DVonMises();

  double logDensity(double x, PDFType type,
                    std::vector<double const *> const &parameters,
                    double const *lower, double const *upper) const;
  double randomSample(std::vector<double const *> const &parameters,
                      double const *lower, double const *upper,
                      RNG *rng) const;
  double typicalValue(std::vector<double const *> const &parameters,
                      double const *lower, double const *upper) const;
  bool checkParameterValue(std::vector<double const *> const &parameters) const;

  double l(std::vector<double const *> const &parameters) const;
  double u(std::vector<double const *> const &parameters) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;

 private:
  static double twoPi();
  static double logTwoPi();
};

}  // namespace vonmises
}  // namespace jags

#endif  /* DVONMISES_H_ */
