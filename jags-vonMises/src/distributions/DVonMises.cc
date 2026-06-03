#include "DVonMises.h"

#include <cmath>
#include <vector>

#include <JRmath.h>
#include <rng/RNG.h>
#include <util/nainf.h>

namespace jags {
namespace vonmises {

namespace {

const double kTwoPi = 6.283185307179586;
const double kLogTwoPi = 1.837877066409345;

double wrapAngle(double x) {
  x = std::fmod(x, kTwoPi);
  if (x < 0.0) {
    x += kTwoPi;
  }
  return x;
}

}  // namespace

DVonMises::DVonMises()
    : ScalarDist("dvonmises", 2, DIST_SPECIAL) {}

double DVonMises::twoPi() { return kTwoPi; }

double DVonMises::logTwoPi() { return kLogTwoPi; }

double DVonMises::logDensity(double x, PDFType type,
                             std::vector<double const *> const &par,
                             double const * /*lower*/,
                             double const * /*upper*/) const {
  const double xw = wrapAngle(x);
  const double mu = wrapAngle(*par[0]);
  const double kappa = *par[1];

  if (kappa < 0.0) {
    return JAGS_NEGINF;
  }

  const double kernel = kappa * std::cos(xw - mu);

  if (type == PDF_PRIOR) {
    return kernel;
  }

  return kernel - logTwoPi() - std::log(bessel_i(kappa, 0.0, 1.0));
}

double DVonMises::l(std::vector<double const *> const & /*parameters*/) const {
  return 0.0;
}

double DVonMises::u(std::vector<double const *> const & /*parameters*/) const {
  return twoPi();
}

double DVonMises::randomSample(std::vector<double const *> const &par,
                               double const * /*lower*/,
                               double const * /*upper*/,
                               RNG *rng) const {
  const double mu = *par[0];
  const double kappa = *par[1];

  if (kappa == 0.0) {
    return twoPi() * rng->uniform();
  }

  double theta = 0.0;

  if (kappa <= 0.5) {
    // Siegerstetter rejection with uniform envelope.
    while (true) {
      const double u1 = twoPi() * rng->uniform();
      const double u2 = rng->uniform();
      if (std::log(u2) <= kappa * (std::cos(u1) - 1.0)) {
        theta = u1;
        break;
      }
    }
  } else {
    // Best-Fisher rejection with wrapped Cauchy envelope.
    const double tau = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    const double rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * kappa);
    const double r = (1.0 + rho * rho) / (2.0 * rho);

    while (true) {
      const double z = std::cos(M_PI * rng->uniform());
      const double f = (1.0 + r * z) / (r + z);
      const double c = kappa * (r - f);
      const double u2 = rng->uniform();

      if ((c * (2.0 - c) - u2) > 0.0 ||
          (std::log(c / u2) + 1.0 - c) > 0.0) {
        const double u3 = rng->uniform();
        theta = (std::floor(u3 + 0.5) * 2.0 - 1.0) * std::acos(f);
        break;
      }
    }
  }

  return wrapAngle(theta + mu);
}

double DVonMises::typicalValue(std::vector<double const *> const &par,
                               double const * /*lower*/,
                               double const * /*upper*/) const {
  return wrapAngle(*par[0]);
}

bool DVonMises::checkParameterValue(
    std::vector<double const *> const &par) const {
  // mu is treated as an angle on the circle and wrapped in logDensity.
  return (*par[1] >= 0.0);
}

bool DVonMises::isSupportFixed(std::vector<bool> const &fixmask) const {
  return fixmask[0] && fixmask[1];
}

}  // namespace vonmises
}  // namespace jags
