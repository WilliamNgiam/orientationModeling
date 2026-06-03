# JAGS von Mises module

Vendored JAGS extension (based on [yeagle/jags-vonmises](https://github.com/yeagle/jags-vonmises)) with local fixes to angle wrapping and sampling. Provides the von Mises distribution for circular data.

## Distribution

```bugs
x ~ dvonmises(mu, precision)
```

| Parameter   | Description                                                  | Support             |
| ----------- | ------------------------------------------------------------ | ------------------- |
| `mu`        | Mean direction (radians); any real value, wrapped internally | $(-\infty, \infty)$ |
| `precision` | Concentration `kappa`                                        | `[0, infinity)`     |

When `precision = 0`, the distribution is uniform on `[0, 2*pi)`.

Density (up to normalization):

```
f(x | mu, kappa) = exp(kappa * cos(x - mu)) / (2*pi * I_0(kappa))
```

Random variates use the Best-Fisher algorithm (Best & Fisher, 1979).

## Build and install

Requires JAGS development headers and libraries.

```bash
cd jags-vonMises
make
sudo make install
```

The module is installed to `$(libdir)/JAGS/modules-4/vonmises.so`.

## Usage

Load the module before defining the model:

```r
library(rjags)

model <- "
model {
  for (i in 1:N) {
    theta[i] ~ dvonmises(mu, kappa)
  }
  mu ~ dunif(0, 6.283185307179586)
  kappa ~ dgamma(0.1, 0.1)
}
"

jags.model(model, data = list(N = 10), modules = "vonmises")
```

From the JAGS command line (script file):

```
load vonmises
```

To load from a build directory without installing, set `JAGS_LIBS` to the directory
containing `vonmises.so`:

```bash
JAGS_LIBS=/path/to/jags-vonMises jags myscript.jag
```

## Example model

See `examples/vonmises.bug`.

## Circular mean estimation (important)

If you estimate a mean direction with a **uniform prior on `[0, 2π)`** (`mu ~ dunif(0, pi2)`), JAGS's default scalar samplers treat that interval as a line with hard endpoints. When the true mean is near 0, the posterior has equivalent mass near 0 and near `2π`, but a single chain can get stuck at one endpoint and never explore the other. This is the issue discussed in [this Cross Validated thread](https://stats.stackexchange.com/questions/459521/jags-circular-distribution-sampling-issues) for the yeagle module; it is a **sampler / parameterization issue**, not a bug in the Best–Fisher random generator.

**Recommended workarounds:**

1. **Initialize `mu` at the sample circular mean** (in radians).
2. Prefer a **wide unbounded prior** such as `mu ~ dnorm(0, 0.001)` instead of `dunif(0, pi2)`. The von Mises likelihood is periodic in `mu`, so this is usually more appropriate.
3. Run **multiple chains** with dispersed starting values.
4. Summarize with the **circular mean** of posterior samples, not the linear mean.

This module improves on the yeagle implementation by:

- Wrapping angles in `logDensity` and `randomSample` (yeagle's bare `fmod` can return negative values).
- Allowing any real `mu` in the density (wrapped internally), so periodic likelihood is consistent.

## References

- Forbes, M., Evans, M., Hastings, N., & Peacock, B. (2011). *Statistical Distributions*, 4th ed.
- Best, D. J., & Fisher, N. I. (1979). Efficient simulation of the von Mises distribution. *J. R. Statist. Soc. Ser. C*, 28, 152–157.
- Wabersich, D., & Vandekerckhove, J. (2014). Extending JAGS: A tutorial on adding custom distributions. *Behavior Research Methods*, 46, 15–28.
