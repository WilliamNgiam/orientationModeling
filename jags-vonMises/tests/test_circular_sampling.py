#!/usr/bin/env python3
"""Regression tests for circular boundary / wrapping issues.

Scenario from:
https://stats.stackexchange.com/questions/459521/jags-circular-distribution-sampling-issues
"""

from __future__ import annotations

import math
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.stats import vonmises

MODULE_DIR = Path(__file__).resolve().parent.parent
JAGS_LIBS = str(MODULE_DIR)
TWO_PI = 2.0 * math.pi


def run_jags(
    model: str,
    data_lines: list[str],
    *,
    n_iter: int = 5000,
    burnin: int = 500,
    inits_lines: list[str] | None = None,
    monitor: str = "mu_d",
    stem: str = "out",
) -> np.ndarray:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        model_file = tmp_path / "model.bug"
        data_file = tmp_path / "data.txt"
        script_file = tmp_path / "script.jag"
        model_file.write_text(model)
        data_file.write_text("\n".join(data_lines) + "\n")

        init_block = ""
        if inits_lines:
            inits_file = tmp_path / "inits.txt"
            inits_file.write_text("\n".join(inits_lines) + "\n")
            init_block = f'parameters in "{inits_file}"\n'

        script = f"""load vonmises
model in "{model_file}"
data in "{data_file}"
compile
{init_block}initialize
update {burnin}
monitor {monitor}
update {n_iter}
coda {monitor}, stem("{tmp_path / stem}")
"""
        script_file.write_text(script)
        env = os.environ.copy()
        env["JAGS_LIBS"] = JAGS_LIBS
        result = subprocess.run(
            ["jags", str(script_file)],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "JAGS failed\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )

        chain_file = tmp_path / f"{stem}chain1.txt"
        samples = []
        for line in chain_file.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    samples.append(float(parts[1]))
                except ValueError:
                    pass
        return np.array(samples)


def circular_mean_deg(samples_deg: np.ndarray) -> float:
    rad = np.deg2rad(samples_deg)
    return math.degrees(math.atan2(np.mean(np.sin(rad)), np.mean(np.cos(rad))))


def mass_near_zero_deg(samples_deg: np.ndarray, threshold_deg: float = 45.0) -> tuple[float, float]:
    low = np.sum(samples_deg <= threshold_deg) / len(samples_deg)
    high = np.sum(samples_deg >= (360.0 - threshold_deg)) / len(samples_deg)
    return low, high


def sample_circular_mean_rad(x: np.ndarray) -> float:
    return math.atan2(np.mean(np.sin(x)), np.mean(np.cos(x)))


def main() -> None:
    rng = np.random.default_rng(42)
    x = vonmises.rvs(4.0, loc=0.0, size=50, random_state=rng)
    x = np.mod(x, TWO_PI)
    data_lines = ['"N" <- 50', '"x" <- c(' + ", ".join(f"{v:.16g}" for v in x) + ")"]

    print("=== 1. randomSample quality (x ~ dvonmises(0, 4)) ===")
    prior_samples = run_jags(
        """model {
  for (i in 1:N) {
    x[i] ~ dvonmises(0, 4)
  }
}""",
        ['"N" <- 2000'],
        n_iter=2000,
        burnin=100,
        monitor="x",
        stem="prior",
    )
    invalid = np.sum((prior_samples < 0) | (prior_samples >= TWO_PI))
    prior_deg = np.rad2deg(prior_samples)
    low, high = mass_near_zero_deg(prior_deg, 30.0)
    print(f"  invalid samples: {invalid} / {len(prior_samples)}")
    print(f"  near 0°: {low:.2%}  near 360°: {high:.2%}  (scipy ~34% each)")
    assert invalid == 0, "randomSample produced out-of-range values (yeagle fmod bug)"
    assert 0.20 < low < 0.50 and 0.20 < high < 0.50, "asymmetric prior simulation"
    print("  OK")

    so_model = """model {
  pi2 <- 3.141592653589793 * 2
  mu_d ~ dunif(0, pi2)
  for (i in 1:N) {
    x[i] ~ dvonmises(mu_d, 4)
  }
}"""

    print("\n=== 2. SO scenario: dunif prior, default init (problematic) ===")
    mu_bad = run_jags(so_model, data_lines, n_iter=5000, burnin=500)
    mu_bad_deg = np.rad2deg(mu_bad)
    low, high = mass_near_zero_deg(mu_bad_deg, 45.0)
    print(f"  near 0°: {low:.2%}  near 360°: {high:.2%}")
    print(f"  circular mean: {circular_mean_deg(mu_bad_deg):.1f}°")
    one_sided = (low > 0.15 and high < 0.05) or (high > 0.15 and low < 0.05)
    if one_sided:
        print("  EXPECTED: one-sided boundary sticking with dunif + default init")
    else:
        print("  NOTE: both sides explored (may depend on RNG)")

    print("\n=== 3. dunif prior, init at sample circular mean (recommended) ===")
    cm = sample_circular_mean_rad(x)
    mu_good = run_jags(
        so_model,
        data_lines,
        inits_lines=[f'"mu_d" <- {cm:.16g}'],
        n_iter=5000,
        burnin=500,
    )
    mu_good_deg = np.rad2deg(mu_good)
    low, high = mass_near_zero_deg(mu_good_deg, 45.0)
    cm_post = circular_mean_deg(mu_good_deg)
    print(f"  near 0°: {low:.2%}  near 360°: {high:.2%}")
    print(f"  circular mean: {cm_post:.1f}° (data ~ {math.degrees(cm):.1f}°)")
    assert abs(cm_post) < 30 or abs(cm_post - 360) < 30, "circular mean far from 0"
    print("  OK")

    print("\n=== 4. wide normal prior (recommended) ===")
    mu_norm = run_jags(
        """model {
  mu_d ~ dnorm(0, 0.001)
  for (i in 1:N) {
    x[i] ~ dvonmises(mu_d, 4)
  }
}""",
        data_lines,
        inits_lines=[f'"mu_d" <- {cm:.16g}'],
        n_iter=5000,
        burnin=500,
    )
    cm_norm = circular_mean_deg(np.mod(np.rad2deg(mu_norm), 360.0))
    print(f"  circular mean: {cm_norm:.1f}°")
    assert abs(cm_norm) < 30 or abs(cm_norm - 360) < 30
    print("  OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
