# TAMW-codes
MATLAB scripts for the geometric analysis accompanying the TAMW manuscript.
# TAMW: Analysis of the Geometric Parameter λ

This repository contains MATLAB scripts used for the numerical and graphical analysis of the geometric parameter λ in the Total Active Myocardial Work (TAMW) framework.

The analysis explores the dependence of λ on left ventricular wall volume and end-diastolic volume, as well as the influence of their statistical distributions and correlation on the resulting variability of λ.

## Mathematical background

The geometric parameter is defined as:

```math
\lambda = 1 + \frac{3V_{\mathrm{ED}}}{V_w}
```

where:
- $V_{\mathrm{ED}}$ is the left ventricular end-diastolic volume.
- $V_w$ is the left ventricular wall volume.
The analysis considers these parameters as normally distributed random variables, initially assuming independence and subsequently introducing different degrees of correlation.

The probability calculations are restricted to specified physiological parameter ranges, with normalization over the resulting domain.

## Repository contents

### Figure 1 — Geometric and probabilistic visualization

Located in `figure1/`.

Four independent MATLAB scripts generate the subfigures of manuscript Figure 1:

- `lambda_3d_figure1_gaussian_marginals.m` — Joint Gaussian probability density and its marginal distributions.
- `lambda_3d_figure2_gaussian_shadow.m` — Projection of the joint probability density onto the parameter plane.
- `lambda_3d_figure3_lambda_surface.m` — Three-dimensional representation of the geometric parameter λ.
- `lambda_3d_figure4_lambda_band_projection.m` — Visualization of a relative λ interval, its geometric boundaries, and its projection onto the probability distribution.

For visualization purposes, the fourth script exaggerates vertical separation without changing the underlying definition of λ.

### Figure 2 — Relative variability of λ

Located in `figure2/`.

`lambda_relative_error_80.m`

This script evaluates the probability associated with relative intervals around λ under the assumption of independent Gaussian distributions of ventricular wall volume and end-diastolic volume.

It identifies the value λ_peak maximizing the probability of a ±5% relative interval and determines the relative interval width δ₈₀ required to encompass 80% of the modeled probability within the specified parameter domain.

The script generates graphical representations of the probability function, the dependence of interval probability on relative width, and the resulting 80% interval.

### Figure 3 — Influence of correlation

Located in `figure3/`.

`lambda_relative_error_80_correlated_figure3.m`

This script extends the preceding analysis to correlated bivariate Gaussian distributions.

Calculations are performed for correlation coefficients:

$$
\rho = 0,\ 0.25,\ 0.50,\ 0.75
$$

The script calculates the corresponding λ_peak and δ₈₀ values and generates a comparison of the resulting 80% intervals.

## Parameters and numerical methods

The calculations use the following distribution parameters and bounds:

| Parameter | Value |
|---|---|
| Mean LV wall volume | 97 mL |
| SD of LV wall volume | 30 mL |
| Mean LV end-diastolic volume | 139 mL |
| SD of LV end-diastolic volume | 31 mL |
| LV wall volume bounds | 19–195 mL |
| LV end-diastolic volume bounds | 57–245 mL |

Numerical evaluation uses MATLAB integration routines, bounded optimization, and root-finding methods.

The probability analysis is performed within the specified bounded parameter domain. Its results are conditional on the assumed distributions, correlations, and parameter ranges.

The theoretical derivation and rationale for the parameter choices are described in the accompanying manuscript.

## Software requirements

- MATLAB R2026a

No external datasets are required. Input parameters are defined directly within the scripts.

The scripts are standalone and can be executed individually in MATLAB. Figures are generated in the MATLAB graphical environment; the scripts do not automatically export image files.

## Reproducibility

To reproduce a figure:

1. Download the corresponding MATLAB script.
2. Open the file in MATLAB R2026a.
3. Run the script.
4. Compare the resulting output with the corresponding manuscript figure.

The numerical scripts also print relevant calculated values in the MATLAB Command Window.

## Code development

The MATLAB implementation was developed with coding assistance from OpenAI Codex. The mathematical formulation, analytical design, and scientific interpretation were determined by the authors.

## Related publication

This repository accompanies the manuscript presenting the theoretical and numerical analysis of TAMW.

A publication reference and DOI will be added when available.
