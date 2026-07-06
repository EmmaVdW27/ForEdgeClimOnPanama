# ForEdgeClimOnPanama

**ForEdgeClimOnPanama** is an R package for modelling forest microclimates at high spatial resolution (~1 m). The model simulates key surface energy fluxes, including shortwave and longwave radiation, latent and sensible heat fluxes, and ground heat fluxes.

The package is based on **ForEdgeClim v1.0**, a 3D process-based microclimate model originally developed to simulate microclimate gradients along forest edge-to-core transects:

> Van de Walle et al. (2026). *ForEdgeClim v1.0: a 3D process-based microclimate model incorporating vertical and lateral radiative and thermal fluxes to simulate forest edge-to-core transitions*. Geoscientific Model Development. https://gmd.copernicus.org/articles/19/4661/2026/

This repository contains a modified version of the original model that is specifically designed for simulations in **forest interiors**, without explicitly accounting for forest edges.

## Installation

Install the package directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("EmmaVdW27/ForEdgeClimOnPanama")
```

## Getting started

An example demonstrating how to run **ForEdgeClimOnPanama** is provided in:

```text
testing_foredgeclim.R
```


