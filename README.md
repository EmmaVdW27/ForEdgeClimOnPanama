# ForEdgeClim

ForEdgeClim is an R package for modelling microclimates in forests, including radiation processes and heat transfer.
The model is initially written to simulate microclimate gradients along transect lines from a forest’s core towards its edge.
This repository contains the code used in the manuscript 'ForEdgeClim v1.0: a 3D process-based microclimate model incorporating vertical and lateral energy fluxes to simulate forest edge-to-core transitions' and a more detailed overview with the applied physical equations can be found therein (submitted to Geoscientific Model Development).

In this version of ForEdgeClim, ForEdgeClimOnPanama, some changes have been made to the original (v1.0) version of ForEdgeClim, for the model only to run on forest interiors, not taking edges into account.

## Installation
Install the package directly from GitHub:
```r
# Using devtools
install.packages("devtools")
devtools::install_github("EmmaVdW27/ForEdgeClimOnPanama")
```

## Example functions
An example script to show you how the microclimate model ForEdgeClimOnPanama can be run is given in 'testing_foredgeclim.R'.


