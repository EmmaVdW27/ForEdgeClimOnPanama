# ForEdgeClim

ForEdgeClim is an R package for modelling microclimates in forests, including radiation processes and heat transfer.
The model is initially written to simulate microclimate gradients along transect lines from a forest’s core towards its edge. However it can also be applied on forest interiors only.
This repository is build upon the code used in the manuscript 'Simulating forest edge-to-core microclimate with ForEdgeClim v1.0: a process-based model incorporating vertical and lateral energy fluxes' and a more detailed overview with the applied physical equations can be found therein (DOI:).

In this version of ForEdgeClim, ForEdgeClimOnPanama, some changes have been made to the original (v1) version of ForEdgeClim, for the model only to run on forest interiors, not taking edges into account. The changes are enlisted below:

- no sensible heat exchange (diffusion/convection/turbulence) with the edge macroenvironment in the function calculate_D()
- update weighting factors in the function run_foredgeclim() for the linearization of the air temperature
- 2D shortwave two-stream RTM is updated towards 1D in shortwave_two_stream_RTM(), parameters linked to lateral radiation are removed - same for longwave RTM
- change ForEdgeClim to ForEdgeClimOnPanama in NAMESPACE, init.R, .dll file in src map en RcppExports.cpp
- copy the DESCRIPTION file and write ForEdgeClimOnPanama instead of ForEdgeClim in package name
- remove _2D in testing_foredgeclim.R
- remove _h related variables in plots_sw(), plots_lw(), plots_temp() en plots_flux() -> mainly captions
- update import_observations() with input dataframes of interest
- update create_input_drivers()


