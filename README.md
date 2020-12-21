# Run simulations testing variance propagation

In order to setup the simulations, several pieces of code are used to generate data, build models and then finally run the simulation. Several `RData` files are generated during the setup process.

The main simulation `scrubsim.R` will run all simulations in parallel where possible. This is *extremely* memory/processor hungry. We do not recommend just running all of these scripts in order without some thought about the potential memory footprint.

Due to parallelisation (and it's interaction with Nimble) some of the code here is written in a slightly convoluted way. We apologise for this!

File/folder descriptions:

- `proposed_detfcts.R` saves and plots the proposed detection functions to be used in the simulation. Outputs `df_pars.RData`
- `nimble_dsm.R` generates the objects needed by Nimble to fit the model.
- `scrubsim.R` runs the simulations (see above caveats about parallel code).
- `support_functions.R` includes various support functions for the simulation.
- `results_csvs/` results used to generate plots etc from our runs of the simulation, output from the above scripts happens in the top level directory, these results are kept separate to avoid being overwritten.

