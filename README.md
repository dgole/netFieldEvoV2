# netFieldEvoV2

This code evolves the net-magnetic-flux and surface density of an accretion disk while self-consistently solving a static disk model at each point in time.  For a full description of the model and it's implementation, see chapter 3 of my thesis work: https://dgole.github.io/pages/thesis.html.   

# Organization
fmatrix: Contains grid geometry files and tables of pre-calculated geometry terms (the "f matricies").  

mathematica: Contians code to set up grids and calculate the f matricies.  Also contains lots of random notebooks plotting and calculating various disk models.  

python: The python code-base to run the model and plot the results.  

scripts_run: Scripts (bash) that run the models.  

scripts_analysis: Scripts (bash and python) that make plots.  

# Example Plots
Time evolution plot:  
<img src="./demo_spaceTime.png" width="651" height="578" />  

Array of light curves for various disk/star parameters:  
<img src="./demo_lightCurves.png" width="560" height="597" />  
