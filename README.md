# Python Metacommunity Models (pyMCME)

This library aims to implement the metacommunity model outlined within Leibhold and Chase (2018) in an efficient, python-based, framework to later study community dynamics using the Haar flucutation wavelet outlined by Lovejoy (2015).

This model numerically simulates a network of habitat patches (or local communities) that are connected by dispersal. The metacommunity framework includes all of the normal archetypes as described by Leibhold and Chase (2018), with the addition of trait evolution and speciation.  

![Key Metacommunity Model Processes](./main/images/eco-evo-fig.png)

## Model Framework

At the core of this general simulation model is an array of spatially explicit habitat patches with a dynamic environment that can change over a pre-defined temporal extent and resolution. As such, the user must provide information regarding the number of patches, their spatial orientation, and changes in it's climatic conditions. This model was built to handle HadCM3 climate data imported as a raster stack. However, data can be imported from a NetCDF4 file, but assurances need to be made that the spatial resolution matches that of a standard terrestrial HadCM3 simulation (3.75 by 2.5).   

The model, which will be referred to as pyMCME, operates through four major processes: 1) species interactions; 2) environmental (abiotic) effects; 3) dispersal processes; and 4) speciation and evolution.

### Population Dynamics

#### Species Interactions

Species interactions are determined according the the Beaverton-Holt equation as described in Thompson et al. (2020). Species interaction strengths are assigned using the MCME_Functions.initalize_aij function whereby the interaction category needs to be defined (stabalizing, equal, mixed, neutral) as well as the numerical range of the interactions for the stabalizing and mixed interaction categories. For the mixed category, a percentage of competitively superior species needs to be provided, which is ignored in all of the other catagories. During speciation events, new species interaction coefficients are redrawn from the range provided within the initialization stage. Furthermore, interaction strengths are multiplied by 0.05 to increase interaction rates.

```math
N_{ix}(t+1) = \frac{1}{1 + \sum{a_{ix}N_{jx}(t)}}
```

Future work will be to include in tradeoffs along multiple trait axes, including a trait axes related to dispersal ability.

<img src="./main/images/git_interactions_fig.png" width=50% height=50%>

#### Environmental (abiotic) Effects

Growth rate within pyMCME is determined by the match between a species niche optimum and the current environmental conditions. Currently, a Gaussian function is used to determine the relationship between the environment and growth. For this function a maximum growth rate needs to be selected, which is defaulted to 5 to speed up processes within the model. Additionally, a sigma value needs to be provided which describes the niche breadth for a species (see figure 2).

* Future work will be to convert the Gaussian function to a Gompertz Gaussian function. I also aim to incorporate the metabolic theory of ecology (Brown et al.) into the model, and an effect of temperature on dispersal.

** Growth Rate Exponent **
```math
r_{ix}(t) = r_{max}e^{-(\frac{z_{i} - env_{x}(t)}{2\sigma_{i}})^2}
```

<img src="./main/images/git_growth_fig.png" width=50% height=50%>

#### Demographic Stochasticity

Following Thompson et al. (2020), demographic stochasticity is introduced into the model via a random population integer value drawn from a Poisson distribution prior to dispersal process by the following equation:

```math
N_{ix}(t+1) = Poisson(max{\frac{1}{1 + \sum{a_{ix}N_{jx}(t)}}, 0})
```

### Dispersal

After each time-step individuals from within a habitat patch can disperse across the landscape based upon their density within each patch and a pre-defined dispersal rate. Dispersal distance follows a Poisson distribution according to the below equation. A probability of dispersal is then assigned based upon the density of species within a patch.

```math
I_{ix}(t) = \frac{\sum{E_{iy}(t)^{-L_{i}d_{x}}}}{\sum{E_{ix}(t)}}
```

Future work will look to move beyond a stochastic dispersal process and instead link dispersal probability to the environment as well as to the strength of species interactions. This may be an ideal process to include as a tradeoff.    

The full equation that determines population abundance within each patch is as follows, with the a Poisson value drawn from the right side of the equation, excluding immigration and emigration.

```math
N_{ix}(t+1) = r_{ix}N_{ix}(t)\frac{1}{1 + \sum{a_{ix}N_{jx}(t)}}-E_{ix}(t) + I_{ix}(t)

```

<img src="./main/images/git_dispersal_fig.png" width=50% height=50%>

### Trait Evolution and Speciation

Trait evolution within pyMCME is very simple whereby a trait will evolve via random drift based upon a defined standard deviation for a random value drawn from a uniform distribution. The default standard deviation for trait evolution is 0.05.

** Trait Evolution **
```math
z_{ix}(t+1) = z_{ix}(t) + U(0.05)
```

* Future work could look to alter the rate of trait evolution as a tradeoff with dispersal ability or species interactions.

The most complicated aspect of this model, and the area by which it varies the most from other models, is via speciation. To be clear, pyMCME does not adhere to strict allopatic speciation based upon patch isolation from the main species pool. Instead, speciation occurs when a local communities (single population within a patch) mean niche optimum extends beyond the standard deviation of the entire species (whole metacommunity excluding the patch of interest) plus some threshold value for speciation. This formulation of speciation allows for gene flow to occur due to dispersal processes, but also allows for speciation to occur via random drift.

** Speciation **
```math
\mu_{ix}^{M} > \sigma_{ix} + \omega_{x}
```

* Future work will aim to include increased speciation rates that are tied to environmental conditions. Increased speciation due to warming temperatures.

<img src="./main/images/git_speciation_fig.png" width=50% height=50% align:center>

## Simulation Initialization

Simulations are initialized at the beginning of each experiment based upon a starting number of patches and species. Currently, habitat patches are derived from the HadCM3 climate model simulations, with graph distances calculated based on each gridcells latitude and longitude positions. Habitat patches are drawn randomly from a total map extent. A user can define the desired patch coordinates for the HadCM3 climate inputs.

* Future work will aim to accommodate the CHELSA-Trace21K climate ensemble data.

Users must define:
- Number of starting species
- Maximum growth rate
- Number of patches
- Geographic positions of patches
- Climate input
- Dispersal rate
- Niche breadth
- Matrix of niche optimums for start of simulation
- Speciation threshold
- Species interaction type
- Alpha matrix
- Simulation time (derived from seed time, burn in time, and total steps provided within the climate input)

Simulation output:
- MxS matrix of patch occupancy
- MxS matrix of niche optimums for species per patch
- MxS matrix of species interaction coefficients
- Phylogeny (anscestor/descendent)
- Divergence time corresponding to speciation occurrences
- Species origination times in patches
