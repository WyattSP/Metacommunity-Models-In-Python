# Python Metacommunity Models (pyMCME)

This library aims to implement the metacommunity model outlined within Leibhold and Chase (2018) in an efficient, python-based, framework to later study community dynamics using the Haar flucutation wavelet outlined by Lovejoy (2015).

This model numerically simulates a network of habitat patches (or local communities) that are connected by dispersal. The metacommunity framework includes all of the normal archetypes as described by Leibhold and Chase (2018), with the addition of trait evolution and speciation.  

## Model Framework

At the core of this general simulation model is an array of spatially explicit habitat patches with a dynamic environment that can change over a pre-defined temporal extent and resolution. As such, the user must provide information regarding the number of patches, their spatial orientation, and changes in it's climatic conditions. This model was built to handle HadCM3 climate data imported as a raster stack. However, data can be imported from a NetCDF4 file, but assurances need to be made that the spatial resolution matches that of a standard terrestrial HadCM3 simulation (3.75 by 2.5).   

The model, which will be referred to as pyMCME, operates through four major processes: 1) species interactions; 2) environmental (abiotic) effects; 3) dispersal processes; and 4) speciation and evolution.

### Species Interactions

Species interactions are determined according the the Beaverton-Holt equation as described in Thompson et al. (2020). Species interaction strengths are assigned using the MCME_Functions.initalize_aij function whereby the interaction category needs to be defined (stabalizing, equal, mixed, neutral) as well as the numerical range of the interactions for the stabalizing and mixed interaction categories. For the mixed category, a percentage of competitively superior species needs to be provided, which is ignored in all of the other catagories. During speciation events, new species interaction coefficients are redrawn from the range provided within the initialization stage. Furthermore, interaction strengths are multiplied by 0.05 to increase interaction rates.

* Future work will be to include in tradeoffs along multiple trait axes, including a trait axes related to dispersal ability.

** Beaverton-Holt Equation **
```math
N_{ix}(t+1) = \frac{1}{1 + \sum_{j = 1}^{S}{a_{ix}N_{jx}(t)}}
```

### Environmental (abiotic) Effects

Growth rate within pyMCME is determined by the match between a species niche optimum and the current environmental conditions. Currently, a Gaussian function is used to determine the relationship between the environment and growth. For this function a maximum growth rate needs to be selected, which is defaulted to 5 to speed up processes within the model. Additionally, a sigma value needs to be provided which describes the niche breadth for a species (see figure 2).

* Future work will be to convert the Gaussian function to a Gompertz Gaussian function. I also aim to incorporate the metabolic theory of ecology (Brown et al.) into the model, and an effect of temperature on dispersal.

** Growth Rate Exponent **
```math
r_{ix}(t) = r_{max}e^{-(\frac{z_{i} - env_{x}(t)}{2\sigma_{i}})^2}
```

### Dispersal

After each time-step individuals from within a habitat patch can disperse across the landscape based upon their density within each patch and a pre-defined dispersal rate. Dispersal distance follows a Poisson distribution according to the below equation. A probability of dispersal is then assigned based upon the density of species within a patch.

* Future work will look to move beyond a stochastic dispersal process and instead link dispersal probability to the environment as well as to the strength of species interactions. This may be an ideal process to include as a tradeoff.    

** Dispersal Distance **
```math

```

The full equation that determines population abundance within each patch is as follows:
```math
N_{ix}(t+1) = r_{ix}N_{ix}(t)\frac{1}{1 + \sum_\{j = 1}^{S}{a_{ix}N_{jx}(t)}}-E_{ix}(t) + I_{ix}(t)
```


### Trait Evolution and Speciation

Trait evolution within pyMCME is very simple whereby a trait and evolve via random drift based upon a defined standard deviation for a random value to be drawn from a uniform distribution. The default standard deviation for trait evolution is 0.05.

* Future work could look to alter the rate of trait evolution as a tradeoff with dispersal ability or species interactions.

The most complicated aspect of this model, and the area by which it varies the most from other models, is the way in which species speciate. To be clear, pyMCME does not adhere to strict allopatic speciation based upon isolation time from the main species population. Instead, speciation occurs when a local communities (single population within a patch) mean niche optimum extends beyond the standard deviation of the entire species (whole metacommunity excluding the patch of interest) plus some threshold value for speciation. This formulation of speciation allows for gene flow to occur due to dispersal processes, but also allows for speciation to occur via random drift.

* Future work will aim to include increased speciation rates that are tied to environmental conditions. Increased speciation due to warming temperatures.

## Simulation Initialization
User functionality:
- User selected number of local communities
- User selected inclusion of a metacommunity
- User selected distance between local communities

Input:
- User functionality is input via raster with correct landscape geometries or through a 2-dimensional matrix or array

Internal:
- All raster files will be converted to sparse matrixes for efficient calculations
- Algorithm description:
  1) User provides a raster or matrix with demarcated habitat versus non-habitat. This is converted to a binary matrix for now. Could include elevation at some point.
  2) A distance is calculated between all habitat patches. Distance is a single value to different islands stored within a vector indexed by local community index.
  3) User can assign a edge direction off map for a metacommunity, it can be randomly assigned, or it can not be included.
  4)
