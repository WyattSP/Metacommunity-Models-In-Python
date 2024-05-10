# Metacommmunity Models in Python

This library aims to implement metacommunity models outlined within Leibhold and Chase (2018) in an efficient, python-based, framework to later study community dynamics using the Haar flucutation wavelet outlined by Lovejoy (2015).

This model aims to be able to numerically simulate a network of habitat patches (or local communities) that are either connected by dispersal to each other or to a metacommunity that can evolve through time. The metacommunity frameworks include: 1) Patch Dynamics, 2) Species Sorting, 3) Mass Effect, and 4) Neutral Dynamics.

- Diagram illustrating model geometry

## Model Framework

### Landscape Initialization
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


### Simulation Loop
