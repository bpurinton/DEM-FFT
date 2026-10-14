---
title: "Fourier Frequency Analysis of DEM"
author: "Ben Purinton ([purinton@uni-potsdam.de](purinton@uni-potsdam.de))"
date: "May 2020"
---

# Python code for DEM noise analysis using 2D DFT

This repository is intended for the spectral analysis of gridded topographic digital elevation models (DEMs) as presented in:

* Purinton, B., and Bookhagen, B.: Validation of digital elevation models (DEMs) and geomorphic metrics on the southern Central Andean Plateau, Earth Surface Dynamics, 2017. [https://doi.org/10.5194/esurf-5-211-2017](https://doi.org/10.5194/esurf-5-211-2017)

A walk through example is provided in the ipython notebook ```example_analysis.ipynb```

*Note: This code was originally developed as a MATLAB<sup>TM</sup> repository at: [https://github.com/bpurinton/DEM_fourier_noise](https://github.com/bpurinton/DEM_fourier_noise)*

# Original source

This code was largely translated from Taylor Perron's ```2DSpecTools``` package for MATLAB<sup>TM</sup> available at [http://web.mit.edu/perron/www/downloads.html](http://web.mit.edu/perron/www/downloads.html)

For background on this spectral analysis procedure and the paper that spurred this analysis refer to: 

* Perron, J. T., Kirchner, J. W., and Dietrich, W. E.: Spectral signatures of characteristic spatial scales and nonfractal structure in landscapes, Journal of Geophysical Research, 113, 2008. [https://doi.org/10.1029/2007JF000866](https://doi.org/10.1029/2007JF000866)


# Spatial recurrence lacunarity

As a complement to the Fourier analysis, ```recurrence_lacunarity_analysis.ipynb``` applies the normalized recurrence lacunarity (nRL) of spatial recurrence plots to slope maps of the same two DEMs, using the functions in ```recurrence_lacunarity.py```. The method follows:

* Bielig, C., Braun, T., Palmero, M. S., Rheinwalt, A., and Marwan, N.: Recurrence lacunarity of spatial data, The European Physical Journal Special Topics, 2026. [https://doi.org/10.1140/epjs/s11734-026-02569-4](https://link.springer.com/article/10.1140/epjs/s11734-026-02569-4)
