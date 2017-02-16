# Rainmaker


[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![Build Status](https://travis-ci.org/granttremblay/rainmaker.svg?branch=master)](https://travis-ci.org/granttremblay/rainmaker)
[![Coverage Status](https://coveralls.io/repos/github/granttremblay/rainmaker/badge.svg?branch=master)](https://coveralls.io/github/granttremblay/rainmaker?branch=master)

[Dr. Grant R. Tremblay](www.granttremblay.com) | Einstein Fellow | [Yale University](www.yale.edu) | *grant.tremblay @ yale.edu*
___
### Predicting precipitation in galaxy clusters
`rainmaker` is a python implementation of an `IDL` code written by Prof. G. Mark Voit and Prof. Megan Donahue at Michigan State University.

It works by fitting log temperature and density profiles to the [ACCEPT](http://www.pa.msu.edu/astro/MC2/accept/) table of cluster entropy profiles, originally compiled by [Cavagnolo et al. (2009)](https://ui.adsabs.harvard.edu/?#abs/2009ApJS..182...12C).

Because
$t_\mathrm{cool} \sim \frac{kT}{n\Lambda(T)}$

and

$t_\mathrm{FreeFall} \sim \frac{r}{\left( kT \times \left|\frac{d\ln{P}}{d\ln{r}}\right|\right)^{1/2}},$

the ratio between the two goes as

$\frac{t_\mathrm{cool}}{t_\mathrm{FreeFall}} \sim \left(\frac{1}{nr} \right ) \times \left( kT \right)^{3/2} \times [\Lambda(T)]^{-1} \times \left|\frac{d\ln{P}}{d\ln{r}}\right|^{1/2}.$

At radii of several tens of kpc in a cool-core cluster,
the product of density times radius is roughly constant,
so the value of this ratio will track the temperature gradient.
However, we find that the density profile at smaller radii flattens
out in our deprojected profiles, which causes the product nr to
grow with radius and the timescale ratio to drop with radius.

So the key thing to look at is the behavior of the product of
radius and density at small radii.
