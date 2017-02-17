# RainMaker


[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![Build Status](https://travis-ci.org/granttremblay/rainmaker.svg?branch=master)](https://travis-ci.org/granttremblay/rainmaker)
[![Coverage Status](https://coveralls.io/repos/github/granttremblay/rainmaker/badge.svg)](https://coveralls.io/github/granttremblay/rainmaker)

[Dr. Grant R. Tremblay](www.granttremblay.com) | Einstein Fellow | [Yale University](www.yale.edu) | *grant.tremblay @ yale.edu*
___
### Predicting precipitation in galaxy clusters
`rainmaker` is a python implementation of an `IDL` code written by Prof. G. Mark Voit and Prof. Megan Donahue at Michigan State University.

It works by fitting log temperature and density profiles to the [ACCEPT](http://www.pa.msu.edu/astro/MC2/accept/) table of cluster entropy profiles, originally compiled by [Cavagnolo et al. (2009)](https://ui.adsabs.harvard.edu/?#abs/2009ApJS..182...12C).
