# uot4ens
Unbalanced Optimal Transport for Ensemble Forecast

This library contains data and Python scripts used to generate figures for the article:

Geometry of rainfall ensemble means: from arithmetic averages to Gaussian-Hellinger barycenters in unbalanced optimal transport

by Le Duc and Yohei Sawada in Journal of the Meteorological Society of Japan.

There are two types of Python scripts:
1. The Sinkhorn-Knoop algorithms are implemented into: uot1D_sinkhornlog.py for one-dimensional distributions; uot2D_sinkhornlog.py for two-dimensional distributions
2. Figure generation Python scripts: they have names 01 to 08 *py corresponding to 8 figures in the article.

There are two types of data: here 02 means LETKF1000, 04 means LETKF100, 100 means MEPS.
1. Rainfall time series: Ichifusa*nc files
2. Rainfall distributions: Kyushu*nc files

