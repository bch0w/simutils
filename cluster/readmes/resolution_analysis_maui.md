# Resolution Analysis on Maui

Bryant Chow 
5/15/23

The following document describes how to run the resolution analysis performed in Chow et al. (2022a,b) 
using SPECFEM3D_Cartesian on Maui. These analyses include a Zeroth Moment test and Point Spread
Function test, outlined in various publications from Andreas Fichtner et al. (see background).

## Background

The full theory and mathematical background of these resolution studies can be found in 
[Fichter et al. (2011)](https://academic.oup.com/gji/article/187/3/1604/616815)
and [Fichter et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012106).

The general idea here is that we are approximating the Hessian (derivative of the gradient) to get
some information about resolution of given features.

"The Hessian represents our blurred perception of a point-localized perturbation at position 
y in a linearized tomographic inversion. The effect of the off-diagonal elements Hij|i≠j is 
to introduce unwanted updates of model parameters mj|i≠j that have initially not been perturbed." 
-- Fichter et al. (2011)

Mathematically, we are calculating the "action Hessian", which is the finite difference between two
gradients. That is, the Hessian `H(m)` can be calculated as `H(m) = g(m + dm) - g(m)`, where `g(m)`
is the gradient of your final model, and `g(m + dm)` is the gradient calculated using a perturbed
version of your final model. The choice of perturbation `dm` determines what type of test we 
are running.

Practically, this involves N additional simulations, where N is the total number of events. 
You should already have `g(m)` from your final simulation run. That means we just need to 
calculate `g(m+dm)` by perturbing our final model `m`, and calculating the gradient 
`g(m+dm)`. We then take the difference of the two gradients to retrieve `H(m)`.


## Setup

We will use SeisFlows to run the test. Start a **new** working directory (resolution) and:

1) Copy gradient files from final model to: *resolution/global/eval1/gradient*
2) Create empty directory: *resolution/global/eval2/gradient*
3) 

## Zeroth Moment Test (Volumentric Point Spread Function)

"The zeroth moment, M(0)(x), is equal to the integral over the PSF for a point perturbation at position x. 
It follows that M(0)(x) is small when the misfit χ is nearly unaffected by a point perturbation at x. In 
contrast, M(0)(x) is large when the PSF has a high amplitude, a large spatial extent, or both. In this sense, 
the spatial variability of M(0)(x) reflects the relative weight of neighbouring PSFs. Information on the 
resolution length is not contained in the zeroth moment."

In the Zeroth Moment test our perturbation `dm` is a homogeneous offset to all values in our final model `m`.

**Amplitude of Perturbation**: "...travel times measured by cross correlation are linearly related to seismic velocity 
perturbations of up to 10% [Mercerat and Nolet, 2013], meaning that **random velocity perturbations in the percent range 
would be a meaningful choice.**


1. Place your gradient `g(m)` (scratch/evalgrad/gradient/\*.bin) 


## Point Spread Function

**Size of PSF**: "Most of the requirements imposed by numerical methods can be met by choosing the width of the random 
perturbations to be **around half a wavelength**... Wider perturbations, for instance, in response to numerical 
restrictions, will generally increase the estimated resolution lengths, thus making them more conservative."

