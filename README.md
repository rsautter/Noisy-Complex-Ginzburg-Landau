# Stochastic Complex Ginzburg-Landau

Implementation of the Stochastic Complex Ginzburg-Landau (SCGL) using pseudospectral method and Runge-Kuta-Fehlberg 4-5 method.
The additive SCGL is:
    $$\partial_t A = (1+ib) \nabla^2 A + A  - (1+ic) |A|^2A + \sigma \partial_t(\eta_\beta)$$
The multiplicative SCGL is:   
    $$\partial_t A = (1+ib) \nabla^2 A + A  - (1+ic) |A|^2A + \sigma A \partial_t(\eta_\beta)$$
where $A$ is a complex number, and $\eta_\beta$ is a colored noise.
This method was implemented for a multidimensional context, meaning 1D, 2D and 3D examples are presented.

## Files

The implementation file is [this](https://github.com/rsautter/Noisy-Complex-Ginzburg-Landau/blob/main/NCGL.py). 

Some example of iteration is presented [here](https://github.com/rsautter/Noisy-Complex-Ginzburg-Landau/blob/main/CGL_Example.ipynb).

We are pleased to offer snapshots of the simulation upon request. 

## Noise Examples


## SCGL


## Multidimension CGL

Also available on [my youtube channel](https://youtu.be/tc63qzBIM7I)

https://user-images.githubusercontent.com/14216783/224053529-af028b73-6f58-47f4-a0a4-0108392fa167.mp4

## Paper link
Under submission
