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

The following video shows the  traditional Complex Ginzburg-Landau:

https://user-images.githubusercontent.com/14216783/224129416-69f13958-8244-4741-9cf3-d6b9e51144e2.mp4


## SCGL
The additive noise iteration example is:

https://user-images.githubusercontent.com/14216783/224137213-3b0a5406-c470-48c6-b131-187e2313415a.mp4


## Noise Examples


## Multidimension CGL

The multidimensional case is under development. The solution for the simplest case is:

https://user-images.githubusercontent.com/14216783/224053529-af028b73-6f58-47f4-a0a4-0108392fa167.mp4

## Paper link
Under submission
