# Stochastic Complex Ginzburg-Landau

Implementation of the Stochastic Complex Ginzburg-Landau (SCGL) using pseudospectral method and Runge-Kuta-Fehlberg 4-5 method.
The additive SCGL is:
    $$\partial_t A = (1+ib) \nabla^2 A + A  - (1+ic) |A|^2A + \sigma \partial_t(\eta_\beta)$$
The multiplicative SCGL is:   
    $$\partial_t A = (1+ib) \nabla^2 A + A  - (1+ic) |A|^2A + \sigma A \partial_t(\eta_\beta)$$
where $A$ is a complex number, and $\eta_\beta$ is a colored noise.
This method was implemented for a multidimensional context, meaning 1D, 2D and 3D examples are presented.

## Files

The implementation of SCGL is given by  [NCGL.py](https://github.com/rsautter/Noisy-Complex-Ginzburg-Landau/blob/main/NCGL.py). 

The noise generation algorithm is given by the [cNoise.py](https://github.com/rsautter/Noisy-Complex-Ginzburg-Landau/blob/main/cNoise.py)

The examples are presented [here](https://github.com/rsautter/Noisy-Complex-Ginzburg-Landau/blob/main/CGL_Example.ipynb).

The following video shows the  traditional Complex Ginzburg-Landau:

https://user-images.githubusercontent.com/14216783/224129416-69f13958-8244-4741-9cf3-d6b9e51144e2.mp4


## SCGL
An example of SCGL with additive noise is:

https://user-images.githubusercontent.com/14216783/224137213-3b0a5406-c470-48c6-b131-187e2313415a.mp4

An example of SCGL with multiplicative noise is:


https://user-images.githubusercontent.com/14216783/224166543-ab332e49-3d7c-4383-862e-034d7d6e3227.mp4


## Noise Examples


## Multidimension CGL

The multidimensional case is under development. The solution for the simplest case is:

1D:

<img src="https://raw.githubusercontent.com/rsautter/Noisy-Complex-Ginzburg-Landau/main/CGl1D.png" width=60% height=60%>

3D:

https://user-images.githubusercontent.com/14216783/224053529-af028b73-6f58-47f4-a0a4-0108392fa167.mp4

## Paper link
Under submission
