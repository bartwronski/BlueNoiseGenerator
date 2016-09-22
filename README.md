# BlueNoiseGenerator

This is my attempt of an implementation of a Siggraph 2016 paper "Blue-noise Dithered Sampling" by Iliyan Georgiev and Marcos Fajardo from Solid Angle.
You can find the link to the paper abstract here: https://www.solidangle.com/research/dither_abstract.pdf 

## Motivation 

Main motivation for using blue noise is avoiding "clumping" of noise and by using only high frequency of signal, getting perceptually better sample distribution. 
Blue noise can be used for dithering, offsetting samples, covering space (for example uniform disk sampling).

Mikkel Gjol and Mikken Svendsen from Playdead have shown some excellent examples of its use in their GDC presentation "Low complexity, high fidelity - rendering of INSIDE":
http://www.gdcvault.com/play/1023002/Low-Complexity-High-Fidelity-INSIDE 
http://f0716f2bff707a1b9e85-36c178e006d3d30c5b9c8dd905f8236a.r70.cf2.rackcdn.com/rendering_inside.pdf 

## Implementation 

Implementation is pretty straightforward implementation of the Siggraph 2016 paper with my attempts to make it "generalized" for more dimensions and larger vector values.
More than 2 dimenstions are useful for either volumetric effects, or treating time as next, 3rd or 4th dimension.
Vectors of variables can be used for things like vectors of many independent random variables (radius and angle in polar coordinates).
NOTE: I tested so far only 2 and 3 dimensions. :) 

It's implemented in standard C++ style, without any extravaganza.
Single file, compile it with clang or MSVS as a console app.

Also added a Mathematica helper visualization notebook https://github.com/bartwronski/BlueNoiseGenerator/blob/master/blue_noise.nb 

## Results

Given some 16 by 16 intial distribution:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/initial.png)

With following white-noise-like Fourier distribution:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/initial_fourier.png)

Algorithm generates following pattern:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/blue_noise.png)

With following Fourier distribution (less low frequencies):

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/blue_noise_fourier.png)

It works great under repetition / tiling!

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/blue_noise_repeated.png)


Here is an example of it working on a 2D vector:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/vector2.png)

Slice 1 of this vector:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/vector2a.png)

Slice 2 of this vector:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/vector2b.png)

Slightly more clumpy, but still pretty decent, here is Fourier of 1st slice and the whole image:

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/vector2_fourier.png)

![Image](https://github.com/bartwronski/BlueNoiseGenerator/blob/master/Images/vector2a_fourier.png)

## TODO

Support different sized of different dimensions sizes (for example for 3D noise having less depth slices).
