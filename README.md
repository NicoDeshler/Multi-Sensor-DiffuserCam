# Multi-Sensor-DiffuserCam
The Multi-Sensor DiffuserCam project explores diffuser-based imaging systems with multiple sensors. This repository contains MATLAB scripts that interface with Zemax OpticStudio (here we are using version 17.5) allowing users to run analyses available in the OpticStudio suite on a lens file of their choice.


## Overview of Programmatic Interfacing with Zemax
Zemax offers two ways of programmatically interfacing with OpticStudio suite: *Interactive Applications* and *Standalone Applications* (These are described in greater detail below). Put simply, these applications serve as channels between Zemax and MATLAB, allowing users to operate OpticStudio via code rather than through the GUI. Consequently, users can quickly perform multiple iterations of a batch analysis on their optical designs and collect simulated data. Here we have made an interface with MATLAB available and offer a guide on how to use it.

### Interactive Applications
Interactive applications *actually open* a live session of OpticStudio when you run your MATLAB code. Thus, you should expect to see your analyses appearing in the OpticStudio suite as they are executed in the code. One primary advantage of having a live session of Zemax open when you run your MATLAB code is it offers a visual sanity check. You can verify that your code is calling the proper analyses and operating as intended on your lens file. The disadvantage is that the GUI slows down the speed of the analyses as Zemax renders images of the results. We've also found that the program will sometimes crash if the analysis and number of iterations is too computationally demanding. Hence, it is advisable to begin testing the functionality of your programmatic interface in an Interactive Application with few iterations. Once things are running smoothly, you can move your code over to a standalone application and run more demanding tasks with more iterations.   

### Standalone Applications
Standalone applications *do not open* a live sessin of OpticStudio when you run your MATLAB code. 


## Workflow

## The ZOS-API





