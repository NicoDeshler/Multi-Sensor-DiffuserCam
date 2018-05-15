# Multi-Sensor-DiffuserCam
The Multi-Sensor DiffuserCam project explores diffuser-based imaging systems with multiple sensors. This repository contains MATLAB scripts that interface with Zemax OpticStudio (version 13 and above) allowing users to run analyses available in the OpticStudio suite on a lens file of their choice.


## Overview of Programmatic Interfacing with Zemax
Zemax offers two ways of programmatically interfacing with OpticStudio suite: *Interactive Applications* and *Standalone Applications* (These are described in greater detail below). Put simply, these applications serve as channels between Zemax and MATLAB, allowing users to operate OpticStudio via code rather than through the GUI. Consequently, users can quickly perform multiple iterations of a batch analysis on their optical designs and collect simulated data. Here we have written both an Interactive Application and a Standlone Application in MATLAB and offer a guide on how to use them.

### The Interactive Application
The Interactive Application *actually opens* a live session of OpticStudio when you run your MATLAB code. Thus, you should expect to see results from your analyses appearing in the OpticStudio suite as they are executed. One primary advantage of having a live session of Zemax open while your code runs is that it offers a visual sanity check. You can verify that your code is calling the proper analyses and operating as intended on your lens file. However, the GUI does slow down the speed of the analyses as Zemax automatically renders images of the results. We've also found that the program will sometimes crash if the number of iterations is too high. Hence, it is advisable to begin testing the functionality of your programmatic interface in an Interactive Application with few iterations. Once things are running smoothly, you can move your code over to a standalone application and run more demanding tasks with more iterations.   

### The Standalone Application
The Standalone Application *does not open* a live session of OpticStudio when you run your MATLAB code. Instead, the optical analyses available in Zemax are accessed independently from the GUI and run on your lens file. This improves performance and is the ideal way to programmatically interact with Zemax.   

## Writing your own Analysis and The ZOS-API
### Setting Up Your Interactive Application - MATLAB FILE

### Setting Up Your Standalone Application - MATLAB File
See this [tutorial](http://customers.zemax.com/zmx/webinars/opticstudio-recordings/matlab-zos-api-net)
###

Zemax has its own API for accessing methods and object attributes. While we refrain from going into too much detail here, the 
-Intellisense
-Accessing documentation

## Setting Up



## Example 






