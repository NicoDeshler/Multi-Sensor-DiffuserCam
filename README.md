# Multi-Sensor DiffuserCam
The Multi-Sensor DiffuserCam project explores diffuser-based imaging systems with multiple sensors. This repository contains MATLAB scripts that interface with Zemax OpticStudio (version 13 and above) allowing users to run analyses available in the OpticStudio suite on a lens file of their choice.


## Overview of Programmatic Interfacing with Zemax
Zemax offers two ways of programmatically interfacing with OpticStudio suite: *Interactive Extensions* and *Standalone Applications* (These are described in greater detail below). Put simply, these modes serve as channels between Zemax and MATLAB, allowing users to operate OpticStudio via code rather than through the GUI. Consequently, users can quickly perform multiple iterations of a batch analysis on their optical designs and collect simulated data. Here we have written both an Interactive Application and a Standlone Application in MATLAB and offer a guide on how to use them.

### The Interactive Extension
The Interactive Extension *actually opens* a live session of OpticStudio when you run your MATLAB code. Thus, you should expect to see results from your analyses appearing in the OpticStudio suite as they are executed. One primary advantage of having a live session of Zemax open while your code runs is that it offers a visual sanity check. You can verify that your code is calling the proper analyses and operating as intended on your lens file. However, the GUI does slow down the speed of the analyses as Zemax automatically renders images of the results. We have also found that the program will sometimes crash if the number of iterations is too high. Hence, it is advisable to begin testing the functionality of your programmatic interface in an Interactive Extension with few iterations. Once things are running smoothly, you can move your code over to a Standalone Application and run demanding tasks with more iterations.   

### The Standalone Application
The Standalone Application *does not open* a live session of OpticStudio when you run your MATLAB code. Instead, the optical analyses available in Zemax are accessed independently from the GUI yet still run on your lens file. This improves performance and is the ideal way to programmatically interact with Zemax.   

## Writing your own Analysis and ZOS-API
Zemax uses a custom API known as the ZOS-API.NET. You will be using this API to call methods, run analyses, and get object attributes for the surfaces in your lens file. While we refrain from providing step-by-step instructions on setting up a MATLAB file that utilizes the API and interfaces with Zemax, we have listed a collection of informative articles and video tutorials from the Zemax Knowledgebase that should help you get started with this.

*Two vital resources for successfully implementing the ZOS-API in your code are the API documentation and the Intellisense feature.* 

You can find the API documentation by navigating through the following tabs in an open session of OpticStudio: `Help/ZOS-API/Documentation`. The documentation also contains some informative snippets of example code as a syntactic guide. Intellisense is an autocompletion feature that comes with the ZOS-API. As you write your code, Intellisense provides a drop down list of accessible methods and attributes for a given object instance. This feature will drastically streamline the process of writing functional code. 


### Articles
- [ZOS-API Overview](http://customers.zemax.com/os/resources/learn/knowledgebase/zos-api-net-an-overview)
- [Support Material for ZOS-API Users](http://customers.zemax.com/os/resources/learn/knowledgebase/support-material-for-zos-api-users)
- [Setting Up an Interactive Extension in MATLAB](http://customers.zemax.com/os/resources/learn/knowledgebase/how-to-connect-to-the-zos-api-with-the-interactive)
- [Sample MATLAB Code using ZOS-API: Some Common Analyses](http://customers.zemax.com/os/resources/learn/knowledgebase/zosapi-using-matlab)
- [Creating a User Analysis with ZOS-API](http://customers.zemax.com/os/resources/learn/knowledgebase/how-to-create-a-user-analysis-using-zos-api)

### Tutorials
- [General Interfacing and ZOS-API Overview](http://customers.zemax.com/os/opticstudio/opticstudio/user-interface/zos-api)
- [Setting Up a Standalone Application in MATLAB](http://customers.zemax.com/zmx/webinars/opticstudio-recordings/matlab-zos-api-net)

## Example 






