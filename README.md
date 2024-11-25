# SlopeEXCL - Slope-based Exclusion

Creating a .tif-file excluding areas above a defined threshold. 

This code excludes zones, which are unsuitable for wind and solar energy placements. 
For wind turbines, a uniform threshold for all aspects is considered. 
For solar panels, a distinct threshold can be applied for south, and the north-east-west aspect. 

## Getting ready
This repository requires the digital elevation model in resolution 3s from Hydrosheds, available here: https://www.hydrosheds.org/hydrosheds-core-downloads

## Run it


## Repository structure

* `scripts`: contains the Python source code
* `notebooks`: contains the Python code as jupyter-notebook
* `data`: place for raw data
* `output`: will contain created .tif-files
