# Codes_Microgeometry_PS
Matlab codes for microgrometry and RGB albedo estimation using photometric stereo without demosaicing.

## Introduction

These Matlab codes implement the method described in [1]. Given a set of non-demosaiced RAW photometric stereo images of a still scene acquired from a (pinhole or orthographic camera) camera, and dense lighting fields, this algorithm estimates depth, normals, and demosaiced RGB albedo. It can output a colored mesh in the .obj format. 

[1] "Microgeometry capture and RGB albedo estimation by photometric stereo without demosaicing", Yvain Quéau et al., Proceedings of the International Conference on Quality Control by Artificial Vision (QCAV 2017).

Please cite the above work if using the provided codes and/or datasets for your own research. 

Author: Yvain Quéau, Technical University Munich, yvain.queau@tum.de

Credits: Data acquisition was made possible by a technological transfer between the IRIT lab and the Pixience company, supervised by Toulouse Tech Transfer.

## Datasets

- The `Datasets/` folder contains three datasets: a 10euro banknote with microgeometry details hidden by a somewhat complicated albedo, and two Euro coins. The `Datasets/Calib/` folder contails the dense lighting fields associated with each image, and the camera's intrinsics. 

## Usage

The main fuction is `Toolbox/microgeometry.m` (see header for details).

Outputs:
- a gridded point cloud (the third dimension represents depth)
- normal map
- albedo

Inputs: 

- a data structure such that:
  * data.I contains the 3D or 4D image stack (**REQUIRED**)
  * data.mask contains a binary 2D mask
  
- a calib structure such that:
  * calib.S contains the dense lighting fields (**REQUIRED**)
  * calib.K contains the camera's intrinsics
  
- a param structure containing optional parameters:
  * params.z0 is the initial depth. We advise to roughly (visually) estimate the distance from camera to object in mm, and set z0 to a constant matrix with this rough distance 
  * params.ratio downsamples images by a factor of ratio
  
## Demo

The three demo files can be used to reproduce Figs. 2 and 4 in [1]. 
