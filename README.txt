HSICS_TCI
Simulations for testing Multispectral Compressive Imaging Strategies
====================================================================

### Description:
This Matlab package comes as supplementary material for the paper 
K. Degraux, V. Cambareri, Bert Geelen, Laurent Jacques, Gauthier Lafruit,
Multispectral Compressive Imaging Strategies using Fabry-Pérot Filtered 
Sensors, Submitted to IEEE Transactions on Computational Imaging, 2018
It can be accessed and downloaded at https://github.com/kevd42/hsics_tci

### Size
Total size: 108.3 MB

### Packing list: 
code.zip 1.6MB
images.zip 4MB
matfiles.zip 102.7MB


### Instructions:
To generate the figures from the paper:
- Unzip code.zip
- Unzip matfiles.zip
- Unzip the toolboxes (spot-master and tight_subplot).
- Launch hsics_setup
- Launch hsics_psnr_figure
- Launch hsics_patches_figure

To reproduce, from the data, the results presented in the paper:
- Unzip images.7z (or images.zip) in the root folder.
- Download the CAVE dataset at
http://www.cs.columbia.edu/CAVE/databases/multispectral/
- Unzip it in the images folder.
- Unzip the toolboxes (spot-master and tight_subplot).
- Launch hsics_setup, then hsics_main(jobid) with jobid between 0 and 705.
For more information on the jobid  and the output files, see hsics_jobid.

The function HSinspector is useful to display and explore 3D arrays.

The code was tested on MATLAB R2016b on Mac OS X 10.10.5 and on 
MATLAB R2015a compiled with MCC and ran with MCR on a Linux workstation.

Copyright (C) 2014-2018 Kevin Degraux and Valerio Cambareri

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Corresponding author: kevin.degraux@uclouvain.be, kevin.degraux@gmail.com

Copyright (C) 2018 UCLouvain 
