# Cluster_IP - version 0.3

## Content of the repository

### License
The library is a product of a MSc project of Muliukov A.R. You can use any content of the code in this repository under GPL-3.0 license.

### Content

This work is devoted to modernization of a standard method for solving parametic inverse problems by building a hierarchical data structure over a pre-computed signal database. The proposed method was tested to solve the problem of determining platelet parameters based on their signals measured on a scanning flow cytometer.

The classic method is based on simulating the signal for particles with a given shape model (e.g., spheroid or biconcave spheroid) determined by several parameters, using a third-party software package (e.g., for platelets, it can be [ADDA](https://github.com/adda-team/adda)). A large database is build once by varying the particle parameters in predefined ranges. Then, any measured signal is compared against all database elements (simulated signals) to find the nearest one in terms of L2 (Euclidean) metric. The parameters of this simulated signals are assigned to the measured one, i.e. the nearest-neighbor interpolation is used ([Moskalensky et al., _J. Biomed. Opt._, 2013](http://doi.org/10.1117/1.JBO.18.1.017001)) 

The method implemented in this library allows one to significantly speed up the search over the database by using KD-tree. A distinctive feature of the work is the additional output of distances to all considered database elements or clusters of them, which allows one to get more information about the solution by calcualting statistical estimates of the parameters (such as mathematical expectation and standard deviation), see ([Strokotov et al., _J. Biomed. Opt._, 2009](http://doi.org/10.1117/1.3275471)) for details.

The library is close to classic realisation of kd-tree, with addition of the ability to return the distances to considered clusters. Also one could vary the depth of search to set appropriate compromise between speed and accuracy of the statistical estimates (see the choice of ratio parameter in function `c_probs`).

The solution is designed as C++ code in Microsoft Visual Studio. The code is compiled into a single DLL file that has two key functions. Their detailed declarations can be found in the [`cluster_ip.cpp`](actual_code/cluster_ip_dll/cluster_ip/cluster_ip.cpp), here there are only short explanations:

1) `build_tree` - allows you to load an array of data into the DLL, which will be converted into a binary pseudo-tree and stored in the internal memory of the DLL for subsequent use. 

2) `c_probs`  - a function that finds the nearest element in the tree, as well as outputs the considered and discarded clusters (and distances to them) along the way. And example of using the latter distance for further calculations is given as LabView program in the folder [`actual_code/LabView`](actual_code/LabView).

The effectiveness of this implementation was tested on a real biological problem (you can find the tests in the folder [`Python_tests`](Python_tests)). The corresponding results were presented at various conferences (e.g., at [ELS-17](https://www.giss.nasa.gov/staff/mmishchenko/ELS-XVII/)), and as part of the MSc thesis of Muliukov A.R.

## Presented files

Here you can find all significant parts of the work. Text version (and more deep description) of the work is given in the [thesis](Extra_info-Publications/Diploma_final_version.pdf) (in Russian). Less detailed but quite informative is [this poster](Extra_info-Publications/poster_clustering_ELS17.pdf) (in English) from [ELS-17](https://www.giss.nasa.gov/staff/mmishchenko/ELS-XVII/) (Texas, USA) or [this presentation](Extra_info-Publications/ISSC2020.pdf) (in Russian) from ISSC 2020 (Novosibirsk, Russia).

All current code is in the folder [`actual_code`](actual_code). There you can find two folders - [`LabView`](actual_code/LabView) and [`cluster_ip_dll`](actual_code/cluster_ip_dll). 

[`LabView`](actual_code/LabView) - with an example of DLL usage in IDE LabVIEW. 
The code itself in the file [`Experiment_treatment1_ArtemVer1KG.vi`](actual_code/LabView/Experiment_treatment1_ArtemVer1KG.vi) and in [`Experiment_treatment1_ArtemVer2_FastKG.vi`](actual_code/LabView/Experiment_treatment1_ArtemVer2_FastKG.vi) (for two functions for processing of experimental light-scattering patterns, with and without return of detailed clusters information). How to connect the DLL, one can see at page 4 of the code window in the LabVIEW .vi files.

[`cluster_ip_dll`](actual_code/cluster_ip_dll) is a complete project for Microsoft Visual Studio. You can reassemble it with any needed changes, the code itself is in a single file - [`cluster_ip.cpp`](actual_code/cluster_ip_dll/cluster_ip/cluster_ip.cpp). Also you may need to copy the file  [`cluster_ip.def`](actual_code/cluster_ip_dll/cluster_ip/cluster_ip.def) to correctly declare functions for you linker.

MVS project may be opened (with all needed DLL installations) by running a file [`cluster_ip.sln`](actual_code/cluster_ip_dll/cluster_ip.sln), so the code is intentionally distributed with all needed config files. But you can also use the DLL out of the box, see [`cluster_ip.dll`](actual_code/cluster_ip_dll/x64/Release/cluster_ip.dll).

Also you can find Jupyter notebooks used for tests in Python or for constructing figures (for conferences and the thesis) in the folder [`Python_tests`](Python_tests). Short information about its content you can find in the local [readme file](Python_tests/readme.txt).
