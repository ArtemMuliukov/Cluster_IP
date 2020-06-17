# Cluster_IP - version 0.3

## Content of the repository

### License
The libriary is a product of work under master thesys of Muliukov AR. You could use any content of the code in this repository under GPL-3.0 license.

### Content

This work is devoted to modernization of a standard method for solving inverse problems by building a hierarchical data structure over a pre-computed signal database. The proposed method was also tested to solve the problem of determining platelet parameters based on their signals measured on a scanning flow cytometer.

The classical method is based on modeling the signal using a third party software package (for example, for platelets, it can be [ADDA](https://github.com/adda-team/adda)) and comparing the measured signal with each of the available data in this database to find the nearest by L2 metric (euclidean). The method proposed here allows you to significantly speed up the search among the database objects by using KD-tree. A distinctive feature of the work is the output of discarded distances, which allows you to evaluate the solution over the entire space of possible solutions and use this knowledge to calculate statistical estimates of the resulting solution. The statictical estimations were presened in articles [1](https://www.spiedigitallibrary.org/journals/Journal-of-Biomedical-Optics/volume-19/issue-08/085004/Additivity-of-light-scattering-patterns-of-aggregated-biological-particles/10.1117/1.JBO.19.8.085004.full?SSO=1) and [2](https://www.spiedigitallibrary.org/journals/Journal-of-Biomedical-Optics/volume-14/issue-06/064036/Is-there-a-difference-between-T--and-B-lymphocyte/10.1117/1.3275471.full).

The libriary is close to classical realisation of kd-tree, except the possibility of flexible search with returning varing number of  seen clusters to enrich the solution with it's statistical estimation. Also you could vary the depth of search wich will help you to set convinient ratio speed/accuracy for the counted statistical estimations (see choose of ratio parameter in c_probs function).

The solution is designed as C++ code in Microsoft Visual Studio. The code is compiled into a single dll file that has 2 key functions. It's detailed declaration you an find in the code file[code file](actual_code/cluster_ip_dll/test_dll_creation/cluster_ip_code.cpp), here there are just short explications.

1) "build_tree" - allows you to load an array of data into the dll, which will be converted into a binary pseudo-tree and stored in the internal memory of the DLL for subsequent operation of the program. 

2) "c_probs"  - a function that finds the nearest element in the tree, as well as outputs the viewed and discarded classes and distances to them along the way, which can be used for subsequent calculations (an example can be found in the folder with tests in LabVIew)

The effectiveness of this implementation was tested on a real biological problem (you can find the tests in the ["Python_tests"](Python_tests) folder) and presented at various conferences(for ex. [here](https://www.giss.nasa.gov/staff/mmishchenko/ELS-XVII/)), as well as the work was presented as qualification master's thesis work.

## Presented files

Here you can find all significant parts of the work. Text version (and more deep description) of the work you can reed in the file [Dipoma_final_version](Extra_info-Publications/Diploma_final_version.pdf) (in Russian)

Less detaled but quite informative will be [this poster](Extra_info-Publications/poster_clustering_ELS18.pdf) (in English) from ELS18 (Texas, USA) or [this presentation](Extra_info-Publications/ISSC2020.pdf) (in Russian) from ISSC 2020 (Novosibirsk, Russia).

All current code is in the folder ["actual_code"](actual_code). There you can find 2 folders - ["LabView"](actual_code/LabView) and ["cluster_ip_dll"](actual_code/cluster_ip_dll). 

["LabView"](actual_code/LabView) - with an example of dll usage in IDE LabVIEW. 
The code itself in the file ["Experiment_treatment1_ArtemVer1KG.vi"](actual_code/LabView/Experiment_treatment1_ArtemVer1KG.vi) and in ["Experiment_treatment1_ArtemVer2_FastKG.vi"](actual_code/LabView/Experiment_treatment1_ArtemVer2_FastKG.vi) (for 2 functions of indicatrices analysis, with and without return of detailed clusters information).
How to connect the dll you can see at page 4 of code window in the LabVIEW .vi files.

["cluster_ip_dll"](actual_code/cluster_ip_dll) is a full project for Microsoft Visual Studio. You can reassemble it with any needed changes, all the code is in one file - ["cluster_ip_code.cpp"](actual_code/cluster_ip_dll/test_dll_creation/cluster_ip_code.cpp). Also you may need to copy the file  ["cluster_ip_code.def"](actual_code/cluster_ip_dll/test_dll_creation/cluster_ip_code.def) to correctly declare functions for you linker.

MVS project may be opened (with all needed dll installations) by runnig a file ["cluster_ip.sln"](actual_code/cluster_ip_dll/cluster_ip.sln), so by perpose the code is distributed with all needed config files. But if just want to use the dll, you can find it compitely prepared here: ["cluster_ip.dll"](actual_code/cluster_ip_dll/x64/Release/cluster_ip.dll).

Also you can find Jupyter notebooks used for tests in Python or constructing of graphs (for conferences or the diploma work itself) in folder ["Python_tests"](Python_tests). Short information about its content you can find in local readme file.
