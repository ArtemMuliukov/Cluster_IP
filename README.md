# Cluster_IP - version 0.3

## Content of the repository

### License
The libriary is a product of work under master thesys of Muliukov AR. You could use any content of the code in this repository under GPL-3.0 license.

### Content

Here you could find

## Presented files

Here you can find all significant parts of the work. Text version (and more deep description) of the work you can reed in the file [Dipoma_final_version](Diploma_final_version.pdf) (in Russian)

All current code is in the folder ["actual_code"](actual_code). There you can find 2 folders - "Labview" and "dll". 

["Labview"](actual_code/Labview) - with an example of dll usage in IDE LabVIEW. 
The code itself in the file ["Experiment_treatment1_ArtemVer1KG.vi"](actual_code/Labview/Experiment_treatment1_ArtemVer1KG.vi) and in  ["Experiment_treatment1_ArtemVer2_FastKG.vi"](actual_code/Labview/Experiment_treatment1_ArtemVer2_FastKG.vi) (for 2 functions of indicatrices analysis, with and without return of detailed clusters information).
How to connect the dll you can see at page 4 of code window in the LabVIEW .vi files.

["dll"](actual_code/dll) is a full project for Microsoft Visual Studio. You can reassemble it with any needed changes, all the code is in one file - ["cluster_ip_code.cpp"](NSU_Diplom_Clustering/actual_code/dll/test_dll_creation/cluster_ip_code.cpp). Also you may need to copy the file  ["cluster_ip_code.def"](NSU_Diplom_Clustering/actual_code/dll/test_dll_creation/cluster_ip_code.def) to correctly declare functions for you linker.

MVS project may be opened (with all needed dll installations) by runnig a file "test_dll_creation.sln", so by perpose the code is distributed with all needed config files.

Either you can find a prepared dll here: ["cluster_ip.dll"](NSU_Diplom_Clustering/actual_code/cluster_ip_dll/x64/Release/cluster_ip.dll).

Also you can find Jupyter notebooks used for tests in Python or constructing of graphs (for conferences or the diploma work itself) in folder ["Python tests"](Python tests). Short information about it's content you can find in local readme file.
