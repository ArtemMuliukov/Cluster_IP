# Cluster_IP

Files used in diploma work of Muliukov AR.

Here you can find all significant parts of the work. Text version (and more deep description) 
of the work you can reed in the file Dipoma_final_version.pdf

All current code is in the folder 'actual_code'. There you can find 2 folders - "Labview" and "dll". 

"Labview" - with an example of dll usage in IDE LabVIEW. 
The code itself in the file "Experiment_treatment1_ArtemVer1KG.vi" and in "Experiment_treatment1_ArtemVer2_FastKG.vi" (for 2 functions of
indicatrices analysis, with and without return of detailed clusters information).
You can find the connection of the dll in page 4 of code window.

"dll" is a full project for Microsoft Visual Studio. You can reassemble it with any needed changes, all needed code is in one file, 
by the adress "NSU_Diplom_Clustering/actual_code/dll/test_dll_creation/test.cpp". MVS project may be opened (with all needed dll installations) by runnig a file "test_dll_creation.sln". All needed files are supplied.

The prepeared dll you can already find on adress "NSU_Diplom_Clustering/actual_code/dll/Release/test_dll_creation.dll".


Also you can find Jupyter notebooks used for tests in Python or constructing of graphs in folder "Python tests".
Short information about it's usage you can find in local readme file.
