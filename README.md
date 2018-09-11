# Node-Parcellation
This is the documentation file for "Spatially-constrained exemplar-based node parcellation" of the human brain at the individual- and state-specific level.

The parcellation algorithm is implemented in C++, and is part of the open-source BioImage Suite Project. To run these files, you need to download BioImageSuite source code from:

After downloading BioImage Suite source code (under the name "bioimagesuite32_0b1_src"), follow the below steps:

1- Move the file 'vtkbisIndividualizeParcellation.cpp' to '/bioimagesuite32_0b1_src/Connectivity/'

2- Move the file 'bis_individualizeconnectivity.tcl' to '/bioimagesuite32_0b1_src/bioimagesuite/bis_algorithm/'

3- Build the package according to the BioImageSuite mannual (see http://bioimagesuite.yale.edu/manual/index.aspx)

4- Make sure you set the environment path, by running the following command:
source bioimagesuite32_0b1_src/build/setpaths.csh

5- Finally, the individualzied parcellation algorithm can be called using the following command:
bis_individualizeconnectivity.tcl -inp *INPUT NII FILE (VOXEL-LEVEL TIME SERIES)* -inp2 *INITIAL GROUP-LEVEL PARCELLATION*  -indiv_group *1 TO GENERATE INDIVIDUALIZED FILE* -blursigma *SMOOTHIN KERNEL BANDWIDTH* -num_exemplar *NUMBER OF NODES*


The number of nodes should match the group-level parcellation. Here for *Shen* parcellation, the number of nodes is 268 for whole brain analysis, and 188 for cortical analysis.

