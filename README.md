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
**bis_individualizeconnectivity.tcl -inp** *Input* **-inp2** *Parc*  **-indiv_group** *1* **-blursigma** *BW* **-num_exemplar** *K*

- *Input* is the .nii file containing the voxel-level time series.

- *Parc* is the initial group-level parcellation that the algorihtm starts from. Please refer to the manuscript for more details.

- indiv_group option should be set to *1* in order to generate individualized parcellations, as well as the functional connectivity matrices using the individualize partcellations. If it is set to *0*, it will only output the functional connectivity matrices, using initial group-level parcellation.

- *BW* is the smoothing kernel's bandwidth. In the original work, we set *BW=4*.

- *K* is the number of nodes in the parcellation. *K* should match the group-level parcellation. Here for *Shen* parcellation, *K=268* for whole brain analysis, and *K=188* for cortical analysis.


