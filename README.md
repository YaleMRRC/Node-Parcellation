# Node-Parcellation
This is the documentation file for "priority-based exemplar parcellation" of human brain. 

The node-level individualized parcellation is implemented in C++, and is part of the open-source BioImage Suite Project. To run these files, you need to download BioImageSuite source code from:

After downloading BioImage Suite source code, in the folder named "bioimagesuite32_0b1_src"
move the file 'vtkbisIndividualizeParcellation.cpp' to 

and '/bioimagesuite/bis_algorithm/bis_individualizeconnectivity.tcl'


How to?

1- Build the package according to the BioImageSuite mannual (see http://bioimagesuite.yale.edu/manual/index.aspx)
2- Make sure you set the path:
source bioimagesuite32_0b1_src/build/setpaths.csh

3- you can run the individualzied parcellation with the following command:

bis_individualizeconnectivity.tcl -inp individual_fmri_data.nii.gz -inp2 group_parcellation_scheme.hdr -num_exemplar K

where K is the number of exempalrs (or nodes) that you want to parcellate to. K should match the group-level parcellation. Here for the Shen parcellation K=268 for whole brain analysis and K=188 for the cortical regions.

The .m file calculate the homogeneity and DB indices for the group-level and individualized parcellations (see the paper).
For the .m file, note that you need to have NIFTI package in your MATLAB path to run the .m file.


1- First, download t
