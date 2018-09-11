#!/bin/sh
# the next line restarts using wish \
    exec vtk "$0" -- "$@"

#BIOIMAGESUITE_LICENSE  ---------------------------------------------------------------------------------
#BIOIMAGESUITE_LICENSE  This file is part of the BioImage Suite Software Package.
#BIOIMAGESUITE_LICENSE  
#BIOIMAGESUITE_LICENSE  X. Papademetris, M. Jackowski, N. Rajeevan, H. Okuda, R.T. Constable, and L.H
#BIOIMAGESUITE_LICENSE  Staib. BioImage Suite: An integrated medical image analysis suite, Section
#BIOIMAGESUITE_LICENSE  of Bioimaging Sciences, Dept. of Diagnostic Radiology, Yale School of
#BIOIMAGESUITE_LICENSE  Medicine, http:#www.bioimagesuite.org.
#BIOIMAGESUITE_LICENSE  
#BIOIMAGESUITE_LICENSE  This program is free software; you can redistribute it and/or
#BIOIMAGESUITE_LICENSE  modify it under the terms of the GNU General Public License version 2
#BIOIMAGESUITE_LICENSE  as published by the Free Software Foundation.
#BIOIMAGESUITE_LICENSE  
#BIOIMAGESUITE_LICENSE  This program is distributed in the hope that it will be useful,
#BIOIMAGESUITE_LICENSE  but WITHOUT ANY WARRANTY; without even the implied warranty of
#BIOIMAGESUITE_LICENSE  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#BIOIMAGESUITE_LICENSE  GNU General Public License for more details.
#BIOIMAGESUITE_LICENSE  
#BIOIMAGESUITE_LICENSE  You should have received a copy of the GNU General Public License
#BIOIMAGESUITE_LICENSE  along with this program; if not, write to the Free Software
#BIOIMAGESUITE_LICENSE  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#BIOIMAGESUITE_LICENSE  See also  http:#www.gnu.org/licenses/gpl.html
#BIOIMAGESUITE_LICENSE  
#BIOIMAGESUITE_LICENSE  If this software is modified please retain this statement and add a notice
#BIOIMAGESUITE_LICENSE  that it had been modified (and by whom).  
#BIOIMAGESUITE_LICENSE 
#BIOIMAGESUITE_LICENSE  -----------------------------------------------------------------------------------

lappend auto_path [ file dirname [ info script ]]
lappend auto_path [file join [file join [ file dirname [ info script ]] ".." ] base]
lappend auto_path [file join [file join [ file dirname [ info script ]] ".." ] apps]

package require bis_dualimagetransformationalgorithm 1.0
package require bis_smoothimage 1.0
package provide bis_individualizeconnectivity 1.0
package require bis_matrixcorrelation 1.0
package require bis_roimean 1.0

itcl::class bis_individualizeconnectivity {

    inherit bis_dualimagetransformationalgorithm

     constructor { } {
	 $this Initialize
     }

    public method Initialize { }
    public method Execute { }

    public method GetGUIName { } { return "Individualize"}
}

# -----------------------------------------------------------------------------------------
# Initialize
# ----------------------------------------------------------------------------------------

itcl::body bis_individualizeconnectivity::Initialize { } {

    PrintDebug "bis_individualizeconnectivity::Initialize" 

    set options {
	{ num_exemplar   "Number of exemplars" "Start T"  { integer default   } 268 { 1   1000000 }  9 }	
	{ blursigma "kernel size [mm/voxel] of FWHM filter size" "Filter Size"  { real triplescale 100 } 6.0 { 0.0 20.0 }  0 }
	{ indiv_group "Individualized 1 or group 0 parcellation" "Do Indiv" boolean  1 {0 1} 0 }
	{ hemisphere "Whether to do each hemisphere separately 1 or the whole brain together 0" "do single hemi" boolean 1 {0 1} 0 }
	{ lambda "The relative importance of functional and spatial distances (dist = lambda*fxn + (1-lambda)*xyz)"  "lambda" { real triplescale 100 } 1 { 0 1 }  0 }
    }
	
    set defaultsuffix { "_individ" }   
    set scriptname bis_individualizeconnectivity

    set description "Calculate individualized parcellations."
    set description2 "Calculate individualized parcellations starting from a group-level parcellation."
    set backwardcompatibility ""
    set authors "mehraveh.salehi@yale.edu"

    set category "Functional Imaging"


    $this InitializeDualImageTransformationAlgorithm
}

# ----------------------------------------------------------------------------------------
# Execute
# ----------------------------------------------------------------------------------------

itcl::body bis_individualizeconnectivity::Execute {  } {
    
    set ok [ pxtclvtkpxcontrib::ConditionalLoadLibrary  vtkbisConnectivityTCL vtkbisROICorrelation 0  ]
    if { $ok == 0 } {
	set errormessage "Failed to load library vtkbisConnectivityTCL"
	return 0
    }

    set num_exemplar    [ $OptionsArray(num_exemplar) GetValue ]
    set image_in [ $InputsArray(input_image) GetObject ]
    set blursigma    [ $OptionsArray(blursigma) GetValue ]
    set hemisphere [ $OptionsArray(hemisphere) GetValue ]
    set lambda [ $OptionsArray(lambda) GetValue ]
    set indiv_group     [ $OptionsArray(indiv_group) GetValue ] 
    set indiv    [ $this GetOutput ] 
    set group    [ $this GetSecondInput ]


    set tmpname0 [ $InputsArray(input_image) GetFileName ]
#    puts stderr "extension = [file extension $tmpname0]"
    if { [ file extension $tmpname0 ] == ".gz" } { set tmpname0 [ file rootname $tmpname0 ] } 
    if { [ file extension $tmpname0 ] == ".nii" } { set tmpname0 [ file rootname $tmpname0 ] }
#    puts stderr "extension = [file extension $tmpname0]"
   
 
    set tmpname ${tmpname0}_WHOLEBRAIN
    if { $num_exemplar == 188 } { set tmpname ${tmpname0}_CORTICAL }
   
    if { $indiv_group ==  1 } { set ig "indiv" }
    if { $indiv_group ==  0 } { set ig "group" }

    if { $hemisphere ==  1 } { set hem "" }
    if { $hemisphere ==  0 } { set hem "WBtogether" }

############################## From here gets commented to generate group parcellation connectivities ##################################

    if { $indiv_group == 1 } {
	    # Step 1 ... Invoke bis_smoothimage -----
	    set smooth_alg [bis_smoothimage [pxvtable::vnewobj]]
	    $smooth_alg InitializeFromContainer $this
	    $smooth_alg SetInput $image_in 
	    $smooth_alg SetOptionValue blursigma $blursigma
	    $smooth_alg SetOptionValue unit "mm"
	    $smooth_alg Execute
	    
	    if { $blursigma > 0 } {
	    	set smimage [  [ $smooth_alg GetOutput ] GetImage ]
	    }
	    set individualize [ vtkbisIndividualizeParcellation [ pxvtable::vnewobj ]]
	    if { $blursigma > 0 } {
		$individualize SetFMRIImage $smimage
	    }
	    if { $blursigma == 0 } {
		puts stdout "No smoothing!"
		$individualize SetFMRIImage [$image_in GetImage]
	    }
	    $individualize SetInput [ $group GetObject ]
	    $individualize SetNumberOfExemplars $num_exemplar
	    $individualize SetHemisphere $hemisphere
	    $individualize SetIndivGroup $indiv_group
	    $individualize SetLambda $lambda
	    $individualize Update

	    $indiv ShallowCopyImage [ $individualize GetOutput ]
	    $indiv CopyImageHeader [ $image_in GetImageHeader ]

	    set comment [ format " [ $this GetCommandLine full ]" ]
	    [ $indiv GetImageHeader ] AddComment "$comment $Log" 0

	    $indiv SetFileName ${tmpname}_parcellation_${ig}_Priority_AllNorm1_lambda${lambda}.nii.gz
	    puts stderr "Indiv Description = [ $indiv GetDescription ]"

	    $individualize Delete
	    itcl::delete object $smooth_alg
    }
############################### Until here gets commented to generate group parcellation connectivities ##############################

    set roimean_alg  [ bis_roimean \#auto ]
    $roimean_alg InitializeFromContainer 0 $this
    $roimean_alg SetInput $image_in
    if { $indiv_group == 1 } { $roimean_alg SetSecondInput $indiv }
    if { $indiv_group == 0 } { $roimean_alg SetSecondInput $group }
    $roimean_alg SetOptionValue dotextfile 1
    $roimean_alg SetOptionValue filename ${tmpname}_roimean_Priority_AllNorm1_lambda${lambda}.txt
    $roimean_alg SetOptionValue filename2 ${tmpname}_xyzCoordinate_Priority_AllNorm1_lambda${lambda}.txt
    $roimean_alg Execute

#    if { $dosave == 1 } {
    set tmpname1 ${tmpname}_roimean_Priority_AllNorm1_lambda${lambda}.nii.gz
    [ $roimean_alg GetOutput ] Save $tmpname1
    puts stdout "ROI Mean saved in $tmpname1"
#    }


    set corr [ bis_matrixcorrelation \#auto ]
    $corr InitializeFromContainer 0 $this
    $corr SetInput [ $roimean_alg GetOutput ]
    $corr SetOptionValue dotextfile 1
    $corr SetOptionValue filename ${tmpname}_matrix_Priority_AllNorm1_lambda${lambda}.txt
    $corr SetOptionValue filename2 ${tmpname}_outputForConnectivityViewer_Priority_AllNorm1_lambda${lambda}.txt
    $corr Execute
    


    return 1
}


# -----------------------------------------------------------------------------------------
#  This checks if executable is called (in this case bis_individualizeconnectivity.tcl) if it is execute
# ----------------------------------------------------------------------------------------
 
if { [ file rootname $argv0 ] == [ file rootname [ info script ] ] } {
    # this is essentially the main function
 

    set alg [bis_individualizeconnectivity [pxvtable::vnewobj]]
    $alg MainFunction 
}

