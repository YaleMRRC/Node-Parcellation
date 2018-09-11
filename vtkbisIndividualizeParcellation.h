//BIOIMAGESUITE_LICENSE  ---------------------------------------------------------------------------------
//BIOIMAGESUITE_LICENSE  This file is part of the BioImage Suite Software Package.
//BIOIMAGESUITE_LICENSE  
//BIOIMAGESUITE_LICENSE  X. Papademetris, M. Jackowski, N. Rajeevan, H. Okuda, R.T. Constable, and L.H
//BIOIMAGESUITE_LICENSE  Staib. BioImage Suite: An integrated medical image analysis suite, Section
//BIOIMAGESUITE_LICENSE  of Bioimaging Sciences, Dept. of Diagnostic Radiology, Yale School of
//BIOIMAGESUITE_LICENSE  Medicine, http://www.bioimagesuite.org.
//BIOIMAGESUITE_LICENSE  
//BIOIMAGESUITE_LICENSE  This program is free software; you can redistribute it and/or
//BIOIMAGESUITE_LICENSE  modify it under the terms of the GNU General Public License version 2
//BIOIMAGESUITE_LICENSE  as published by the Free Software Foundation.
//BIOIMAGESUITE_LICENSE  
//BIOIMAGESUITE_LICENSE  This program is distributed in the hope that it will be useful,
//BIOIMAGESUITE_LICENSE  but WITHOUT ANY WARRANTY; without even the implied warranty of
//BIOIMAGESUITE_LICENSE  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//BIOIMAGESUITE_LICENSE  GNU General Public License for more details.
//BIOIMAGESUITE_LICENSE  
//BIOIMAGESUITE_LICENSE  You should have received a copy of the GNU General Public License
//BIOIMAGESUITE_LICENSE  along with this program; if not, write to the Free Software
//BIOIMAGESUITE_LICENSE  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//BIOIMAGESUITE_LICENSE  See also  http://www.gnu.org/licenses/gpl.html
//BIOIMAGESUITE_LICENSE  
//BIOIMAGESUITE_LICENSE  If this software is modified please retain this statement and add a notice
//BIOIMAGESUITE_LICENSE  that it had been modified (and by whom).  
//BIOIMAGESUITE_LICENSE 
//BIOIMAGESUITE_LICENSE  -----------------------------------------------------------------------------------


#ifndef __vtkbisIndividualizeParcellation_h
#define __vtkbisIndividualizeParcellation_h

#include "vtkSimpleImageToImageFilter.h"
#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkIdList.h"





#include "Eigen/Dense"
using Eigen::VectorXi;
using Eigen::MatrixXd;



class vtkbisIndividualizeParcellation : public vtkSimpleImageToImageFilter
{
public:
  vtkTypeMacro(vtkbisIndividualizeParcellation,vtkSimpleImageToImageFilter);
  static vtkbisIndividualizeParcellation *New();

  // Mask Image
  vtkSetObjectMacro(FMRIImage,vtkImageData);
  vtkGetObjectMacro(FMRIImage,vtkImageData);

  // Number of Exemplars
  vtkSetClampMacro(NumberOfExemplars,int,1,1000000);
  vtkGetMacro(NumberOfExemplars,int);

  // Do Hemisphere?
  vtkSetClampMacro(Hemisphere,int,0,1);
  vtkGetMacro(Hemisphere,int);

  // Do Individualzied?
  vtkSetClampMacro(IndivGroup,int,0,1);
  vtkGetMacro(IndivGroup,int);

  // Lambda
  vtkSetClampMacro(Lambda,float,0,1);
  vtkGetMacro(Lambda,float);

protected:

  vtkbisIndividualizeParcellation();
  virtual ~vtkbisIndividualizeParcellation();

  // Execute Function -- main function
  virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);
  virtual int  ComputeMRFIncrements(vtkImageData* img,int incr[]);
//  virtual int  ComputeMRFIncrements_Diag(vtkImageData* img,int incr[]);
  void ClusterAssignment(vtkImageData* input,vtkImageData* output, int* Voxel_indices, MatrixXd& distvSoptR, MatrixXd& distvSoptL);
  bool ismember(VectorXi V, int p);
private:
	vtkbisIndividualizeParcellation(const vtkbisIndividualizeParcellation& src) {};
	vtkbisIndividualizeParcellation& operator=(const vtkbisIndividualizeParcellation& rhs) {};
	int NumberOfExemplars;
	int Hemisphere;
	int IndivGroup;
	vtkImageData* FMRIImage;
	float Lambda;
};

#endif



