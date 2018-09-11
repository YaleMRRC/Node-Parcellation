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



#include "vtkbisIndividualizeParcellation.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "pxisinf.h"
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_map>
#include <functional>


#include "Eigen/Dense"
#include "Eigen/Sparse"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Map;
#include "igl/slice.h"
#include "igl/colon.h"
#include "igl/mat_min.h"
#include "igl/mat_max.h"
#include "igl/find.h"
using igl::slice;
using igl::colon;
using igl::mat_min;
using igl::mat_max;
using igl::find;


// This creates the "New" Function
vtkStandardNewMacro(vtkbisIndividualizeParcellation);

// Construct object to extract all of the input data.
vtkbisIndividualizeParcellation::vtkbisIndividualizeParcellation()
{
  this->NumberOfExemplars=268;
  this->FMRIImage=NULL;
  this->Hemisphere=1;
  this->Lambda=1;
}
vtkbisIndividualizeParcellation::~vtkbisIndividualizeParcellation()
{
  this->FMRIImage=NULL;
}

// ---------------------------------------------------------------------------
//input=group level parcellation
//output=individualize parcellation

bool isAnyFalse(std::vector<bool> &v)
{
  for (int i=0; i<v.size(); i++)
	if (v[i] == false)
		return true;
  return false;
}


int  vtkbisIndividualizeParcellation::ComputeMRFIncrements(vtkImageData* img,int incr[6])
{
  int dim[3]; img->GetDimensions(dim);
  int slicesize=dim[0]*dim[1];
  int index=0;
  double d[3];

  for (int ic=-1;ic<=1;ic++)
    {
      for (int ib=-1;ib<=1;ib++)
        {
          for (int ia=-1;ia<=1;ia++)
            {
              if ((ia+ib+ic==1 || ia+ib+ic==-1) && (ia==0 || ib==0 || ic==0))
                {
//		  fprintf(stderr,"(%d, %d, %d)\n", ia, ib, ic);
                  incr[index]=ia+ib*dim[0]+ic*slicesize;
                  ++index;
                }
            }
        }
    }


  return 1;
}

bool vtkbisIndividualizeParcellation::ismember(VectorXi V, int p)
{
	for (int i=0;i<V.size();i++)
	{
		if (V(i)==p)
			return 1;
	}
	return 0;
}

void vtkbisIndividualizeParcellation::SimpleExecute(vtkImageData* input, vtkImageData* output)
{


  const clock_t begin_time = clock();
  fprintf(stderr,"Start!\n");
  int dim[3], dim2[3];
  if (this->FMRIImage!=NULL)
    {
      this->FMRIImage->GetDimensions(dim);
      input->GetDimensions(dim2);
	  fprintf(stderr,"fMRI Image dim = %dx%dx%d\n",dim[0],dim[1],dim[2]);
	  fprintf(stderr,"input Image dim = %dx%dx%d\n",dim2[0],dim2[1],dim2[2]);
      int sum=0;
      for (int ia=0;ia<=2;ia++)
      	sum+=abs(dim[ia]-dim2[ia]);
      if (sum>0)
	{
	  fprintf(stderr,"Bad FMRI Input to vtkbisIndividualizeParcellation SimpleExecute - sum = %d\n",sum);
	  return;
	}
    } 
  else 
    {
      fprintf(stderr,"Bad FMRI Input to vtkbisIndividualizeParcellation SimpleExecute - Null image\n");
      return;
    }
  
   
  vtkDataArray* group=input->GetPointData()->GetScalars(); 
  int N=group->GetNumberOfTuples(); //nt

  vtkDataArray* indiv=output->GetPointData()->GetScalars();
  indiv->FillComponent(0,0.0);

  vtkDataArray* fmri=this->FMRIImage->GetPointData()->GetScalars();
  const int t=this->FMRIImage->GetNumberOfScalarComponents(); //nc

  double range[2]; 
  group->GetRange(range);
  
  int Pmax = this->NumberOfExemplars; // Whole brain = 268, Cortical = 188
  double lambda = this->Lambda;

  if (Pmax != range[1]){
	fprintf(stderr,"Bad Group Parcellation Input to vtkbisIndividualizeParcellation SimpleExecute - pmax = %d, range[1] = %d\n",Pmax,range[1]);
	return;	
  }

  fprintf(stderr,"number of frames is %d\n",t);
  fprintf(stderr,"number of voxels is %d x %d x %d = %d\n",dim[0],dim[1],dim[2],N);
  fprintf(stderr,"number of exemplars is %d\n",Pmax);
  fprintf(stderr,"lambda is %f\n",lambda);

  fprintf(stderr,"START: the elapsed time is = %f s \n",float(clock()-begin_time)/CLOCKS_PER_SEC);

  

// Copying data to matrixXd and VectorXd

  const clock_t timebegin1 = clock();
  int count=0;

  for (int voxel=0;voxel<N;voxel++)
	if(group->GetComponent(voxel,0)>0)
		count++;

  MatrixXd X(t,count);
  count=0;
  double position[3];

  for (int voxel=0;voxel<N;voxel++)
    {

	if(group->GetComponent(voxel,0)>0)
	  {
		input->GetPoint(voxel,position);
		for (int frame=0;frame<t;frame++)
			X(frame,count) = fmri->GetComponent(voxel,frame);
		count++;
	  }
    }


  
  VectorXd parcel(count);
  std::unordered_map<int,int> ntoNvoxel;
  std::unordered_map<int,int> Ntonvoxel;
  count=0;
  for (int voxel=0;voxel<N;voxel++)
	if (group->GetComponent(voxel,0)>0)
	  {
		parcel(count) = group->GetComponent(voxel,0); // this is the group label for all the nonzero voxels
		ntoNvoxel.insert( std::make_pair<int,int>(count,voxel) );
		Ntonvoxel.insert( std::make_pair<int,int>(voxel,count) );
		count++;
	  }

  int n = count;

  fprintf(stderr,"COPYING DATA TO EIGEN MATRIX & VECTOR: the elapsed time is = %f s \n",float(clock()-timebegin1)/CLOCKS_PER_SEC);
  fprintf(stderr,"number of non-zero voxels n = %d\n",n);



  const clock_t timebegin11 = clock();
  
  VectorXd mean_subtract(t);
  mean_subtract = (X.rowwise().sum())/double(n);


// Normalizing data points to 0 mean
  double newValue;
  MatrixXd V = MatrixXd(t,n);
  V = X.colwise()-mean_subtract; // mean of V is all 0 [VERIFIED]

  

// Calculating the l2-norm
  VectorXd twoNorm(n);
  twoNorm =V.colwise().norm();
  VectorXd inverse_twoNorm(n);
  inverse_twoNorm = twoNorm.array().inverse();

////////////////////////////////////////////////  1- Dividing with the maximum norm 
//// finding the maximum value in the array
//  double maxValue = twoNorm.maxCoeff();
////  Normalizing data points to a unit ball sphere
//  MatrixXd v = V/maxValue;
//  fprintf(stderr,"NORMALIZATION INTO UNIT SPHERE (DIVIDE BY MAX NORM): the elapsed time is = %f s \n",float(clock()-timebegin11)/CLOCKS_PER_SEC);

////////////////////////////////////////////////  2- Normalizing to the unit norm (all vectors norm = 1)
  MatrixXd v = V.array().rowwise()* inverse_twoNorm.transpose().array();
  twoNorm = v.colwise().norm(); // twoNorm is all 1 [VERIFIED]

//  VectorXd test_ii(t);
//  test_ii = (V.rowwise().sum())/double(n);
//  for (int ii=0; ii<t; ii++)
//	  fprintf(stderr,"V norm = %f\n",twoNorm(ii)); 

  fprintf(stderr,"NORMALIZATION ONTO UNIT SPHERE (ALL NORM=1): the elapsed time is = %f s \n",float(clock()-timebegin11)/CLOCKS_PER_SEC);


////////// Finding the voxels within each parcel ///////////////

// Finding the voxels within each parcel
  const clock_t timebegin2 = clock();
  std::vector< std::vector<int> > indice_p;
  std::vector<int> p_vector;
  for (int p=0;p<Pmax;p++)
    {	
	for (int voxel=0;voxel<n;voxel++)
		if (parcel(voxel) == p+1)
			p_vector.push_back(voxel);
	indice_p.push_back(p_vector);
	p_vector.clear(); // p_vector's size is 0 [VERIFIED]
    }

  fprintf(stderr,"FINDING VOXELS: the elapsed time is = %f s \n",float(clock()-timebegin2)/CLOCKS_PER_SEC); 






// Calculating the squared distance matrix between voxels within each parcel
  const clock_t timebegin3 = clock();
  std::vector<MatrixXd> sqrDist;
  MatrixXd D;
  VectorXi R(t);

  colon(0,1,t-1,R);


  for (int p=0;p<Pmax;p++)
    {
	int psize = indice_p[p].size();
	MatrixXd sqrMatrix(psize,psize);
	int* ptr = &indice_p[p][0];
	Map<VectorXi> C(ptr,psize);
 	MatrixXd vP(t,psize);
	slice(v,R,C,vP);		
	sqrMatrix = ((vP.transpose()*vP*-2).colwise() + vP.colwise().squaredNorm().transpose()).rowwise() + vP.colwise().squaredNorm();


        sqrDist.push_back( sqrMatrix );
    }

fprintf(stderr,"SQUARED DISTANCES: the elapsed time is = %f s \n",float(clock()-timebegin3)/CLOCKS_PER_SEC);


// Calculating the distance between auxiliary exemplar and the rest of the voxels
  const clock_t timebegin4 = clock();

  std::vector<VectorXd> e0sqrDist;
  VectorXd e0 = VectorXd::Zero(t);
  e0(0) = 3;

 

  for (int p=0;p<Pmax;p++)
    {
        int psize = indice_p[p].size();
	VectorXd sqrArray(psize);
	int* ptr = &indice_p[p][0];
	Map<VectorXi> C(ptr,psize);
	MatrixXd vP(t,psize);
   	slice(v,R,C,vP);		
	sqrArray = ((vP.transpose()*e0*-2).colwise() + vP.colwise().squaredNorm().transpose()).rowwise() + e0.colwise().squaredNorm();	


        e0sqrDist.push_back( sqrArray );


//	delete [] sqrArray;
   }
//	for (int pp1=0;pp1<indice_p[p].size();pp1++)
//		for (int pp2=0;pp2<indice_p[p].size();pp2++)
//			fprintf(stderr,"p=%d, pp1=%d, pp2=%d, sqrDist[p][pp1][pp2] = %f\n",p,pp1,pp2,sqrDist[p][pp1][pp2]);
  
  
fprintf(stderr,"AUXILIARY DISTANCES: the elapsed time is = %f s \n",float(clock()-timebegin4)/CLOCKS_PER_SEC);

//  Calculating the exemplar within each parcel
  const clock_t timebegin5 = clock();
  VectorXi Sopt(Pmax);
  VectorXi SoptN(Pmax);
  double loss;
  for (int p=0;p<Pmax;p++)
    {
	int psize = indice_p[p].size();
	double sumd0 = e0sqrDist[p].sum();
	MatrixXd::Index maxFindex;
	VectorXd sumD(psize);	
	sumD = sqrDist[p].colwise().sum();

	VectorXd pFunc(psize);
	pFunc = sumd0 - sumD.array();  // we should divide by n but does not matter as it does not change the maximum!
	double maxF = pFunc.maxCoeff(&maxFindex);
	Sopt(p) = indice_p[p][maxFindex];
	std::unordered_map<int,int>::const_iterator voxelN = ntoNvoxel.find (Sopt(p));
	if (voxelN != ntoNvoxel.end())
	   SoptN(p) = voxelN->second;
    }

fprintf(stderr,"EXEMPLAR IDENTIFICATION: the elapsed time is = %f s \n",float(clock()-timebegin5)/CLOCKS_PER_SEC);

//  for (int p=0; p<Pmax; p++)
//	fprintf(stderr,"Sopt(%f) = %d, ", group->GetComponent(SoptN(p),0),SoptN(p));
 
// Assigning each voxel to the closest exemplar using the priority queue algorithm

  	

  


  const clock_t timebegin6 = clock();
  std::vector<double> label(n,-1);

// minDistIndexR and minDistIndexL need to be double, otherwise the result differs from MATLAB. It is perhaps because of the SetComponent() command.

  MatrixXd vSopt(t,Pmax);
  slice(v,R,Sopt,vSopt);

  MatrixXd distvSopt(n,Pmax);

  distvSopt = ((v.transpose()*vSopt*-2).colwise() + v.colwise().squaredNorm().transpose()).rowwise() + vSopt.colwise().squaredNorm();	

  dim[3]; input->GetDimensions(dim);
  const unsigned int neighbors = 6;
  int incr[neighbors]; this->ComputeMRFIncrements(input,incr);

  int sumVisited = 0;
  std::vector<int> VISITED(n,0);


// Assigning labels to exemplars and marking them as visited
  for (int p=0; p<Pmax; p++){
	label[Sopt(p)] = p;
	VISITED[Sopt(p)] = 1;
	sumVisited ++;		
  }
  
  std::vector<std::priority_queue<std::pair<double, int> , std::vector<std::pair<double,int> >, std::greater<std::pair<double, int> > > > exemplar_min_heaps (Pmax);
  for (int p=0; p<Pmax; p++)
  {
	int exemplarN = SoptN(p);
//	if (NONBOUNDARY[exemplarN] == 1)
   	   for (int ia=0;ia<neighbors;ia++)
	   {	
		int currVoxN = exemplarN+incr[ia]; 
		if (group->GetComponent(currVoxN,0)>0)
		{
			std::unordered_map<int,int>::const_iterator voxeln = Ntonvoxel.find (currVoxN);
			if (voxeln != Ntonvoxel.end())
			{
			   int currVox = voxeln->second; 
			   exemplar_min_heaps[p].push(std::make_pair(distvSopt(currVox,p),currVox));
			}
		}
	   }  
  }
  
 
  while (sumVisited < n)
  {
    int min_idx = -1;
    double min_val = 10000000;
    for (int p=0; p<Pmax; p++)
    {
	if (!exemplar_min_heaps[p].empty())
	{
		std::pair<double,int> curNode = exemplar_min_heaps[p].top();		
		if (curNode.first < min_val)
		{
			min_val = curNode.first; // all distances <=4 [VERIFIED]
		  	min_idx = p;
		}
	}	
	
    }
    // selected exemplar queue (min_idx) is correct [VERIFIED]
    // fprintf(stderr,"\n(selected queue,min_val)=(%d,%f)\n",min_idx,min_val);
  
    if (min_idx >=0)
    {
      std::pair<double,int> chosenNode = exemplar_min_heaps[min_idx].top();
      exemplar_min_heaps[min_idx].pop();
      int chosenVoxel = chosenNode.second;
      if (VISITED[chosenVoxel] == 0)
      {
        label[chosenVoxel] = min_idx;
        VISITED[chosenVoxel] = 1;
        sumVisited ++;
//      fprintf(stderr,"number of visited = %d\n", sumVisited);
        std::unordered_map<int,int>::const_iterator voxelN = ntoNvoxel.find (chosenVoxel);
        if (voxelN != ntoNvoxel.end())
        {
	  int chosenVoxelN = voxelN->second; 
//      if(NONBOUNDARY[chosenVoxelN] == 1)
	  for (int ia=0;ia<neighbors;ia++)
	  {
		int currVoxN = chosenVoxelN + incr[ia];
		if (group->GetComponent(currVoxN,0)>0)
		{
			std::unordered_map<int,int>::const_iterator voxeln = Ntonvoxel.find (currVoxN);
			if (voxeln != Ntonvoxel.end())
			{
				int currVox = voxeln->second;
				if (VISITED[currVox] == 0){
					exemplar_min_heaps[min_idx].push(std::make_pair(distvSopt(currVox,min_idx),currVox));
				}
			}
		}
	  }
        }
      }

    } 
	/* testing the order of elements in the priority queue is from smallest dist to largest dist [VERIFIED]
	std::queue<std::pair<double, int> > tteesstt;	
	for (int ii=0;ii<exemplar_min_heaps[min_idx].size();ii++)
	{
		std::pair<double,int> testNode = exemplar_min_heaps[min_idx].top();
		tteesstt.push(testNode);
		exemplar_min_heaps[min_idx].pop();
		fprintf(stderr,"(%d,%f), ",testNode.second,testNode.first);
	}
	fprintf(stderr,"next...\n");
	*/

    if (min_idx < 0)
	break;
  }

  fprintf(stderr,"ASSIGNING VOXELS TO EXEMPLARS: the elapsed time is = %f s \n",float(clock()-timebegin6)/CLOCKS_PER_SEC);
  
  const clock_t timebegin7 = clock();
  int Voxel_indices[N];
  count = 0;

  for (int voxel=0;voxel<N;voxel++)
    {
	if(group->GetComponent(voxel,0)>0)
	  {
		std::unordered_map<int,int>::const_iterator voxeln = Ntonvoxel.find (voxel);
		if (voxeln != Ntonvoxel.end())
        	{
	        	int finalvoxeln = voxeln->second; 
			indiv->SetComponent(voxel,0,label[finalvoxeln]+1);
		}
	  }
	else
	   {
		indiv->SetComponent(voxel,0,0);	
	   }

    }
  fprintf(stderr,"WRITING IMAGE: the elapsed time is = %f s \n",float(clock()-timebegin7)/CLOCKS_PER_SEC);


}
