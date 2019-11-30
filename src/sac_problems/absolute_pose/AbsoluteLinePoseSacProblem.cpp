/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/


#include <opengv/sac_problems/absolute_pose/AbsoluteLinePoseSacProblem.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <iostream>

bool
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::computeModelCoefficients(
    const std::vector< int >& indices, model_t& outModel) const
{
  transformations_t solutions;

  return false;
}

bool
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::compute2PModelCoefficients(
    const std::vector< int >& indices, model_t& outModel1, model_t& outModel2)
{
  transformations_t solutions;
  
  switch(_algorithm)
  {
  
  case TWOEMixed:
  {
	solutions = opengv::absolute_pose::twoe2p(_adapter,_mpitch,_mroll,indices);
	break; 
  }
  case MCMixed:
  {
	  solutions = opengv::absolute_pose::mc2p(_adapter,_mpitch,_mroll,_myaw,indices);
	  break;
	}
	  
  }
  
  if( solutions.size() <= 0 )
  {
    return false;
  }
  
  if (solutions.size() == 2)
  {
    outModel1 = solutions[0];
    outModel2 = solutions[1];
    return true;
  }
  
  return false;

}


bool
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::compute1P1LModelCoefficients(
    const std::vector< int >& indices, model_t& outModel1, model_t& outModel2)
{
  transformations_t solutions;
  
  switch(_algorithm)
  {
    
  case TWOE1P1L:
  {
	solutions = opengv::absolute_pose::twoe1p1l(_adapter,_mpitch,_mroll,indices);
	break; 
  }
  case TWOEMixed:
  {
	solutions = opengv::absolute_pose::twoe1p1l(_adapter,_mpitch,_mroll,indices);
	break; 
  }
  case MC1P1L:
  {
	  solutions = opengv::absolute_pose::mc1p1l(_adapter,_mpitch,_mroll,_myaw,indices);
	  break; 
  }
  case MCMixed:
  {
	  solutions = opengv::absolute_pose::mc1p1l(_adapter,_mpitch,_mroll,_myaw,indices);
	  break;
	}
	  
  }
  
  if( solutions.size() <= 0 )
  {
    return false;
  }
  
  if (solutions.size() == 2)
  {
    outModel1 = solutions[0];
    outModel2 = solutions[1];
    return true;
  }
  
  return false;

}


void
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::getSelectedDistancesToModel(
    const model_t & model,
    const std::vector<int> & indices,
    std::vector<double> & scores) const
{
  //compute the reprojection error of all points

  //compute inverse transformation
  model_t inverseSolution;
  inverseSolution.block<3,3>(0,0) = model.block<3,3>(0,0).transpose();
  inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*model.col(3);

  Eigen::Matrix<double,4,1> p_hom;
  p_hom[3] = 1.0;

  for(size_t i = 0; i < indices.size(); i++)
  {
    //get point in homogeneous form
    p_hom.block<3,1>(0,0) = _adapter.getPoint(indices[i]);

    //compute the reprojection (this is working for both central and
    //non-central case)
    point_t bodyReprojection = inverseSolution * p_hom;
    point_t reprojection =
        _adapter.getCamRotation(indices[i]).transpose() *
        (bodyReprojection - _adapter.getCamOffset(indices[i]));
    reprojection = reprojection / reprojection.norm();

    //compute the score
    scores.push_back(
        1.0 - (reprojection.transpose() * _adapter.getBearingVector(indices[i])));
  }
}

void
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::getSelectedDistancesToModel(
      const model_t & model,
      const std::vector<int> & indices,
      const std::vector<int> & line_indices,
      std::vector<double> & scores,
	  std::vector<double> & line_scores) const
{
	//compute the reprojection error of all points

  //compute inverse transformation
  model_t inverseSolution;
  inverseSolution.block<3,3>(0,0) = model.block<3,3>(0,0).transpose();
  inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*model.col(3);

  Eigen::Matrix<double,4,1> p_hom;
  p_hom[3] = 1.0;

  for(size_t i = 0; i < indices.size(); i++)
  {
    //get point in homogeneous form
    p_hom.block<3,1>(0,0) = _adapter.getPoint(indices[i]);

    //compute the reprojection (this is working for both central and
    //non-central case)
    point_t bodyReprojection = inverseSolution * p_hom;
    point_t reprojection =
        _adapter.getCamRotation(indices[i]).transpose() *
        (bodyReprojection - _adapter.getCamOffset(indices[i]));
    reprojection = reprojection / reprojection.norm();

    //compute the score
    scores.push_back(
        1.0 - (reprojection.transpose() * _adapter.getBearingVector(indices[i])));
  }
  
  Eigen::Matrix<double,4,1> sp_hom;
  sp_hom[3] = 1.0;
  Eigen::Matrix<double,4,1> ep_hom;
  ep_hom[3] = 1.0;
  for(size_t i=0; i<line_indices.size();++i)
  {
	  sp_hom.block<3,1>(0,0)=_adapter.getLinestartPoint(line_indices[i]);
	  ep_hom.block<3,1>(0,0)=_adapter.getLineendPoint(line_indices[i]);
	  point_t body_ptstart = inverseSolution*sp_hom;
	  point_t body_ptend = inverseSolution*ep_hom;
	  point_t ptstart = _adapter.getCamRotation(line_indices[i]).transpose() * 
	  (body_ptstart - _adapter.getCamOffset(line_indices[i]));
	  point_t ptend = _adapter.getCamRotation(line_indices[i]).transpose() *
        (body_ptend - _adapter.getCamOffset(line_indices[i]));
	  ptstart = ptstart/ptstart[2];
	  ptend = ptend/ptend[2];
	  double mx=ptend[1]-ptstart[1];
	  double my=ptstart[0]-ptend[0];
	  double mz=-mx*ptstart[0]-my*ptstart[1];
	  if (mx*mx+my*my<0.0001)
	  {
		  line_scores.push_back(10000);
		  continue;
	}
	  Eigen::Vector3d m(mx,my,mz);
	  point_t image_ptstart = _adapter.getLinestartBearingVector(line_indices[i]);
	  image_ptstart = image_ptstart/image_ptstart[2];
	  point_t image_ptend = _adapter.getLineendBearingVector(line_indices[i]);
	  image_ptend = image_ptend/image_ptend[2];
	  Eigen::Matrix2Xd A;
	  A.block<1,3>(0,0)=image_ptstart.transpose();
	  A.block<1,3>(1,0)=image_ptend.transpose();
	  Eigen::Matrix2d B;
	  B << 1,0.5,0.5,1;
	  B=B/(3*(mx*mx+my*my));
	  line_scores.push_back(m.transpose()*(A.transpose()*B*A)*m);
}
		
}


void
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::optimizeModelCoefficients(
    const std::vector<int> & inliers,
    const model_t & model,
    model_t & optimized_model)
{
  _adapter.sett(model.col(3));
  _adapter.setR(model.block<3,3>(0,0));
  optimized_model = opengv::absolute_pose::optimize_nonlinear(_adapter,inliers);
}

void
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::optimizeModelCoefficients(
    const std::vector<int> & inliers,
    const std::vector<int> & line_inliers,
    const model_t & model,
    model_t & optimized_model)
{
  _adapter.sett(model.col(3));
  _adapter.setR(model.block<3,3>(0,0));
  optimized_model = opengv::absolute_pose::optimize_nonlinear(_adapter,inliers); // only use points
//   optimized_model = opengv::absolute_pose::optimize_nonlinear(_adapter,inliers,line_inliers);
}

int
opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::getSampleSize() const
{
  int sampleSize = 3;
  return sampleSize;
}


void opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::get1P1LSamples(
      int &iterations, std::vector<int> &samples)
{
	samples.clear();
	int index0,index1;
	int pNum = _adapter.getNumberCorrespondences();
	int lNum = _adapter.getLineNumberCorrespondences();
	getRandomIndex(pNum,index0);
	getRandomIndex(lNum,index1);
	samples.push_back(index0);
	samples.push_back(index1);
	return;
}


void opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::getmixedSamples(
      int &iterations, std::vector<int> &samples)
{
	samples.clear();
	int pNum = _adapter.getNumberCorrespondences();
	int lNum = _adapter.getLineNumberCorrespondences();
	int index0,index1;
	getRandomIndex(pNum,index0);
	getRandomIndex(pNum+lNum,index1);
	while(index1==index0)
	{
		getRandomIndex(pNum+lNum,index1);
	}
	samples.push_back(index0);
	samples.push_back(index1);
	return;
}

int opengv::sac_problems::
    absolute_pose::AbsoluteLinePoseSacProblem::getPointsNum()
	{
		return _adapter.getNumberCorrespondences();
	}



