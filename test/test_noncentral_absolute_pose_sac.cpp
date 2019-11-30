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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <Eigen/Eigen>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/LineNoncentralAbsoluteAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac/Lmeds.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/sac_problems/absolute_pose/AbsoluteLinePoseSacProblem.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;

Eigen::Vector3d rotationMatrixToEulerAngles(Eigen::Matrix3d R)
{
    
    //     assert(cv::isRotationMatrix(R));

	float sy = sqrt(R(0,0) * R(0,0) +  R(1,0) * R(1,0) );
    
	bool singular = sy < 1e-6; // If
    
	float x, y, z;
	if (!singular)
	{
	    x = atan2(R(2,1) , R(2,2));
	    y = atan2(-R(2,0), sy);
	    z = atan2(R(1,0), R(0,0));
	}
	else
	{
	    x = atan2(-R(1,2), R(1,1));
	    y = atan2(-R(2,0), sy);
	    z = 0;
	}
	return Eigen::Vector3d(x, y, z);
	
}

int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  //set experiment parameters
  double noise = 0.5;
  double outlierFraction = 0.1;
//   double noise = 0;
//   double outlierFraction = 0;
  size_t numberPoints = 20;
  size_t numberLines = numberPoints;
  int numberCameras = 2;

  //create a random viewpoint pose
  translation_t position = generateRandomTranslation(2.0);
  rotation_t rotation = generateRandomRotation(0.5);
  double mypitch, mxroll, mzyaw;
  Eigen::Vector3d Eulerxyz=rotationMatrixToEulerAngles(rotation);
	mxroll = Eulerxyz(0);
	mypitch = Eulerxyz(1);
	mzyaw = Eulerxyz(2);
  
  //create a random camera-system
  translations_t camOffsets;
  rotations_t camRotations;
  generateRandomCameraSystem( numberCameras, camOffsets, camRotations );
  
  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors;
  points_t points;
  bearingVectors_t linestart_bearingVectors;
  points_t linestart_points;
  bearingVectors_t lineend_bearingVectors;
  points_t lineend_points;
  std::vector<int> camCorrespondences; //unused in the central case!
  std::vector<int> Line_camCorrespondences; //unused in the central case!
  Eigen::MatrixXd gt(3,numberPoints);
  Eigen::MatrixXd lgt(3,2*numberLines);
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberLines, noise, outlierFraction,
      linestart_bearingVectors, linestart_points,
	  lineend_bearingVectors, lineend_points, Line_camCorrespondences, lgt );

  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );

  //create a non-central absolute adapter
  absolute_pose::NoncentralAbsoluteAdapter adapter(
      bearingVectors,
      camCorrespondences,
      points,
      camOffsets,
      camRotations);
  
      //create a central absolute adapter
  absolute_pose::LineNoncentralAbsoluteAdapter adapter3(
      bearingVectors,
      camCorrespondences,
      points,
      linestart_bearingVectors,
	  lineend_bearingVectors,
      Line_camCorrespondences,
      linestart_points,
	  lineend_points,
      camOffsets,
      camRotations);
  adapter3.setR(rotation);
  adapter3.sett(position);


  //Create a AbsolutePoseSacProblem and Ransac
  //The method is set to GP3P
  sac::Ransac<
      sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
      new sac_problems::absolute_pose::AbsolutePoseSacProblem(
      adapter,
      sac_problems::absolute_pose::AbsolutePoseSacProblem::MC2P,mypitch,mxroll,mzyaw));
  ransac.sac_model_ = absposeproblem_ptr;
  ransac.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac.max_iterations_ = 50;

  //Run the experiment
  struct timeval tic;
  struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac.compute2EModel();
  gettimeofday( &toc, 0 );
  double ransac_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the MC2P ransac results is: " << std::endl;
  std::cout << ransac.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
  std::cout << ransac_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
    std::cout << ransac.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  
  
  //Create a AbsolutePoseSacProblem and Ransac
  //The method is set to GP3P UPNP GPNP
  sac::Ransac<
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> ransac4;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> absposeproblem_ptr4(
      new sac_problems::absolute_pose::AbsoluteLinePoseSacProblem(
      adapter3,
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::MC1P1L,mypitch,mxroll,mzyaw));
  ransac4.sac_model_ = absposeproblem_ptr4;
  ransac4.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac4.max_iterations_ = 50;

  //Run the experiment
//   struct timeval tic;
//   struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac4.compute1P1LModel();
  gettimeofday( &toc, 0 );
  double ransac4_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the MC1P1L ransac results is: " << std::endl;
  std::cout << ransac4.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac4.iterations_ << " iterations and ";
  std::cout << ransac4_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac4.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac4.inliers_.size(); i++)
    std::cout << ransac4.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  
  //Create a AbsolutePoseSacProblem and Ransac
  //The method is set to GP3P UPNP GPNP
  sac::Ransac<
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> ransac5;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> absposeproblem_ptr5(
      new sac_problems::absolute_pose::AbsoluteLinePoseSacProblem(
      adapter3,
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::MCMixed,mypitch,mxroll,mzyaw));
  ransac5.sac_model_ = absposeproblem_ptr5;
  ransac5.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac5.max_iterations_ = 50;

  //Run the experiment
//   struct timeval tic;
//   struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac5.computemixedModel();
  gettimeofday( &toc, 0 );
  double ransac5_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the MCMixed ransac results is: " << std::endl;
  std::cout << ransac5.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac5.iterations_ << " iterations and ";
  std::cout << ransac5_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac5.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac5.inliers_.size(); i++)
    std::cout << ransac5.inliers_[i] << " ";
  std::cout << std::endl << std::endl;

  
  
  
  sac::Ransac<
      sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac3;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr3(
      new sac_problems::absolute_pose::AbsolutePoseSacProblem(
      adapter,
      sac_problems::absolute_pose::AbsolutePoseSacProblem::GPNP));
  ransac3.sac_model_ = absposeproblem_ptr3;
  ransac3.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac3.max_iterations_ = 50;

  //Run the experiment
//   struct timeval tic;
//   struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac3.computeModel();
  gettimeofday( &toc, 0 );
  double ransac3_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the GPNP ransac results is: " << std::endl;
  std::cout << ransac3.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac3.iterations_ << " iterations and ";
  std::cout << ransac3_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac3.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac3.inliers_.size(); i++)
    std::cout << ransac3.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  
  

//   // Create LMedS
//   sac::Lmeds<sac_problems::absolute_pose::AbsolutePoseSacProblem> lmeds;
//   lmeds.sac_model_ = absposeproblem_ptr;
//   lmeds.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
//   lmeds.max_iterations_ = 50;
// 
//   //Run the experiment
//   gettimeofday( &tic, 0 );
//   lmeds.computeModel();
//   gettimeofday( &toc, 0 );
//   double lmeds_time = TIMETODOUBLE(timeval_minus(toc,tic));
// 
//   //print the results
//   std::cout << "the lmeds results is: " << std::endl;
//   std::cout << lmeds.model_coefficients_ << std::endl << std::endl;
//   std::cout << "Lmeds needed " << lmeds.iterations_ << " iterations and ";
//   std::cout << lmeds_time << " seconds" << std::endl << std::endl;
//   std::cout << "the number of inliers is: " << lmeds.inliers_.size();
//   std::cout << std::endl << std::endl;
//   std::cout << "the found inliers are: " << std::endl;
//   for(size_t i = 0; i < lmeds.inliers_.size(); i++)
//     std::cout << lmeds.inliers_[i] << " ";
//   std::cout << std::endl << std::endl;
}
