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
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/LineCentralAbsoluteAdapter.hpp>
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


int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  //set experiment parameters
  double noise = 0.5;
  double outlierFraction = 0;
  size_t numberPoints = 20;
  size_t numberLines = 20;
  std::vector<int> indices;
  for(int i=0;i<numberPoints;++i)
	  indices.push_back(i);

  //create a random viewpoint pose
  translation_t position = generateRandomTranslation(2.0);
  rotation_t rotation = generateRandomRotation(0.5);
  double mypitch, mxroll, mzyaw;
  Eigen::Vector3d Eulerxyz=opengv::absolute_pose::rotationMatrixToEulerAngles(rotation);
	mxroll = Eulerxyz(0);
	mypitch = Eulerxyz(1);
	mzyaw = Eulerxyz(2);
	std::cout <<"yaw is "<< mzyaw <<"\n,";
  
// 	Matrix3d Rx_Roll;
// 	Rx_Roll <<  1, 0, 0,
// 				0,cos(mxroll),-sin(mxroll),
// 				0,sin(mxroll),cos(mxroll);
// 	Matrix3d Ry_Pitch;
// 	Ry_Pitch << cos(mypitch), 0, sin(mypitch),
// 				0, 1, 0,
// 				-sin(mypitch),0, cos(mypitch);
// 	Matrix3d Ryx=Ry_Pitch*Rx_Roll;
// 	Matrix3d Rz,loc_Rzyx;
// 		Rz << cos(mzyaw),-sin(mzyaw),0,sin(mzyaw),cos(mzyaw),0,0,0,1;
// 		loc_Rzyx = Rz*Ryx;
// 	std::cout <<"loc_Rzyx is \n"<< loc_Rzyx <<"\n";
  //create a fake central camera
  translations_t camOffsets;
  rotations_t camRotations;
  generateCentralCameraSystem( camOffsets, camRotations );
  
  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors;
  points_t points;
  bearingVectors_t linestart_bearingVectors;
  points_t linestart_points;
  bearingVectors_t lineend_bearingVectors;
  points_t lineend_points;
  std::vector<int> camCorrespondences; //unused in the central case!
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberLines, noise, outlierFraction,
      linestart_bearingVectors, linestart_points, camCorrespondences, gt );
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberLines, noise, outlierFraction,
      lineend_bearingVectors, lineend_points, camCorrespondences, gt );


  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );

  //create a central absolute adapter
  absolute_pose::CentralAbsoluteAdapter adapter(
      bearingVectors,
      points,
      rotation);

  //Create an AbsolutePoseSac problem and Ransac
  //The method can be set to KNEIP, GAO or EPNP
  sac::Ransac<sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
      new sac_problems::absolute_pose::AbsolutePoseSacProblem(
      adapter,
      sac_problems::absolute_pose::AbsolutePoseSacProblem::TWOE2P,
		  mypitch, mxroll,indices));
  ransac.sac_model_ = absposeproblem_ptr;
  ransac.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac.max_iterations_ = 50;

  //Run the experiment
  struct timeval tic;
  struct timeval toc;
  gettimeofday( &tic, 0 );
  bool ransac_success = ransac.compute2EModel();
  gettimeofday( &toc, 0 );
  double ransac_time = TIMETODOUBLE(timeval_minus(toc,tic));
  Eigen::Matrix<double, 3, 4> final_model = ransac.model_coefficients_;
   if (ransac_success) {
		// Optional nonlinear model refinement over all inliers.
		absposeproblem_ptr->optimizeModelCoefficients(ransac.inliers_,
														ransac.model_coefficients_,
														final_model);
	}

  //print the results
  std::cout << "the 2E2P ransac results is: " << std::endl;
  std::cout << final_model << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
  std::cout << ransac_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
    std::cout << ransac.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  
    //create a central absolute adapter
  absolute_pose::LineCentralAbsoluteAdapter adapter3(
      bearingVectors,
      points,
      linestart_bearingVectors,
	  lineend_bearingVectors,
      linestart_points,
	  lineend_points,
      rotation);
  adapter3.setR(rotation);
  adapter3.sett(position);

  // for 1p1l
  //Create an AbsolutePoseSac problem and Ransac
  //The method can be set to KNEIP, GAO or EPNP
  sac::Ransac<sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> ransac3;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> absposeproblem_ptr3(
      new sac_problems::absolute_pose::AbsoluteLinePoseSacProblem(
      adapter3,
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::TWOE1P1L,
		  mypitch, mxroll,indices));
  ransac3.sac_model_ = absposeproblem_ptr3;
  ransac3.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac3.max_iterations_ = 50;

  //Run the experiment
//   struct timeval tic;
//   struct timeval toc;
  gettimeofday( &tic, 0 );
  bool ransac_success3 = ransac3.compute1P1LModel();
  gettimeofday( &toc, 0 );
  double ransac3_time = TIMETODOUBLE(timeval_minus(toc,tic));
  Eigen::Matrix<double, 3, 4> final_model3 = ransac3.model_coefficients_;
	  if (ransac_success3) {
		// Optional nonlinear model refinement over all inliers.
// 		absposeproblem_ptr3->optimizeModelCoefficients(ransac3.inliers_,ransac3.line_inliers_,
// 														ransac3.model_coefficients_,
// 														final_model3);
		absposeproblem_ptr3->optimizeModelCoefficients(ransac3.inliers_,
														ransac3.model_coefficients_,
														final_model3);
	}

  //print the results
  std::cout << "the 2E1P1L ransac results is: " << std::endl;
  std::cout << final_model3 << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac3.iterations_ << " iterations and ";
  std::cout << ransac3_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << (ransac3.inliers_.size());
//   std::cout << "the number of line inliers is: " << (ransac3.line_inliers_.size());
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac3.inliers_.size(); i++)
    std::cout << ransac3.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  
  
  //Create an AbsolutePoseSac problem and Ransac
  //The method can be set to KNEIP, GAO or EPNP
  sac::Ransac<sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac2;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr2(
      new sac_problems::absolute_pose::AbsolutePoseSacProblem(
      adapter,
      sac_problems::absolute_pose::AbsolutePoseSacProblem::GAO,indices));
  ransac2.sac_model_ = absposeproblem_ptr2;
  ransac2.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac2.max_iterations_ = 50;

  //Run the experiment
//   struct timeval tic;
//   struct timeval toc;
  gettimeofday( &tic, 0 );
  bool ransac_success2 = ransac2.computeModel();
  gettimeofday( &toc, 0 );
  double ransac2_time = TIMETODOUBLE(timeval_minus(toc,tic));
  Eigen::Matrix<double, 3, 4> final_model2 = ransac2.model_coefficients_;
   if (ransac_success2) {
		// Optional nonlinear model refinement over all inliers.
		absposeproblem_ptr2->optimizeModelCoefficients(ransac2.inliers_,
														ransac2.model_coefficients_,
														final_model2);
	}

  //print the results
  std::cout << "the UPNP ransac results is: " << std::endl;
  std::cout << final_model2 << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac2.iterations_ << " iterations and ";
  std::cout << ransac2_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac2.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac2.inliers_.size(); i++)
    std::cout << ransac2.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  // for 2p1l
  //Create an AbsolutePoseSac problem and Ransac
  //The method can be set to KNEIP, GAO or EPNP
  sac::Ransac<sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> ransac4;
  std::shared_ptr<
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem> absposeproblem_ptr4(
      new sac_problems::absolute_pose::AbsoluteLinePoseSacProblem(
      adapter3,
      sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::TWOEMixed,
	  mypitch, mxroll,indices));
  ransac4.sac_model_ = absposeproblem_ptr4;
  ransac4.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac4.max_iterations_ = 50;

  //Run the experiment
//   struct timeval tic;
//   struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac4.computemixedModel();
  gettimeofday( &toc, 0 );
  double ransac4_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the Mixed ransac results is: " << std::endl;
  std::cout << ransac4.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac4.iterations_ << " iterations and ";
  std::cout << ransac4_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac4.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac4.inliers_.size(); i++)
    std::cout << ransac4.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  
  
  
}
