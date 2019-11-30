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

/**
 * \file AbsoluteLinePoseSacProblem.hpp
 * \brief Functions for fitting an absolute-pose model to a set of
 *        bearing-vector-point correspondences, using different algorithms
 *        (central and non-central one). Used in a sample-consenus paradigm for
 *        rejecting outlier correspondences.
 */

#ifndef OPENGV_SAC_PROBLEMS_ABSOLUTE_POSE_ABSOLUTELINEPOSESACPROBLEM_HPP_
#define OPENGV_SAC_PROBLEMS_ABSOLUTE_POSE_ABSOLUTELINEPOSESACPROBLEM_HPP_

#include <opengv/sac/LineSampleConsensusProblem.hpp>
#include <opengv/types.hpp>
#include <opengv/absolute_pose/AbsoluteLineAdapterBase.hpp>

/**
 * \brief The namespace of this library.
 */
namespace opengv
{
/**
 * \brief The namespace for the sample consensus problems.
 */
namespace sac_problems
{
/**
 * \brief The namespace for the absolute pose methods.
 */
namespace absolute_pose
{

/**
 * Provides functions for fitting an absolute-pose model to a set of
 * bearing-vector to point correspondences, using different algorithms (central
 * and non-central ones). Used in a sample-consenus paradigm for rejecting
 * outlier correspondences.
 */
class AbsoluteLinePoseSacProblem :
    public sac::LineSampleConsensusProblem<transformation_t>
{
public:
  /** The model we are trying to fit (transformation) */
  typedef transformation_t model_t;
//   typedef transformations_t models_t;
  /** The type of adapter that is expected by the methods */
  typedef opengv::absolute_pose::AbsoluteLineAdapterBase adapter_t;

  /** The possible algorithms for solving this problem */
  typedef enum Algorithm
  {
    TWOE1P1L = 0,
	TWOEMixed = 1,
	MC1P1L=2,
	MCMixed=3,
  } algorithm_t;

  /**
   * \brief Constructor.
   * \param[in] adapter Visitor holding bearing vectors, world points, etc.
   * \param[in] algorithm The algorithm we want to use.
   * \param[in] randomSeed Whether to seed the random number generator with
   *            the current time.
   */
  AbsoluteLinePoseSacProblem(adapter_t & adapter, algorithm_t algorithm,
      bool randomSeed = true) :
      sac::LineSampleConsensusProblem<model_t> (randomSeed),
      _adapter(adapter),
      _algorithm(algorithm)
  {
    setUniformIndices(adapter.getNumberCorrespondences());
    setUniformLineIndices(adapter.getLineNumberCorrespondences());
  };
  
    /**
   * \brief Constructor.
   * \param[in] adapter Visitor holding bearing vectors, world points, etc.
   * \param[in] algorithm The algorithm we want to use.
   * \param[in] mpitch The pitch (y) angle with repect to the inertial frame
   * \param[in] mroll The roll (x) angle with repect to the inertial frame
   * \param[in] randomSeed Whether to seed the random number generator with
   *            the current time.
   */
    AbsoluteLinePoseSacProblem(adapter_t & adapter, algorithm_t algorithm, double mpitch, double mroll,
      bool randomSeed = true) :
      sac::LineSampleConsensusProblem<model_t> (randomSeed),
      _adapter(adapter),
      _algorithm(algorithm),
      _mpitch(mpitch),
      _mroll(mroll)
  {
    setUniformIndices(adapter.getNumberCorrespondences());
    setUniformLineIndices(adapter.getLineNumberCorrespondences());
  };
  
    AbsoluteLinePoseSacProblem(adapter_t & adapter, algorithm_t algorithm, double mpitch, double mroll,
      const std::vector<int> & indices, bool randomSeed = true) :
      sac::LineSampleConsensusProblem<model_t> (randomSeed),
      _adapter(adapter),
      _algorithm(algorithm),
      _mpitch(mpitch),
      _mroll(mroll)
  {
    setIndices(indices);
    setUniformLineIndices(adapter.getLineNumberCorrespondences());
  };
  
  

  
    AbsoluteLinePoseSacProblem(adapter_t & adapter, algorithm_t algorithm, double mpitch, double mroll,double myaw,
      bool randomSeed = true) :
      sac::LineSampleConsensusProblem<model_t> (randomSeed),
      _adapter(adapter),
      _algorithm(algorithm),
      _mpitch(mpitch),
      _mroll(mroll),
      _myaw(myaw)
  {
    setUniformIndices(adapter.getNumberCorrespondences());
    setUniformLineIndices(adapter.getLineNumberCorrespondences());
  };
  
    AbsoluteLinePoseSacProblem(adapter_t & adapter, algorithm_t algorithm, double mpitch, double mroll,double myaw,
      const std::vector<int> & indices, bool randomSeed = true) :
      sac::LineSampleConsensusProblem<model_t> (randomSeed),
      _adapter(adapter),
      _algorithm(algorithm),
      _mpitch(mpitch),
      _mroll(mroll),
      _myaw(myaw)
  {
    setIndices(indices);
    setUniformLineIndices(adapter.getLineNumberCorrespondences());
  };

  /**
   * \brief Constructor.
   * \param[in] adapter Visitor holding bearing vectors, world points, etc.
   * \param[in] algorithm The algorithm we want to use.
   * \param[in] indices A vector of indices to be used from all available
   *                    correspondences.
   * \param[in] randomSeed Whether to seed the random number generator with
   *            the current time.
   */
  AbsoluteLinePoseSacProblem(
      adapter_t & adapter,
      algorithm_t algorithm,
      const std::vector<int> & indices,
      bool randomSeed = true) :
      sac::LineSampleConsensusProblem<model_t> (randomSeed),
      _adapter(adapter),
      _algorithm(algorithm)
  {
    setIndices(indices);
    setUniformLineIndices(adapter.getLineNumberCorrespondences());
  };

  /**
   * Destructor.
   */
  virtual ~AbsoluteLinePoseSacProblem() {};

  /**
   * \brief See parent-class.
   */
  virtual bool computeModelCoefficients(
      const std::vector<int> & indices,
      model_t & outModel) const;
  /**
   * \brief See parent-class.
   */
  bool compute1P1LModelCoefficients(
      const std::vector<int> & indices,
      model_t & outModel1, model_t & outModel2);
  bool compute2PModelCoefficients(
      const std::vector<int> & indices,
      model_t & outModel1, model_t & outModel2);
  
  void get1P1LSamples(int &iterations, std::vector<int> &samples);
  void get2P1LSamples(int &iterations, std::vector<int> &samples);
  void getmixedSamples(int &iterations, std::vector<int> &samples);
  int getPointsNum();

  /**
   * \brief See parent-class.
   */
  virtual void getSelectedDistancesToModel(
      const model_t & model,
      const std::vector<int> & indices,
      std::vector<double> & scores) const;
	  
  virtual void getSelectedDistancesToModel(
      const model_t & model,
      const std::vector<int> & indices,
      const std::vector<int> & line_indices,
      std::vector<double> & scores,
	  std::vector<double> & line_scores) const;

  /**
   * \brief See parent-class.
   */
  virtual void optimizeModelCoefficients(
      const std::vector<int> & inliers,
      const model_t & model,
      model_t & optimized_model);
  

  virtual void optimizeModelCoefficients(
      const std::vector<int> & inliers,
      const std::vector<int> & line_inliers,
      const model_t & model,
      model_t & optimized_model);

  /**
   * \brief See parent-class.
   */ 
  virtual int getSampleSize() const;
  

protected:
  /** The adapter holding all input data */
  adapter_t & _adapter;
  /** The algorithm we are using */
  algorithm_t _algorithm;
  double _mpitch;
  double _mroll;
  double _myaw;
};

}
}
}

#endif  //#ifndef OPENGV_SAC_PROBLEMS_ABSOLUTE_POSE_ABSOLUTEPOSESACPROBLEM_HPP_
