/******************************************************************************
 * Authors:  Laurent Kneip & Paul Furgale                                     *
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

//Note: has been derived from ROS

template<typename P>
opengv::sac::Ransac<P>::Ransac(
    int maxIterations, double threshold, double probability) :
    SampleConsensus<P>(maxIterations, threshold, probability)
{}

template<typename P>
opengv::sac::Ransac<P>::~Ransac(){}


template<typename PROBLEM_T>
bool
opengv::sac::Ransac<PROBLEM_T>::computeModel(
    int debug_verbosity_level)
{
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;

  iterations_ = 0;
  int n_best_inliers_count = -INT_MAX;
  double k = 1.0;

  std::vector<int> selection;
  model_t model_coefficients;

  int n_inliers_count = 0;
  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_iterations_ * 10;
  

  // Iterate
  while( iterations_ < k && skipped_count < max_skip )
  {
    // Get X samples which satisfy the model criteria
    sac_model_->getSamples( iterations_, selection );

    if(selection.empty()) 
    {
      fprintf(stderr,
          "[sm::RandomSampleConsensus::computeModel] No samples could be selected!\n");
      break;
    }
    
    // Search for inliers in the point cloud for the current plane model M
    if(!sac_model_->computeModelCoefficients( selection, model_coefficients ))
    {
      //++iterations_;
      ++ skipped_count;
      continue;
    }


    n_inliers_count = sac_model_->countWithinDistance(
        model_coefficients, threshold_ );
// 	mexPrintf("the inlier num is %d \n",n_inliers_count);

    // Better match ?
    if(n_inliers_count > n_best_inliers_count)
    {
      n_best_inliers_count = n_inliers_count;

      // Save the current model/inlier/coefficients selection as being the best so far
      model_              = selection;
      model_coefficients_ = model_coefficients;

      // Compute the k parameter (k=log(z)/log(1-w^n))
      double w = static_cast<double> (n_best_inliers_count) /
          static_cast<double> (sac_model_->getIndices()->size());
      double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
      p_no_outliers =
          (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by -Inf
      p_no_outliers =
          (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by 0.
      k = log(1.0 - probability_) / log(p_no_outliers);
    }

    ++iterations_;
    if(debug_verbosity_level > 1)
      fprintf(stdout,
          "[sm::RandomSampleConsensus::computeModel] Trial %d out of %f: %d inliers (best is: %d so far).\n",
          iterations_, k, n_inliers_count, n_best_inliers_count );
    if(iterations_ > max_iterations_)
    {
      if(debug_verbosity_level > 0)
        fprintf(stdout,
            "[sm::RandomSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
      break;
    }
  }

  if(debug_verbosity_level > 0)
    fprintf(stdout,
        "[sm::RandomSampleConsensus::computeModel] Model: %zu size, %d inliers.\n",
        model_.size(), n_best_inliers_count );

  if(model_.empty())
  {
    inliers_.clear();
// 	mexPrintf("false 1 \n");
    return (false);
  }

  // Get the set of inliers that correspond to the best model found so far
  sac_model_->selectWithinDistance(
      model_coefficients_, threshold_, inliers_, inlier_distances_to_model_);

  return (true);
}


template<typename PROBLEM_T>
bool
opengv::sac::Ransac<PROBLEM_T>::compute2EModel(
    int debug_verbosity_level)
{
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;

  iterations_ = 0;
  int n_best_inliers_count = -INT_MAX;
  double k = 1.0;

  std::vector<int> selection;
  model_t model_coefficients;

  int n_inliers_count = 0;
  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_iterations_ * 10;

  // Iterate
  while( iterations_ < k && skipped_count < max_skip )
  {
// 	  std::cout << "into iterations" << iterations_ << std::endl;
    // Get X samples which satisfy the model criteria
    sac_model_->getSamples( iterations_, selection );

    if(selection.empty()) 
    {
      fprintf(stderr,
          "[sm::RandomSampleConsensus::computeModel] No samples could be selected!\n");
      break;
    }
    
//     std::cout << "the selected point1 and point2 id is " << selection[0] << "," << selection[1] << std::endl;
    model_t model_coefficients1,model_coefficients2;
    if(!sac_model_->compute2PModelCoefficients( selection, model_coefficients1, model_coefficients2 ))
    {
      //++iterations_;
      ++ skipped_count;
      continue;
    }
    int n_inliers_count1 = sac_model_->countWithinDistance(
        model_coefficients1, threshold_ );
    int n_inliers_count2 = sac_model_->countWithinDistance(
        model_coefficients2, threshold_ );
    if(n_inliers_count1>n_inliers_count2)
    {
      model_coefficients=model_coefficients1;
      n_inliers_count=n_inliers_count1;
    }
    else
    {
      model_coefficients=model_coefficients2;
      n_inliers_count = n_inliers_count2;
    }
//     std::cout << "the model_coefficients is \n" << model_coefficients  << std::endl;

    // Better match ?
    if(n_inliers_count > n_best_inliers_count)
    {
      n_best_inliers_count = n_inliers_count;

      // Save the current model/inlier/coefficients selection as being the best so far
      model_              = selection;
      model_coefficients_ = model_coefficients;

      // Compute the k parameter (k=log(z)/log(1-w^n))
      double w = static_cast<double> (n_best_inliers_count) /
          static_cast<double> (sac_model_->getIndices()->size());
      double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
      p_no_outliers =
          (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by -Inf
      p_no_outliers =
          (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by 0.
      k = log(1.0 - probability_) / log(p_no_outliers);
    }

    ++iterations_;
    if(debug_verbosity_level > 1)
      fprintf(stdout,
          "[sm::RandomSampleConsensus::computeModel] Trial %d out of %f: %d inliers (best is: %d so far).\n",
          iterations_, k, n_inliers_count, n_best_inliers_count );
    if(iterations_ > max_iterations_)
    {
      if(debug_verbosity_level > 0)
        fprintf(stdout,
            "[sm::RandomSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
      break;
    }
  }

  if(debug_verbosity_level > 0)
    fprintf(stdout,
        "[sm::RandomSampleConsensus::computeModel] Model: %zu size, %d inliers.\n",
        model_.size(), n_best_inliers_count );

  if(model_.empty())
  {
    inliers_.clear();
    return (false);
  }

  // Get the set of inliers that correspond to the best model found so far
  sac_model_->selectWithinDistance(
      model_coefficients_, threshold_, inliers_, inlier_distances_to_model_);

  return (true);
}


template<typename PROBLEM_T>
bool
opengv::sac::Ransac<PROBLEM_T>::compute1P1LModel( // 2-entity ransac with 1p1l
    int debug_verbosity_level)
{
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;

  iterations_ = 0;
  int n_best_inliers_count = -INT_MAX;
  double k = 1.0;

  std::vector<int> selection;
  model_t model_coefficients;

  int n_inliers_count = 0;
  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_iterations_ * 10;

  // Iterate
  while( iterations_ < k && skipped_count < max_skip )
  {
    // Get X samples which satisfy the model criteria
    sac_model_->get1P1LSamples( iterations_, selection );

    if(selection.empty()) 
    {
      fprintf(stderr,
          "[sm::RandomSampleConsensus::computeModel] No samples could be selected!\n");
      break;
    }
//     std::cout << "the selected point and line id is " << selection[0] << "," << selection[1] << std::endl;
    
    model_t model_coefficients1,model_coefficients2;
    if(!sac_model_->compute1P1LModelCoefficients( selection, model_coefficients1, model_coefficients2 ))
    {
      //++iterations_;
      ++ skipped_count;
      continue;
    }
//     std::cout << "the first model_coefficients1 = \n" << model_coefficients1 << std::endl;
//     std::cout << "the second model_coefficients2 = \n" << model_coefficients2 << std::endl;
    int n_inliers_count1 = sac_model_->countWithinDistance(
        model_coefficients1, threshold_ );
    int n_inliers_count2 = sac_model_->countWithinDistance(
        model_coefficients2, threshold_ );
    if(n_inliers_count1>n_inliers_count2)
    {
      model_coefficients=model_coefficients1;
      n_inliers_count=n_inliers_count1;
    }
    else
    {
      model_coefficients=model_coefficients2;
      n_inliers_count = n_inliers_count2;
    }

//     std::cout << "the n_inliers_count is " << n_inliers_count <<  std::endl;
// 	if(n_inliers_count)
    // Better match ?
    if(n_inliers_count > n_best_inliers_count)
    {
      n_best_inliers_count = n_inliers_count;

      // Save the current model/inlier/coefficients selection as being the best so far
      model_              = selection;
      model_coefficients_ = model_coefficients;

      // Compute the k parameter (k=log(z)/log(1-w^n))
      double w = static_cast<double> (n_best_inliers_count) /
          static_cast<double> (sac_model_->getIndices()->size());
      double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
      p_no_outliers =
          (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by -Inf
      p_no_outliers =
          (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by 0.
      k = log(1.0 - probability_) / log(p_no_outliers);
    }

    ++iterations_;
    if(debug_verbosity_level > 1)
      fprintf(stdout,
          "[sm::RandomSampleConsensus::computeModel] Trial %d out of %f: %d inliers (best is: %d so far).\n",
          iterations_, k, n_inliers_count, n_best_inliers_count );
    if(iterations_ > max_iterations_)
    {
      if(debug_verbosity_level > 0)
        fprintf(stdout,
            "[sm::RandomSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
      break;
    }
  }

  if(debug_verbosity_level > 0)
    fprintf(stdout,
        "[sm::RandomSampleConsensus::computeModel] Model: %zu size, %d inliers.\n",
        model_.size(), n_best_inliers_count );

  if(model_.empty())
  {
    inliers_.clear();
    return (false);
  }

  // Get the set of inliers that correspond to the best model found so far
  sac_model_->selectWithinDistance(
      model_coefficients_, threshold_, inliers_, inlier_distances_to_model_);

  return (true);
}

template<typename PROBLEM_T>
bool
opengv::sac::Ransac<PROBLEM_T>::computemixedModel( // 2-entity ransac with 1p1l
    int debug_verbosity_level)
{
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;

  iterations_ = 0;
  int n_best_inliers_count = -INT_MAX;
  double k = 1.0;

  std::vector<int> selection;
  model_t model_coefficients;

  int n_inliers_count = 0;
  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_iterations_ * 10;

  // Iterate
  while( iterations_ < k && skipped_count < max_skip )
  {
    // Get X samples which satisfy the model criteria
    sac_model_->getmixedSamples( iterations_, selection );

    if(selection.empty()) 
    {
      fprintf(stderr,
          "[sm::RandomSampleConsensus::computeModel] No samples could be selected!\n");
      break;
    }
//     std::cout << "the selected point and line id is " << selection[0] << "," << selection[1] << std::endl;
    int pNum = sac_model_->getPointsNum();
    model_t model_coefficients1,model_coefficients2;
	int secondFeatureId = selection[1];
	if(secondFeatureId<pNum)
	{
		// 2p method
// 		std::cout << "-------------------------------2p method\n";
		if(!sac_model_->compute2PModelCoefficients( selection, model_coefficients1, model_coefficients2 ))
		{
			//++iterations_;
			++ skipped_count;
			continue;
		}
	}
	else
	{
// 		std::cout << "-------------------------------1p1l method\n";
		selection[1] = secondFeatureId-pNum;
		if(!sac_model_->compute1P1LModelCoefficients( selection, model_coefficients1, model_coefficients2 ))
		{
			//++iterations_;
			++ skipped_count;
			continue;
		}
	}
    
//     std::cout << "the first model_coefficients1 = \n" << model_coefficients1 << std::endl;
//     std::cout << "the second model_coefficients2 = \n" << model_coefficients2 << std::endl;
    int n_inliers_count1 = sac_model_->countWithinDistance(
        model_coefficients1, threshold_ );
    int n_inliers_count2 = sac_model_->countWithinDistance(
        model_coefficients2, threshold_ );
    if(n_inliers_count1>n_inliers_count2)
    {
      model_coefficients=model_coefficients1;
      n_inliers_count=n_inliers_count1;
    }
    else
    {
      model_coefficients=model_coefficients2;
      n_inliers_count = n_inliers_count2;
    }

//     std::cout << "the n_inliers_count is " << n_inliers_count <<  std::endl;
// 	if(n_inliers_count)
    // Better match ?
    if(n_inliers_count > n_best_inliers_count)
    {
      n_best_inliers_count = n_inliers_count;

      // Save the current model/inlier/coefficients selection as being the best so far
      model_              = selection;
      model_coefficients_ = model_coefficients;

      // Compute the k parameter (k=log(z)/log(1-w^n))
      double w = static_cast<double> (n_best_inliers_count) /
          static_cast<double> (sac_model_->getIndices()->size());
      double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
      p_no_outliers =
          (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by -Inf
      p_no_outliers =
          (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by 0.
      k = log(1.0 - probability_) / log(p_no_outliers);
    }

    ++iterations_;
    if(debug_verbosity_level > 1)
      fprintf(stdout,
          "[sm::RandomSampleConsensus::computeModel] Trial %d out of %f: %d inliers (best is: %d so far).\n",
          iterations_, k, n_inliers_count, n_best_inliers_count );
    if(iterations_ > max_iterations_)
    {
      if(debug_verbosity_level > 0)
        fprintf(stdout,
            "[sm::RandomSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
      break;
    }
  }

  if(debug_verbosity_level > 0)
    fprintf(stdout,
        "[sm::RandomSampleConsensus::computeModel] Model: %zu size, %d inliers.\n",
        model_.size(), n_best_inliers_count );

  if(model_.empty())
  {
    inliers_.clear();
    return (false);
  }

  // Get the set of inliers that correspond to the best model found so far
  sac_model_->selectWithinDistance(
      model_coefficients_, threshold_, inliers_, inlier_distances_to_model_);

  return (true);
}

