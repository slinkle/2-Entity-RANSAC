library: OpenGV
pages:   http://laurentkneip.github.io/opengv
brief:   OpenGV is a collection of computer vision methods for solving
         geometric vision problems. It contains absolute-pose, relative-pose,
         triangulation, and point-cloud alignment methods for the calibrated
         case. All problems can be solved with central or non-central cameras,
         and embedded into a random sample consensus or nonlinear optimization
         context. Matlab and Python interfaces are implemented as well. The link
         to the above pages also shows links to precompiled Matlab mex-libraries.
         Please consult the documentation for more information.
author:  Laurent Kneip, The Australian National University
contact: kneip.laurent@gmail.com

https://laurentkneip.github.io/opengv/page_installation.html

git clone https://github.com/laurentkneip/opengv

sudo apt-get install build-essential

sudo apt-get install cmake

sudo apt-get install cmake libeigen3-dev

mkdir build && cd build && cmake .. && make

2-entity RANSAC:
author:  Yanmei Jiao, Zhejiang University
contact: ymjiao@zju.edu.cn

for monocular camera pose estimation: 

1. 2P --> sac_problems::absolute_pose::AbsolutePoseSacProblem::TWOE2P
2. 1P1L --> sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::TWOE1P1L
3. Mixed --> sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::TWOEMixed

for multiple camera pose estimaion:

1. MC2P --> sac_problems::absolute_pose::AbsolutePoseSacProblem::MC2P
2. MC1P1L --> sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::MC1P1L
3. MCMixed --> sac_problems::absolute_pose::AbsoluteLinePoseSacProblem::MCMixed


