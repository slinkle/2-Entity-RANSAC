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


#include <opengv/absolute_pose/methods.hpp>
#include <opengv/Indices.hpp>

#include <Eigen/NonLinearOptimization>
#include <Eigen/NumericalDiff>

#include <opengv/absolute_pose/modules/main.hpp>
#include <opengv/absolute_pose/modules/Epnp.hpp>
#include <opengv/OptimizationFunctor.hpp>
#include <opengv/math/cayley.hpp>
#include <opengv/math/quaternion.hpp>
#include <opengv/math/roots.hpp>

#include <iostream>

Eigen::Vector3d opengv::absolute_pose::rotationMatrixToEulerAngles(Eigen::Matrix3d R)
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

Eigen::Matrix3d opengv::absolute_pose::EulerAnglesTorotationMatrix(const double mxRoll, const double myPitch, const double mzYaw)
{
	 Eigen::Matrix3d Rx_Roll;
	Rx_Roll <<  1, 0, 0,
			0,cos(mxRoll),-sin(mxRoll),
			0,sin(mxRoll),cos(mxRoll);
	Eigen::Matrix3d Ry_Pitch;
	Ry_Pitch << cos(myPitch), 0, sin(myPitch),
			0, 1, 0,
			-sin(myPitch),0, cos(myPitch);
	Eigen::Matrix3d Rz_Yaw;
		Rz_Yaw << cos(mzYaw),-sin(mzYaw),0,sin(mzYaw),cos(mzYaw),0,0,0,1;
	Eigen::Matrix3d R=Rz_Yaw*Ry_Pitch*Rx_Roll;
	return R;
}

void opengv::absolute_pose::pose_estimate_3d3d(const vector< Vector3d >& vpts1, const vector< Vector3d >& vpts2, Matrix3d& R, Vector3d& t)
{
	Vector3d p1(0,0,0);
	Vector3d p2(0,0,0);     // center of mass
	int N = vpts1.size();
	for ( int i=0; i<N; i++ )
	{
// 		cout<<"------"<<endl;
// 		std::cout<<vpts1[i]<<endl;
// 		cout<<"------"<<endl;
// 		std::cout<<vpts2[i]<<endl;
		p1 += vpts1[i];
		p2 += vpts2[i];
	}
	p1 = ( (p1) /  N);
	p2 = ( (p2) / N);
// 	cout << "p1=" << p1<<endl;
// 	cout << "p2=" << p2<<endl;
	Matrix3Xd  q1 (3, N );
	Matrix3Xd q2 ( 3,N ); // remove the center
	for ( int i=0; i<N; i++ )
	{
		q1.col(i) = vpts1[i] - p1;
		q2.col(i) = vpts2[i] - p2;
	}
// 	cout << "q1=" << q1<<endl;
// 	cout << "q2=" << q2<<endl;

	// compute q1*q2^T
	Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
	W = q1*q2.transpose();
// 	for ( int i=0; i<N; i++ )
// 	{
// 		W += q1[i]*q2[i].transpose();
// 	}
// 		cout<<"W="<<W<<endl;

	// SVD on W
	Eigen::JacobiSVD<Eigen::Matrix3d> svd ( W, Eigen::ComputeFullU|Eigen::ComputeFullV );
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV();
// 	cout << "U=" << U<<endl;
// 	cout << "V=" << V<<endl;
	
// 	if (U.determinant() * V.determinant() < 0)
// 	{
// 		for (int x = 0; x < 3; ++x)
// 		{
// 			U(x, 2) *= -1;
// 		}
// 	}
	

	R = U* ( V.transpose() );
	t = p1 - R * p2;
// 	cout<<"R="<<R<<endl;
// 	cout<<"t="<<t<<endl;
	return;
}


opengv::translation_t
opengv::absolute_pose::p2p(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return p2p( adapter, indices[0], indices[1] );
}

opengv::translation_t
opengv::absolute_pose::p2p(
    const AbsoluteAdapterBase & adapter,
    size_t index0,
    size_t index1)
{
  Eigen::Vector3d e1 = adapter.getBearingVector(index0);
  Eigen::Vector3d e3 = adapter.getBearingVector(index1);
  e3 = e1.cross(e3);
  e3 = e3/e3.norm();
  Eigen::Vector3d e2 = e3.cross(e1);

  rotation_t T;
  T.row(0) = e1.transpose();
  T.row(1) = e2.transpose();
  T.row(2) = e3.transpose();

  Eigen::Vector3d n1 = adapter.getPoint(index1) - adapter.getPoint(index0);
  n1 = n1/n1.norm();
  Eigen::Vector3d n3;
  if( (fabs(n1[0]) > fabs(n1[1])) && (fabs(n1[0]) > fabs(n1[2])) )
  {
    n3[1] = 1.0;
    n3[2] = 0.0;
    n3[0] = -n1[1]/n1[0];
  }
  else
  {
    if( (fabs(n1[1]) > fabs(n1[0])) && (fabs(n1[1]) > fabs(n1[2])) )
    {
      n3[2] = 1.0;
      n3[0] = 0.0;
      n3[1] = -n1[2]/n1[1];
    }
    else
    {
      n3[0] = 1.0;
      n3[1] = 0.0;
      n3[2] = -n1[0]/n1[2];
    }
  }
  n3 = n3 / n3.norm();
  Eigen::Vector3d n2 = n3.cross(n1);

  rotation_t N;
  N.row(0) = n1.transpose();
  N.row(1) = n2.transpose();
  N.row(2) = n3.transpose();

  Eigen::Matrix3d Q = T * adapter.getR().transpose() * N.transpose();
  Eigen::Vector3d temp1 = adapter.getPoint(index1) - adapter.getPoint(index0);
  double d_12 = temp1.norm();

  Eigen::Vector3d temp2 = adapter.getBearingVector(index1);
  double cos_beta = e1.dot(temp2);
  double b = 1/( 1 - pow( cos_beta, 2 ) ) - 1;

  if( cos_beta < 0 )
    b = -sqrt(b);
  else
    b = sqrt(b);

  double temp3 = d_12 * ( Q(1,0) * b - Q(0,0) );

  translation_t solution = -temp3 * Q.row(0).transpose();
  solution = adapter.getPoint(index0) + N.transpose()*solution;

  if(
    solution(0,0) != solution(0,0) ||
    solution(1,0) != solution(1,0) ||
    solution(2,0) != solution(2,0) )
    solution = Eigen::Vector3d::Zero();

  return solution;
}

opengv::transformations_t 
opengv::absolute_pose::twoe1p1l(
    const AbsoluteLineAdapterBase & adapter,
    const double & mpitch,
    const double & mroll,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return twoe1p1l( adapter, mpitch, mroll, indices[0], indices[1] );
}

opengv::transformations_t 
opengv::absolute_pose::twoe1p1l(
    const AbsoluteLineAdapterBase & adapter,
    const double & myPitch,
    const double & mxRoll,
    size_t index0,
    size_t index1)
{
	transformations_t solutions;
		// compute the pose using 1 point and 1 line 2D-3D matches 
	Eigen::Vector3d e1 = adapter.getBearingVector(index0);
	Eigen::Vector3d D10 = e1/e1(2);
	
	Eigen::Vector3d ls = adapter.getLinestartBearingVector(index1);
	Eigen::Vector3d le = adapter.getLineendBearingVector(index1);
	Eigen::Vector3d D30_ = ls/ls(2);
	Eigen::Vector3d D40_ = le/le(2);
	
	Eigen::Vector3d p3d=adapter.getPoint(index0);
	
	Vector3d L30 = adapter.getLinestartPoint(index1);
	Vector3d L40 = adapter.getLineendPoint(index1);
	
	Eigen::Vector3d l3d = L30.cross(L40);
	if (abs(l3d.dot(p3d))<1e-04)
	{
		return solutions;
// 		std::cout << "point is on the line\n";
	}
	
		
	Vector3d C0(0,0,0);
	
	Vector3d d30 = D30_-C0;
	d30 /= d30.norm();
	Vector3d d40 = D40_-C0;
	d40 /= d40.norm();
	
	Vector3d D30 = C0+d30;
	Vector3d D40 = C0+d40/(d30.transpose()*d40);
	
	Vector3d C(0,0,-1);
	Vector3d D3(0,0,0);
	double cosd34 = d30.transpose()*d40;
	if (cosd34>1)
		cosd34=1;
	else
		if (cosd34<-1)
			cosd34=-1;
		
	Vector3d D4(tan(acos(cosd34)),0,0);
	Matrix3d RC1_C0;
	Vector3d tC1_C0;
	vector<Vector3d> pts1,pts2;
	pts1.push_back(C0);
	pts1.push_back(D30);
	pts1.push_back(D40);
	pts2.push_back(C);
	pts2.push_back(D3);
	pts2.push_back(D4);
	pose_estimate_3d3d(pts1,pts2,RC1_C0,tC1_C0);
// 		cout << "after 3d3d pose estimate :\n";
// 		cout << RC1_C0 << endl;
// 		cout << tC1_C0 << endl;
	
	Vector3d D1 = RC1_C0.inverse()*D10-RC1_C0.transpose()*tC1_C0;
	double a1 = D1(0)/(D1(2)+1);
	double b1 = D1(1)/(D1(2)+1);
	
	//W0 coordinates
	Vector3d P10 = p3d;
	
	// W1 coordinates
	Vector3d P1(0,0,0);
	Matrix3d RW0_W1 = Matrix3d::Identity(3,3);
	Vector3d tW0_W1=P1-P10;
	Vector3d L3=RW0_W1*L30+tW0_W1;
	Vector3d L4=RW0_W1*L40+tW0_W1;
	Matrix4d TW0_W1 = Matrix4d::Identity(4,4);
	Matrix4d TC1_C0 = Matrix4d::Identity(4,4);
	TW0_W1.block(0,0,3,3) = RW0_W1;
	TW0_W1.block(0,3,3,1) = tW0_W1;
	TC1_C0.block(0,0,3,3) = RC1_C0;
	TC1_C0.block(0,3,3,1) = tC1_C0;
	Matrix4d TW1_W0 = TW0_W1.inverse();
	Matrix4d TC0_C1 = TC1_C0.inverse();
	
	double x3=L3(0);double y3=L3(1);double z3=L3(2);
	double x4=L4(0);double y4=L4(1);double z4=L4(2);
	
	Matrix3d Rx_Roll;
	Rx_Roll <<  1, 0, 0,
				0,cos(mxRoll),-sin(mxRoll),
				0,sin(mxRoll),cos(mxRoll);
	Matrix3d Ry_Pitch;
	Ry_Pitch << cos(myPitch), 0, sin(myPitch),
				0, 1, 0,
				-sin(myPitch),0, cos(myPitch);
	Matrix3d Ryx=Ry_Pitch*Rx_Roll;
	
	MatrixX3d h_a(1,3),h_b(1,3),g_a(1,3),g_b(1,3),k(1,3);
	h_a = -Ryx.row(1)*RC1_C0;
	h_b = Ryx.row(0)*RC1_C0;
	g_a = h_b;
	g_b = -h_a;
	k = Ryx.row(2)*RC1_C0;
	
	double T1_a=TW1_W0(0,3)*h_a(0)+TW1_W0(1,3)*g_a(0); //sin
	double T1_b=TW1_W0(0,3)*h_b(0)+TW1_W0(1,3)*g_b(0);
	double T1_c=TW1_W0(2,3)*k(0)+TC0_C1(0,3);



	double R21_a=TW1_W0(0,0)*h_a(1)+TW1_W0(1,0)*g_a(1); //sin
	double R21_b=TW1_W0(0,0)*h_b(1)+TW1_W0(1,0)*g_b(1);
	double R21_c=TW1_W0(2,0)*k(1);

	double R22_a=TW1_W0(0,1)*h_a(1)+TW1_W0(1,1)*g_a(1); //sin
	double R22_b=TW1_W0(0,1)*h_b(1)+TW1_W0(1,1)*g_b(1);
	double R22_c=TW1_W0(2,1)*k(1);

	double R23_a=TW1_W0(0,2)*h_a(1)+TW1_W0(1,2)*g_a(1); //sin
	double R23_b=TW1_W0(0,2)*h_b(1)+TW1_W0(1,2)*g_b(1);
	double R23_c=TW1_W0(2,2)*k(1);

	double T2_a=TW1_W0(0,3)*h_a(1)+TW1_W0(1,3)*g_a(1); //sin
	double T2_b=TW1_W0(0,3)*h_b(1)+TW1_W0(1,3)*g_b(1);
	double T2_c=TW1_W0(2,3)*k(1)+TC0_C1(1,3);
	
	double T3_a=TW1_W0(0,3)*h_a(2)+TW1_W0(1,3)*g_a(2); //sin
	double T3_b=TW1_W0(0,3)*h_b(2)+TW1_W0(1,3)*g_b(2);
	double T3_c=TW1_W0(2,3)*k(2)+TC0_C1(2,3);
	
	double coefb=R21_b*(x3-x4)+R22_b*(y3-y4)+R23_b*(z3-z4);//sin
	double coefa=R21_a*(x3-x4)+R22_a*(y3-y4)+R23_a*(z3-z4);//cos
	double coefc=R21_c*(x3-x4)+R22_c*(y3-y4)+R23_c*(z3-z4);
	
// 	double mzYaw = -0.302449;
// 	cout << "Let's check the gt :\n" ;
// 	cout << coefb*cos(mzYaw)+coefa*sin(mzYaw)+coefc << endl;
// 	cout << "the gt alpha is " << mzYaw << endl;
// 	cout << "coefa = " << coefa << ", coefb = " << coefb << ", coefc = " << coefc << endl;
// 	
	
	double my_alpha;
	double my_alpha2;
	if (coefa != 0)
	{
		if (coefa*coefa+coefb*coefb-coefc*coefc<0)
		{
// 				cout << "can not solve zYaw from this Point and Line\n";
			return solutions;
		}
		double cos_alpha = (-coefb*coefc+abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha=-coefc/coefa-coefb/coefa*cos_alpha;
		my_alpha=atan2(sin_alpha,cos_alpha);
// 		cout << "let us check ours:\n";
// 		cout << coefa*sin(my_alpha)+coefb*cos(my_alpha)+coefc;
// 		cout << endl;
// 		cout << "my_alpha = " << my_alpha << endl;
		
		double cos_alpha2=(-coefb*coefc-abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha2 = -coefc/coefa-coefb/coefa*cos_alpha2;
		my_alpha2=atan2(sin_alpha2,cos_alpha2);
// 		cout << "let us check ours 2:\n";
// 		cout << coefa*sin(my_alpha2)+coefb*cos(my_alpha2)+coefc;
// 		cout << endl;
// 		cout << "my_alpha2 = " << my_alpha2 << endl;
	}
	else
	{
		if (coefb != 0)
		{
			my_alpha = acos(-coefc/coefb);
			my_alpha2 = -my_alpha;
// 			cout << "my_alpha = " << my_alpha << endl;
// 			cout << "my_alpha2 = " << my_alpha2 << endl;
		}
		else
		{
// 				cout << "cannot solve alpha because a==b==0\n";
			return solutions;
		}
	}
	
	vector<double> thetas;
	thetas.push_back(my_alpha);
	thetas.push_back(my_alpha2);
	
	// for each theta(yaw)
	for (int i=0; i < thetas.size(); ++i)
	{
		double alpha = thetas[i];
		double f_d=Ryx(0,0)*cos(alpha)-Ryx(1,0)*sin(alpha);
		double f_e=Ryx(1,0)*cos(alpha)+Ryx(0,0)*sin(alpha);
		double f_f=Ryx(2,0);
		double m_d=Ryx(0,1)*cos(alpha)-Ryx(1,1)*sin(alpha);
		double m_e=Ryx(1,1)*cos(alpha)+Ryx(0,1)*sin(alpha);
		double m_f=Ryx(2,1);
		double n_d=Ryx(0,2)*cos(alpha)-Ryx(1,2)*sin(alpha);
		double n_e=Ryx(1,2)*cos(alpha)+Ryx(0,2)*sin(alpha);
		double n_f=Ryx(2,2);
		
		
		double T1_d=-TC0_C1(0,0)*f_d-TC0_C1(0,1)*m_d-TC0_C1(0,2)*n_d;
		double T1_e=-TC0_C1(0,0)*f_e-TC0_C1(0,1)*m_e-TC0_C1(0,2)*n_e;
		double T1_f=-TC0_C1(0,0)*f_f-TC0_C1(0,1)*m_f-TC0_C1(0,2)*n_f;
		double T1_g=T1_a*sin(alpha)+T1_b*cos(alpha)+T1_c;
		
		double T2_d=-TC0_C1(1,0)*f_d-TC0_C1(1,1)*m_d-TC0_C1(1,2)*n_d;
		double T2_e=-TC0_C1(1,0)*f_e-TC0_C1(1,1)*m_e-TC0_C1(1,2)*n_e;
		double T2_f=-TC0_C1(1,0)*f_f-TC0_C1(1,1)*m_f-TC0_C1(1,2)*n_f;
		double T2_g=T2_a*sin(alpha)+T2_b*cos(alpha)+T2_c;
		
		double T3_d=-TC0_C1(2,0)*f_d-TC0_C1(2,1)*m_d-TC0_C1(2,2)*n_d;
		double T3_e=-TC0_C1(2,0)*f_e-TC0_C1(2,1)*m_e-TC0_C1(2,2)*n_e;
		double T3_f=-TC0_C1(2,0)*f_f-TC0_C1(2,1)*m_f-TC0_C1(2,2)*n_f;
		double T3_g=T3_a*sin(alpha)+T3_b*cos(alpha)+T3_c;
		
		
		
		double coefd1=T2_d;
		double coefe1=T2_e;
		double coeff1=T2_f;
		double coefg1=(R21_a*sin(alpha)+R21_b*cos(alpha)+R21_c)*x3+(R22_a*sin(alpha)+R22_b*cos(alpha)+R22_c)*y3
						+(R23_a*sin(alpha)+R23_b*cos(alpha)+R23_c)*z3+T2_g;
		double coefd2=a1*T2_d-b1*T1_d;double coefe2=a1*T2_e-b1*T1_e;
		double coeff2=a1*T2_f-b1*T1_f;double coefg2=a1*(T2_g)-b1*(T1_g);
		double coefd3=b1*T3_d-T2_d;double coefe3=b1*T3_e-T2_e;
		double coeff3=b1*T3_f-T2_f;double coefg3=b1*T3_g-T2_g+b1;
		
		Matrix3d A;
// 		A << coefd1, coefe1, coeff1, coefd2, coefe2, coeff2, coefd3, coefe3, coeff3;
		A(0,0) = coefd1; A(0,1)=coefe1; A(0,2) = coeff1;
		A(1,0) = coefd2; A(1,1)=coefe2; A(1,2) = coeff2;
		A(2,0) = coefd3; A(2,1)=coefe3; A(2,2) = coeff3;
		Vector3d B;
		B << -coefg1, -coefg2, -coefg3;
		
		
// 		Vector3d loc_t;
// 		loc_t<<1.36075,-0.422468,1.1324;
// 		cout << "Let's check the gt: \n";
// 			cout << A*B << endl;
// 			cout << A.determinant() << endl;
		
		if (abs(A.determinant())<0.00001)
		{
// 			cout << "cout << A is singular can not solve t\n";
// 			cout << A.determinant() << endl;
			return solutions;
		}
		Vector3d loc_tj = A.inverse()*B;
		
		Matrix3d Rz,loc_Rzyx;
		Rz << cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1;
		loc_Rzyx = Rz*Ryx;
		Eigen::Matrix4d T_G_C = Eigen::Matrix4d::Identity(); // camera to world
		T_G_C.block(0,0,3,3) = loc_Rzyx;
		T_G_C.block(0,3,3,1) = loc_tj;
		rotation_t rotation = T_G_C.block(0,0,3,3);
		translation_t translation = T_G_C.block(0,3,3,1);
		transformation_t transformation;
		transformation.block<3,3>(0,0) = rotation;
		transformation.col(3) = translation;
		solutions.push_back(transformation);
		
	}
		

	return solutions;
}


opengv::transformations_t 
opengv::absolute_pose::twoe2p(
    const AbsoluteAdapterBase & adapter,
    const double & mpitch,
    const double & mroll,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return twoe2p( adapter, mpitch, mroll, indices[0], indices[1] );
}

opengv::transformations_t 
opengv::absolute_pose::twoe2p(
    const AbsoluteAdapterBase & adapter,
    const double & myPitch,
    const double & mxRoll,
    size_t index0,
    size_t index1)
{
	transformations_t solutions;
  // 3D point and 2D correspondences idx in mvCorrespondencesIdx
    // the aligned pitch and roll is mxroll mypitch
    // the camera intrinsic parameters is mK (Matrix3d)
	Eigen::Vector3d e1 = adapter.getBearingVector(index0);
	Eigen::Vector3d e2 = adapter.getBearingVector(index1);
	Eigen::Vector3d D1 = e1/e1(2);
	Eigen::Vector3d D2 = e2/e2(2);
	Eigen::Vector3d p3d1=adapter.getPoint(index0);
	Eigen::Vector3d p3d2=adapter.getPoint(index1);
	
	if((p3d1-p3d2).norm()<1e-3)
	{
		return solutions;
	}
  
	double a1=D1[0];double b1=D1[1];
	double a2=D2[0];double b2=D2[1];

	// W0 coordinates
	Vector3d P10=p3d1;
	Vector3d P20=p3d2;
	
	//W1 coordinates
	Vector3d P1(0,0,0);
	Matrix3d RW0_W1=Matrix3d::Identity(3,3);
	Vector3d tW0_W1=P1-P10;
	Vector3d P2=RW0_W1*P20+tW0_W1;
	double x2=P2[0]; double y2=P2[1]; double z2=P2[2];
	
	Matrix4d TW0_W1=Matrix4d::Identity(4,4);
	TW0_W1.block(0,3,3,1)=tW0_W1;
	Matrix4d TW1_W0=TW0_W1.inverse();
	Matrix4d TC0_C1=Matrix4d::Identity(4,4);
	Matrix3d RC1_C0=Matrix3d::Identity(3,3);
	
	Matrix3d Rx_Roll;
	Rx_Roll <<  1, 0, 0,
				0,cos(mxRoll),-sin(mxRoll),
				0,sin(mxRoll),cos(mxRoll);
	Matrix3d Ry_Pitch;
	Ry_Pitch << cos(myPitch), 0, sin(myPitch),
				0, 1, 0,
				-sin(myPitch),0, cos(myPitch);
	Matrix3d Ryx=Ry_Pitch*Rx_Roll;
	
	// R11=coefa*sin(mzYaw)+coefb*cos(mzYaw)+coefc
// 		Vector3d h_a=  //1*3
	MatrixX3d h_a(1,3),h_b(1,3),g_a(1,3),g_b(1,3),k(1,3);
	h_a = -Ryx.row(1)*RC1_C0;
	h_b = Ryx.row(0)*RC1_C0;
	g_a = h_b;
	g_b = -h_a;
	k = Ryx.row(2)*RC1_C0;
// 		cout << "h_a = \n";
// 		cout << h_a << endl;
	
	double R11_a = TW1_W0(0,0)*h_a(0)+TW1_W0(1,0)*g_a(0);
	double R11_b=TW1_W0(0,0)*h_b(0)+TW1_W0(1,0)*g_b(0);
	double R11_c=TW1_W0(2,0)*k(0);

	double R12_a=TW1_W0(0,1)*h_a(0)+TW1_W0(1,1)*g_a(0); //sin
	double R12_b=TW1_W0(0,1)*h_b(0)+TW1_W0(1,1)*g_b(0);
	double R12_c=TW1_W0(2,1)*k(0);

	double R13_a=TW1_W0(0,2)*h_a(0)+TW1_W0(1,2)*g_a(0); //sin
	double R13_b=TW1_W0(0,2)*h_b(0)+TW1_W0(1,2)*g_b(0);
	double R13_c=TW1_W0(2,2)*k(0);

	double T1_a=TW1_W0(0,3)*h_a(0)+TW1_W0(1,3)*g_a(0); //sin
	double T1_b=TW1_W0(0,3)*h_b(0)+TW1_W0(1,3)*g_b(0);
	double T1_c=TW1_W0(2,3)*k(0)+TC0_C1(0,3);



	double R21_a=TW1_W0(0,0)*h_a(1)+TW1_W0(1,0)*g_a(1); //sin
	double R21_b=TW1_W0(0,0)*h_b(1)+TW1_W0(1,0)*g_b(1);
	double R21_c=TW1_W0(2,0)*k(1);

	double R22_a=TW1_W0(0,1)*h_a(1)+TW1_W0(1,1)*g_a(1); //sin
	double R22_b=TW1_W0(0,1)*h_b(1)+TW1_W0(1,1)*g_b(1);
	double R22_c=TW1_W0(2,1)*k(1);

	double R23_a=TW1_W0(0,2)*h_a(1)+TW1_W0(1,2)*g_a(1); //sin
	double R23_b=TW1_W0(0,2)*h_b(1)+TW1_W0(1,2)*g_b(1);
	double R23_c=TW1_W0(2,2)*k(1);

	double T2_a=TW1_W0(0,3)*h_a(1)+TW1_W0(1,3)*g_a(1); //sin
	double T2_b=TW1_W0(0,3)*h_b(1)+TW1_W0(1,3)*g_b(1);
	double T2_c=TW1_W0(2,3)*k(1)+TC0_C1(1,3);


	double R31_a=TW1_W0(0,0)*h_a(2)+TW1_W0(1,0)*g_a(2); //sin
	double R31_b=TW1_W0(0,0)*h_b(2)+TW1_W0(1,0)*g_b(2);
	double R31_c=TW1_W0(2,0)*k(2);

	double R32_a=TW1_W0(0,1)*h_a(2)+TW1_W0(1,1)*g_a(2); //sin
	double R32_b=TW1_W0(0,1)*h_b(2)+TW1_W0(1,1)*g_b(2);
	double R32_c=TW1_W0(2,1)*k(2);

	double R33_a=TW1_W0(0,2)*h_a(2)+TW1_W0(1,2)*g_a(2); //sin
	double R33_b=TW1_W0(0,2)*h_b(2)+TW1_W0(1,2)*g_b(2);
	double R33_c=TW1_W0(2,2)*k(2);

	double T3_a=TW1_W0(0,3)*h_a(2)+TW1_W0(1,3)*g_a(2); //sin
	double T3_b=TW1_W0(0,3)*h_b(2)+TW1_W0(1,3)*g_b(2);
	double T3_c=TW1_W0(2,3)*k(2)+TC0_C1(2,3);
	
	double coefa=(b1-b2)*(R11_a*x2+R12_a*y2+R13_a*z2) + (a2-a1)*(R21_a*x2+R22_a*y2+R23_a*z2) + (a1*b2-a2*b1)*(R31_a*x2+R32_a*y2+R33_a*z2);//sin
	double coefb=(b1-b2)*(R11_b*x2+R12_b*y2+R13_b*z2) + (a2-a1)*(R21_b*x2+R22_b*y2+R23_b*z2) + (a1*b2-a2*b1)*(R31_b*x2+R32_b*y2+R33_b*z2);//cos
	double coefc=(b1-b2)*(R11_c*x2+R12_c*y2+R13_c*z2) + (a2-a1)*(R21_c*x2+R22_c*y2+R23_c*z2) + (a1*b2-a2*b1)*(R31_c*x2+R32_c*y2+R33_c*z2);
	
// 	double mzYaw = -0.302449;
// 	cout << "Let's check the gt :\n" ;
// 	cout << coefb*cos(mzYaw)+coefa*sin(mzYaw)+coefc << endl;
// 	cout << "the gt alpha is " << mzYaw << endl;
// 	cout << "coefa = " << coefa << ", coefb = " << coefb << ", coefc = " << coefc << endl;
// 		
	double my_alpha;
	double my_alpha2;
	if (coefa != 0)
	{
		if (coefa*coefa+coefb*coefb-coefc*coefc<0)
		{
// 				cout << "can not solve zYaw from this 2 Points\n";
			return solutions;
		}
		double cos_alpha = (-coefb*coefc+abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha=-coefc/coefa-coefb/coefa*cos_alpha;
		my_alpha=atan2(sin_alpha,cos_alpha);
// 		cout << "let us check ours:\n";
// 		cout << coefa*sin(my_alpha)+coefb*cos(my_alpha)+coefc;
// 		cout << endl;
// 		cout << "my_alpha = " << my_alpha << endl;
		
		double cos_alpha2=(-coefb*coefc-abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha2 = -coefc/coefa-coefb/coefa*cos_alpha2;
		my_alpha2=atan2(sin_alpha2,cos_alpha2);
// 		cout << "let us check ours 2:\n";
// 		cout << coefa*sin(my_alpha2)+coefb*cos(my_alpha2)+coefc;
// 		cout << endl;
// 		cout << "my_alpha2 = " << my_alpha2 << endl;
	}
	else
	{
		if (coefb != 0)
		{
			my_alpha = acos(-coefc/coefb);
			my_alpha2 = -my_alpha;
// 				cout << "my_alpha = " << my_alpha << endl;
// 				cout << "my_alpha2 = " << my_alpha2 << endl;
		}
		else
		{
// 				cout << "cannot solve alpha because a==b==0\n";
			return solutions;
		}
	}
	
	vector<double> thetas;
	thetas.push_back(my_alpha);
	thetas.push_back(my_alpha2);
	
	// for each theta(yaw)
	for (int i=0; i < thetas.size(); ++i)
	{
		double alpha = thetas[i];
		double f_d=Ryx(0,0)*cos(alpha)-Ryx(1,0)*sin(alpha);
		double f_e=Ryx(1,0)*cos(alpha)+Ryx(0,0)*sin(alpha);
		double f_f=Ryx(2,0);
		double m_d=Ryx(0,1)*cos(alpha)-Ryx(1,1)*sin(alpha);
		double m_e=Ryx(1,1)*cos(alpha)+Ryx(0,1)*sin(alpha);
		double m_f=Ryx(2,1);
		double n_d=Ryx(0,2)*cos(alpha)-Ryx(1,2)*sin(alpha);
		double n_e=Ryx(1,2)*cos(alpha)+Ryx(0,2)*sin(alpha);
		double n_f=Ryx(2,2);
		
		
		double T1_d=-TC0_C1(0,0)*f_d-TC0_C1(0,1)*m_d-TC0_C1(0,2)*n_d;
		double T1_e=-TC0_C1(0,0)*f_e-TC0_C1(0,1)*m_e-TC0_C1(0,2)*n_e;
		double T1_f=-TC0_C1(0,0)*f_f-TC0_C1(0,1)*m_f-TC0_C1(0,2)*n_f;
		double T1_g=T1_a*sin(alpha)+T1_b*cos(alpha)+T1_c;
		
		double T2_d=-TC0_C1(1,0)*f_d-TC0_C1(1,1)*m_d-TC0_C1(1,2)*n_d;
		double T2_e=-TC0_C1(1,0)*f_e-TC0_C1(1,1)*m_e-TC0_C1(1,2)*n_e;
		double T2_f=-TC0_C1(1,0)*f_f-TC0_C1(1,1)*m_f-TC0_C1(1,2)*n_f;
		double T2_g=T2_a*sin(alpha)+T2_b*cos(alpha)+T2_c;
		
		double T3_d=-TC0_C1(2,0)*f_d-TC0_C1(2,1)*m_d-TC0_C1(2,2)*n_d;
		double T3_e=-TC0_C1(2,0)*f_e-TC0_C1(2,1)*m_e-TC0_C1(2,2)*n_e;
		double T3_f=-TC0_C1(2,0)*f_f-TC0_C1(2,1)*m_f-TC0_C1(2,2)*n_f;
		double T3_g=T3_a*sin(alpha)+T3_b*cos(alpha)+T3_c;
		
		double coefd1=a2*T2_d-b2*T1_d;double coefe1=a2*T2_e-b2*T1_e;double coeff1=a2*T2_f-b2*T1_f;
		double coefg1=( (R21_b*cos(alpha)+R21_a*sin(alpha)+R21_c)*x2+(R22_b*cos(alpha)+R22_a*sin(alpha)+R22_c)*y2+(R23_b*cos(alpha)+R23_a*sin(alpha)+R23_c)*z2 )*a2 
		- ( (R11_b*cos(alpha)+R11_a*sin(alpha)+R11_c)*x2+(R12_b*cos(alpha)+R12_a*sin(alpha)+R12_c)*y2+(R13_b*cos(alpha)+R13_a*sin(alpha)+R13_c)*z2 )*b2 + (T2_g)*a2-(T1_g)*b2;
		double coefd2=a1*T2_d-b1*T1_d;double coefe2=a1*T2_e-b1*T1_e;double coeff2=a1*T2_f-b1*T1_f;
		double coefg2=a1*(T2_g)-b1*(T1_g);
		double coefd3=a1*T3_d-T1_d;double coefe3=a1*T3_e-T1_e;double coeff3=a1*T3_f-T1_f;
		double coefg3=a1*(T3_g)-(T1_g);
		
		Matrix3d A;
		A << coefd1, coefe1, coeff1, coefd2, coefe2, coeff2, coefd3, coefe3, coeff3;
		Vector3d B;
		B << -coefg1, -coefg2, -coefg3;
// 		Vector3d loc_t;
// 		loc_t<<1.36075,-0.422468,1.1324;
// 		cout << "Let's check the gt: \n";
// 		cout << A*loc_t-B << endl;
		
		if (abs(A.determinant())<0.00001)
		{
// 				cout << "cout << A is singular can not solve t\n";
			return solutions;
		}
		Vector3d loc_tj = A.inverse()*B;
// 			cout << "loc_t = " << endl;
// 			cout << loc_t << endl;
// 			cout << "the estimated loc_tj = \n";
// 			cout << loc_tj << endl;
		Matrix3d Rz,loc_Rzyx;
		Rz << cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1;
		loc_Rzyx = Rz*Ryx;
		Eigen::Matrix4d T_G_C = Eigen::Matrix4d::Identity(); // camera to world
		T_G_C.block(0,0,3,3) = loc_Rzyx;
		T_G_C.block(0,3,3,1) = loc_tj;
		rotation_t rotation = T_G_C.block(0,0,3,3);
		translation_t translation = T_G_C.block(0,3,3,1);
		transformation_t transformation;
		transformation.block<3,3>(0,0) = rotation;
		transformation.col(3) = translation;
		solutions.push_back(transformation);
		
	}
	return solutions;
}


opengv::transformations_t 
opengv::absolute_pose::twoe2p(
    const AbsoluteLineAdapterBase & adapter,
    const double & mpitch,
    const double & mroll,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return twoe2p( adapter, mpitch, mroll, indices[0], indices[1] );
}

opengv::transformations_t 
opengv::absolute_pose::twoe2p(
    const AbsoluteLineAdapterBase & adapter,
    const double & myPitch,
    const double & mxRoll,
    size_t index0,
    size_t index1)
{
	transformations_t solutions;
  // 3D point and 2D correspondences idx in mvCorrespondencesIdx
    // the aligned pitch and roll is mxroll mypitch
    // the camera intrinsic parameters is mK (Matrix3d)
	Eigen::Vector3d e1 = adapter.getBearingVector(index0);
	Eigen::Vector3d e2 = adapter.getBearingVector(index1);
	Eigen::Vector3d D1 = e1/e1(2);
	Eigen::Vector3d D2 = e2/e2(2);
	Eigen::Vector3d p3d1=adapter.getPoint(index0);
	Eigen::Vector3d p3d2=adapter.getPoint(index1);
	
	if((p3d1-p3d2).norm()<1e-3)
	{
		return solutions;
	}
  
	double a1=D1[0];double b1=D1[1];
	double a2=D2[0];double b2=D2[1];

	// W0 coordinates
	Vector3d P10=p3d1;
	Vector3d P20=p3d2;
	
	//W1 coordinates
	Vector3d P1(0,0,0);
	Matrix3d RW0_W1=Matrix3d::Identity(3,3);
	Vector3d tW0_W1=P1-P10;
	Vector3d P2=RW0_W1*P20+tW0_W1;
	double x2=P2[0]; double y2=P2[1]; double z2=P2[2];
	
	Matrix4d TW0_W1=Matrix4d::Identity(4,4);
	TW0_W1.block(0,3,3,1)=tW0_W1;
	Matrix4d TW1_W0=TW0_W1.inverse();
	Matrix4d TC0_C1=Matrix4d::Identity(4,4);
	Matrix3d RC1_C0=Matrix3d::Identity(3,3);
	
	Matrix3d Rx_Roll;
	Rx_Roll <<  1, 0, 0,
				0,cos(mxRoll),-sin(mxRoll),
				0,sin(mxRoll),cos(mxRoll);
	Matrix3d Ry_Pitch;
	Ry_Pitch << cos(myPitch), 0, sin(myPitch),
				0, 1, 0,
				-sin(myPitch),0, cos(myPitch);
	Matrix3d Ryx=Ry_Pitch*Rx_Roll;
	
	// R11=coefa*sin(mzYaw)+coefb*cos(mzYaw)+coefc
// 		Vector3d h_a=  //1*3
	MatrixX3d h_a(1,3),h_b(1,3),g_a(1,3),g_b(1,3),k(1,3);
	h_a = -Ryx.row(1)*RC1_C0;
	h_b = Ryx.row(0)*RC1_C0;
	g_a = h_b;
	g_b = -h_a;
	k = Ryx.row(2)*RC1_C0;
// 		cout << "h_a = \n";
// 		cout << h_a << endl;
	
	double R11_a = TW1_W0(0,0)*h_a(0)+TW1_W0(1,0)*g_a(0);
	double R11_b=TW1_W0(0,0)*h_b(0)+TW1_W0(1,0)*g_b(0);
	double R11_c=TW1_W0(2,0)*k(0);

	double R12_a=TW1_W0(0,1)*h_a(0)+TW1_W0(1,1)*g_a(0); //sin
	double R12_b=TW1_W0(0,1)*h_b(0)+TW1_W0(1,1)*g_b(0);
	double R12_c=TW1_W0(2,1)*k(0);

	double R13_a=TW1_W0(0,2)*h_a(0)+TW1_W0(1,2)*g_a(0); //sin
	double R13_b=TW1_W0(0,2)*h_b(0)+TW1_W0(1,2)*g_b(0);
	double R13_c=TW1_W0(2,2)*k(0);

	double T1_a=TW1_W0(0,3)*h_a(0)+TW1_W0(1,3)*g_a(0); //sin
	double T1_b=TW1_W0(0,3)*h_b(0)+TW1_W0(1,3)*g_b(0);
	double T1_c=TW1_W0(2,3)*k(0)+TC0_C1(0,3);



	double R21_a=TW1_W0(0,0)*h_a(1)+TW1_W0(1,0)*g_a(1); //sin
	double R21_b=TW1_W0(0,0)*h_b(1)+TW1_W0(1,0)*g_b(1);
	double R21_c=TW1_W0(2,0)*k(1);

	double R22_a=TW1_W0(0,1)*h_a(1)+TW1_W0(1,1)*g_a(1); //sin
	double R22_b=TW1_W0(0,1)*h_b(1)+TW1_W0(1,1)*g_b(1);
	double R22_c=TW1_W0(2,1)*k(1);

	double R23_a=TW1_W0(0,2)*h_a(1)+TW1_W0(1,2)*g_a(1); //sin
	double R23_b=TW1_W0(0,2)*h_b(1)+TW1_W0(1,2)*g_b(1);
	double R23_c=TW1_W0(2,2)*k(1);

	double T2_a=TW1_W0(0,3)*h_a(1)+TW1_W0(1,3)*g_a(1); //sin
	double T2_b=TW1_W0(0,3)*h_b(1)+TW1_W0(1,3)*g_b(1);
	double T2_c=TW1_W0(2,3)*k(1)+TC0_C1(1,3);


	double R31_a=TW1_W0(0,0)*h_a(2)+TW1_W0(1,0)*g_a(2); //sin
	double R31_b=TW1_W0(0,0)*h_b(2)+TW1_W0(1,0)*g_b(2);
	double R31_c=TW1_W0(2,0)*k(2);

	double R32_a=TW1_W0(0,1)*h_a(2)+TW1_W0(1,1)*g_a(2); //sin
	double R32_b=TW1_W0(0,1)*h_b(2)+TW1_W0(1,1)*g_b(2);
	double R32_c=TW1_W0(2,1)*k(2);

	double R33_a=TW1_W0(0,2)*h_a(2)+TW1_W0(1,2)*g_a(2); //sin
	double R33_b=TW1_W0(0,2)*h_b(2)+TW1_W0(1,2)*g_b(2);
	double R33_c=TW1_W0(2,2)*k(2);

	double T3_a=TW1_W0(0,3)*h_a(2)+TW1_W0(1,3)*g_a(2); //sin
	double T3_b=TW1_W0(0,3)*h_b(2)+TW1_W0(1,3)*g_b(2);
	double T3_c=TW1_W0(2,3)*k(2)+TC0_C1(2,3);
	
	double coefa=(b1-b2)*(R11_a*x2+R12_a*y2+R13_a*z2) + (a2-a1)*(R21_a*x2+R22_a*y2+R23_a*z2) + (a1*b2-a2*b1)*(R31_a*x2+R32_a*y2+R33_a*z2);//sin
	double coefb=(b1-b2)*(R11_b*x2+R12_b*y2+R13_b*z2) + (a2-a1)*(R21_b*x2+R22_b*y2+R23_b*z2) + (a1*b2-a2*b1)*(R31_b*x2+R32_b*y2+R33_b*z2);//cos
	double coefc=(b1-b2)*(R11_c*x2+R12_c*y2+R13_c*z2) + (a2-a1)*(R21_c*x2+R22_c*y2+R23_c*z2) + (a1*b2-a2*b1)*(R31_c*x2+R32_c*y2+R33_c*z2);
	
// 	double mzYaw = -0.302449;
// 	cout << "Let's check the gt :\n" ;
// 	cout << coefb*cos(mzYaw)+coefa*sin(mzYaw)+coefc << endl;
// 	cout << "the gt alpha is " << mzYaw << endl;
// 	cout << "coefa = " << coefa << ", coefb = " << coefb << ", coefc = " << coefc << endl;
// 		
	double my_alpha;
	double my_alpha2;
	if (coefa != 0)
	{
		if (coefa*coefa+coefb*coefb-coefc*coefc<0)
		{
// 				cout << "can not solve zYaw from this 2 Points\n";
			return solutions;
		}
		double cos_alpha = (-coefb*coefc+abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha=-coefc/coefa-coefb/coefa*cos_alpha;
		my_alpha=atan2(sin_alpha,cos_alpha);
// 		cout << "let us check ours:\n";
// 		cout << coefa*sin(my_alpha)+coefb*cos(my_alpha)+coefc;
// 		cout << endl;
// 		cout << "my_alpha = " << my_alpha << endl;
		
		double cos_alpha2=(-coefb*coefc-abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha2 = -coefc/coefa-coefb/coefa*cos_alpha2;
		my_alpha2=atan2(sin_alpha2,cos_alpha2);
// 		cout << "let us check ours 2:\n";
// 		cout << coefa*sin(my_alpha2)+coefb*cos(my_alpha2)+coefc;
// 		cout << endl;
// 		cout << "my_alpha2 = " << my_alpha2 << endl;
	}
	else
	{
		if (coefb != 0)
		{
			my_alpha = acos(-coefc/coefb);
			my_alpha2 = -my_alpha;
// 				cout << "my_alpha = " << my_alpha << endl;
// 				cout << "my_alpha2 = " << my_alpha2 << endl;
		}
		else
		{
// 				cout << "cannot solve alpha because a==b==0\n";
			return solutions;
		}
	}
	
	vector<double> thetas;
	thetas.push_back(my_alpha);
	thetas.push_back(my_alpha2);
	
	// for each theta(yaw)
	for (int i=0; i < thetas.size(); ++i)
	{
		double alpha = thetas[i];
		double f_d=Ryx(0,0)*cos(alpha)-Ryx(1,0)*sin(alpha);
		double f_e=Ryx(1,0)*cos(alpha)+Ryx(0,0)*sin(alpha);
		double f_f=Ryx(2,0);
		double m_d=Ryx(0,1)*cos(alpha)-Ryx(1,1)*sin(alpha);
		double m_e=Ryx(1,1)*cos(alpha)+Ryx(0,1)*sin(alpha);
		double m_f=Ryx(2,1);
		double n_d=Ryx(0,2)*cos(alpha)-Ryx(1,2)*sin(alpha);
		double n_e=Ryx(1,2)*cos(alpha)+Ryx(0,2)*sin(alpha);
		double n_f=Ryx(2,2);
		
		
		double T1_d=-TC0_C1(0,0)*f_d-TC0_C1(0,1)*m_d-TC0_C1(0,2)*n_d;
		double T1_e=-TC0_C1(0,0)*f_e-TC0_C1(0,1)*m_e-TC0_C1(0,2)*n_e;
		double T1_f=-TC0_C1(0,0)*f_f-TC0_C1(0,1)*m_f-TC0_C1(0,2)*n_f;
		double T1_g=T1_a*sin(alpha)+T1_b*cos(alpha)+T1_c;
		
		double T2_d=-TC0_C1(1,0)*f_d-TC0_C1(1,1)*m_d-TC0_C1(1,2)*n_d;
		double T2_e=-TC0_C1(1,0)*f_e-TC0_C1(1,1)*m_e-TC0_C1(1,2)*n_e;
		double T2_f=-TC0_C1(1,0)*f_f-TC0_C1(1,1)*m_f-TC0_C1(1,2)*n_f;
		double T2_g=T2_a*sin(alpha)+T2_b*cos(alpha)+T2_c;
		
		double T3_d=-TC0_C1(2,0)*f_d-TC0_C1(2,1)*m_d-TC0_C1(2,2)*n_d;
		double T3_e=-TC0_C1(2,0)*f_e-TC0_C1(2,1)*m_e-TC0_C1(2,2)*n_e;
		double T3_f=-TC0_C1(2,0)*f_f-TC0_C1(2,1)*m_f-TC0_C1(2,2)*n_f;
		double T3_g=T3_a*sin(alpha)+T3_b*cos(alpha)+T3_c;
		
		double coefd1=a2*T2_d-b2*T1_d;double coefe1=a2*T2_e-b2*T1_e;double coeff1=a2*T2_f-b2*T1_f;
		double coefg1=( (R21_b*cos(alpha)+R21_a*sin(alpha)+R21_c)*x2+(R22_b*cos(alpha)+R22_a*sin(alpha)+R22_c)*y2+(R23_b*cos(alpha)+R23_a*sin(alpha)+R23_c)*z2 )*a2 
		- ( (R11_b*cos(alpha)+R11_a*sin(alpha)+R11_c)*x2+(R12_b*cos(alpha)+R12_a*sin(alpha)+R12_c)*y2+(R13_b*cos(alpha)+R13_a*sin(alpha)+R13_c)*z2 )*b2 + (T2_g)*a2-(T1_g)*b2;
		double coefd2=a1*T2_d-b1*T1_d;double coefe2=a1*T2_e-b1*T1_e;double coeff2=a1*T2_f-b1*T1_f;
		double coefg2=a1*(T2_g)-b1*(T1_g);
		double coefd3=a1*T3_d-T1_d;double coefe3=a1*T3_e-T1_e;double coeff3=a1*T3_f-T1_f;
		double coefg3=a1*(T3_g)-(T1_g);
		
		Matrix3d A;
		A << coefd1, coefe1, coeff1, coefd2, coefe2, coeff2, coefd3, coefe3, coeff3;
		Vector3d B;
		B << -coefg1, -coefg2, -coefg3;
// 		Vector3d loc_t;
// 		loc_t<<1.36075,-0.422468,1.1324;
// 		cout << "Let's check the gt: \n";
// 		cout << A*loc_t-B << endl;
		
		if (abs(A.determinant())<0.00001)
		{
// 				cout << "cout << A is singular can not solve t\n";
			return solutions;
		}
		Vector3d loc_tj = A.inverse()*B;
// 			cout << "loc_t = " << endl;
// 			cout << loc_t << endl;
// 			cout << "the estimated loc_tj = \n";
// 			cout << loc_tj << endl;
		Matrix3d Rz,loc_Rzyx;
		Rz << cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1;
		loc_Rzyx = Rz*Ryx;
		Eigen::Matrix4d T_G_C = Eigen::Matrix4d::Identity(); // camera to world
		T_G_C.block(0,0,3,3) = loc_Rzyx;
		T_G_C.block(0,3,3,1) = loc_tj;
		rotation_t rotation = T_G_C.block(0,0,3,3);
		translation_t translation = T_G_C.block(0,3,3,1);
		transformation_t transformation;
		transformation.block<3,3>(0,0) = rotation;
		transformation.col(3) = translation;
		solutions.push_back(transformation);
		
	}
	return solutions;
}

opengv::transformations_t 
opengv::absolute_pose::mc1p1l(
    const AbsoluteLineAdapterBase & adapter,
    const double & mpitch,
    const double & mroll,
    const double & mzYaw,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return mc1p1l( adapter, mpitch, mroll,mzYaw, indices[0], indices[1] );
}

opengv::transformations_t 
opengv::absolute_pose::mc1p1l(
    const AbsoluteLineAdapterBase & adapter,
    const double & myPitch,
    const double & mxRoll,
    const double & mzYaw,
    size_t index0,
    size_t index1)
{
	transformations_t solutions;
	rotation_t R1 = adapter.getCamRotation(index0); // T_B_Cp camera in world
	translation_t t1 = adapter.getCamOffset(index0);
	Eigen::Matrix4d T_B_C1 = Eigen::Matrix4d::Identity();
	T_B_C1.block(0,0,3,3) = R1;
	T_B_C1.block(0,3,3,1) = t1;

	rotation_t R2 = adapter.getLineCamRotation(index1); // T_B_C
	translation_t t2 = adapter.getLineCamOffset(index1);
	Eigen::Matrix4d T_B_C2 = Eigen::Matrix4d::Identity();
	T_B_C2.block(0,0,3,3) = R2;
	T_B_C2.block(0,3,3,1) = t2;

	Eigen::Matrix4d T_C1_C2 = T_B_C2.inverse()*T_B_C1; // C1 in C2
	Eigen::Matrix3d R_C1_C2 = T_C1_C2.block<3,3>(0,0);
	Eigen::Vector3d t_C1_C2 = T_C1_C2.block<3,1>(0,3);
	Matrix3d lR = R_C1_C2;
	Vector3d lt = t_C1_C2;
	Matrix3d pR = Matrix3d::Identity();
	Vector3d pt(0,0,0);
	
	
	Eigen::Matrix3d Rx_Roll;
	Rx_Roll <<  1, 0, 0,
		0,cos(mxRoll),-sin(mxRoll),
		0,sin(mxRoll),cos(mxRoll);
	Eigen::Matrix3d Ry_Pitch;
	Ry_Pitch << cos(myPitch), 0, sin(myPitch),
		0, 1, 0,
		-sin(myPitch),0, cos(myPitch);
	Eigen::Matrix3d Rz_Yaw;
	Rz_Yaw << cos(mzYaw),-sin(mzYaw),0,sin(mzYaw),cos(mzYaw),0,0,0,1;
	Eigen::Matrix3d R_W_B=Rz_Yaw*Ry_Pitch*Rx_Roll;

	//  Eigen::Matrix3d R_W_B; // restore from eular
	Eigen::Matrix3d Rgt_W_C = R_W_B*R1; // B to world * point cam to B -> point cam to world
	//  Eigen::Matrix3d Rgt_W_C = Tgt_W_C.block<3,3>(0,0);
	//  rotationMatrixToEulerAngles()
	Eigen::Vector3d Eulerxyz=rotationMatrixToEulerAngles(Rgt_W_C);
	double mxroll = Eulerxyz(0);
	double mypitch = Eulerxyz(1);
// 	double mxroll = mxRoll;
// 	double mypitch = myPitch;
	
// 	Matrix3d gtBR = adapter.getR();
// 	Vector3d gtBt = adapter.gett();
	
	
	Vector3d C0(0,0,0);
	
		// compute the pose using 1 point and 1 line 2D-3D matches 
	Eigen::Vector3d e1 = adapter.getBearingVector(index0);
	Eigen::Vector3d D10 = e1/e1(2);
	
	Eigen::Vector3d ls = adapter.getLinestartBearingVector(index1);
	Eigen::Vector3d le = adapter.getLineendBearingVector(index1);
	Eigen::Vector3d D30_ = ls/ls(2);
	Eigen::Vector3d D40_ = le/le(2);
	
	Eigen::Vector3d p3d=adapter.getPoint(index0);
	
	Vector3d L30 = adapter.getLinestartPoint(index1);
	Vector3d L40 = adapter.getLineendPoint(index1);
	
	Eigen::Vector3d l3d = L30.cross(L40);
	if (abs(l3d.dot(p3d))<1e-04)
	{
		return solutions;
// 		std::cout << "point is on the line\n";
	}
	
	Vector3d d30 = D30_-C0;
	d30 /= d30.norm();
	Vector3d d40 = D40_-C0;
	d40 /= d40.norm();
	
	Vector3d D30 = C0+d30;
	Vector3d D40 = C0+d40/(d30.transpose()*d40);
	
	Vector3d C(0,0,-1);
	Vector3d D3(0,0,0);
	double cosd34 = d30.transpose()*d40;
	if (cosd34>1)
		cosd34=1;
	else
		if (cosd34<-1)
			cosd34=-1;
		
	Vector3d D4(tan(acos(cosd34)),0,0);
	Matrix3d RCl_C0;
	Vector3d tCl_C0;
	vector<Vector3d> pts1,pts2;
	pts1.push_back(C0);
	pts1.push_back(D30);
	pts1.push_back(D40);
	pts2.push_back(C);
	pts2.push_back(D3);
	pts2.push_back(D4);
	pose_estimate_3d3d(pts1,pts2,RCl_C0,tCl_C0);
// 		cout << "after 3d3d pose estimate :\n";
// 		cout << RCl_C0 << endl;
// 		cout << tCl_C0 << endl;
	
	
	Matrix3d RCp_C0 = Matrix3d::Identity(3,3);
	Vector3d tCp_C0(0,0,1);
	Vector3d D1 = RCp_C0.inverse()*D10-RCp_C0.transpose()*tCp_C0;
	double a1 = D1(0)/(D1(2)+1);
	double b1 = D1(1)/(D1(2)+1);
	
	//W0 coordinates
	Vector3d P10 = p3d;
	// W1 coordinates
	Vector3d P1(0,0,0);
	Matrix3d RW0_W1 = Matrix3d::Identity(3,3);
	Vector3d tW0_W1=P1-P10;
	Vector3d L3=RW0_W1*L30+tW0_W1;
	Vector3d L4=RW0_W1*L40+tW0_W1;
	Matrix4d TW0_W1 = Matrix4d::Identity(4,4);
	TW0_W1.block(0,0,3,3) = RW0_W1;
	TW0_W1.block(0,3,3,1) = tW0_W1;
	Matrix4d TW1_W0 = TW0_W1.inverse();
	
	Matrix4d TCl_C0 = Matrix4d::Identity(4,4);
	TCl_C0.block(0,0,3,3) = RCl_C0;
	TCl_C0.block(0,3,3,1) = tCl_C0;
	Eigen::Matrix4d TC1_C2 = Eigen::Matrix4d::Identity();
	TC1_C2.block(0,0,3,3) = lR;
	TC1_C2.block(0,3,3,1) = lt;
	Matrix4d TC0_Cl = TCl_C0.inverse()*TC1_C2;
	TCl_C0 = TC0_Cl.inverse();
	RCl_C0 = TCl_C0.block(0,0,3,3);
	
	Matrix4d TCp_C0 = Matrix4d::Identity(4,4);
	TCp_C0.block(0,0,3,3) = RCp_C0;
	TCp_C0.block(0,3,3,1) = tCp_C0;
	TC1_C2.block(0,0,3,3) = pR;
	TC1_C2.block(0,3,3,1) = pt;
	Matrix4d TC0_Cp = TCp_C0.inverse()*TC1_C2;
	TCp_C0 = TC0_Cp.inverse();
	RCp_C0 = TCp_C0.block(0,0,3,3);
	
	double x3=L3(0);double y3=L3(1);double z3=L3(2);
	double x4=L4(0);double y4=L4(1);double z4=L4(2);
	
// 	Matrix3d Rx_Roll;
	Rx_Roll <<  1, 0, 0,
				0,cos(mxroll),-sin(mxroll),
				0,sin(mxroll),cos(mxroll);
// 	Matrix3d Ry_Pitch;
	Ry_Pitch << cos(mypitch), 0, sin(mypitch),
				0, 1, 0,
				-sin(mypitch),0, cos(mypitch);
	Matrix3d Ryx=Ry_Pitch*Rx_Roll;
	
	MatrixX3d h_a(1,3),h_b(1,3),g_a(1,3),g_b(1,3),k(1,3);
	h_a = -Ryx.row(1)*RCl_C0;
	h_b = Ryx.row(0)*RCl_C0;
	g_a = h_b;
	g_b = -h_a;
	k = Ryx.row(2)*RCl_C0;
	
	MatrixX3d ph_a(1,3),ph_b(1,3),pg_a(1,3),pg_b(1,3),pk(1,3);
	ph_a = -Ryx.row(1)*RCp_C0;
	ph_b = Ryx.row(0)*RCp_C0;
	pg_a = ph_b;
	pg_b = -ph_a;
	pk = Ryx.row(2)*RCp_C0;
	
	double pT1_a=TW1_W0(0,3)*ph_a(0)+TW1_W0(1,3)*pg_a(0); //sin
	double pT1_b=TW1_W0(0,3)*ph_b(0)+TW1_W0(1,3)*pg_b(0);
	double pT1_c=TW1_W0(2,3)*pk(0)+TC0_Cp(0,3);



	double R21_a=TW1_W0(0,0)*h_a(1)+TW1_W0(1,0)*g_a(1); //sin
	double R21_b=TW1_W0(0,0)*h_b(1)+TW1_W0(1,0)*g_b(1);
	double R21_c=TW1_W0(2,0)*k(1);

	double R22_a=TW1_W0(0,1)*h_a(1)+TW1_W0(1,1)*g_a(1); //sin
	double R22_b=TW1_W0(0,1)*h_b(1)+TW1_W0(1,1)*g_b(1);
	double R22_c=TW1_W0(2,1)*k(1);

	double R23_a=TW1_W0(0,2)*h_a(1)+TW1_W0(1,2)*g_a(1); //sin
	double R23_b=TW1_W0(0,2)*h_b(1)+TW1_W0(1,2)*g_b(1);
	double R23_c=TW1_W0(2,2)*k(1);

	double T2_a=TW1_W0(0,3)*h_a(1)+TW1_W0(1,3)*g_a(1); //sin
	double T2_b=TW1_W0(0,3)*h_b(1)+TW1_W0(1,3)*g_b(1);
	double T2_c=TW1_W0(2,3)*k(1)+TC0_Cl(1,3);
	
	double pT2_a=TW1_W0(0,3)*ph_a(1)+TW1_W0(1,3)*pg_a(1); //sin
	double pT2_b=TW1_W0(0,3)*ph_b(1)+TW1_W0(1,3)*pg_b(1);
	double pT2_c=TW1_W0(2,3)*pk(1)+TC0_Cp(1,3);
	
	double pT3_a=TW1_W0(0,3)*ph_a(2)+TW1_W0(1,3)*pg_a(2); //sin
	double pT3_b=TW1_W0(0,3)*ph_b(2)+TW1_W0(1,3)*pg_b(2);
	double pT3_c=TW1_W0(2,3)*pk(2)+TC0_Cp(2,3);
	
	double coefb=R21_b*(x3-x4)+R22_b*(y3-y4)+R23_b*(z3-z4);//sin
	double coefa=R21_a*(x3-x4)+R22_a*(y3-y4)+R23_a*(z3-z4);//cos
	double coefc=R21_c*(x3-x4)+R22_c*(y3-y4)+R23_c*(z3-z4);
	
// 	double mzYaw = -0.302449;
// 	cout << "Let's check the gt :\n" ;
// 	cout << coefb*cos(mzYaw)+coefa*sin(mzYaw)+coefc << endl;
// 	cout << "the gt alpha is " << mzYaw << endl;
// 	cout << "coefa = " << coefa << ", coefb = " << coefb << ", coefc = " << coefc << endl;
// 	
	
	double my_alpha;
	double my_alpha2;
	if (coefa != 0)
	{
		if (coefa*coefa+coefb*coefb-coefc*coefc<0)
		{
// 				cout << "can not solve zYaw from this Point and Line\n";
			return solutions;
		}
		double cos_alpha = (-coefb*coefc+abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha=-coefc/coefa-coefb/coefa*cos_alpha;
		my_alpha=atan2(sin_alpha,cos_alpha);
// 		cout << "let us check ours:\n";
// 		cout << coefa*sin(my_alpha)+coefb*cos(my_alpha)+coefc;
// 		cout << endl;
// 		cout << "my_alpha = " << my_alpha << endl;
		
		double cos_alpha2=(-coefb*coefc-abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
		double sin_alpha2 = -coefc/coefa-coefb/coefa*cos_alpha2;
		my_alpha2=atan2(sin_alpha2,cos_alpha2);
// 		cout << "let us check ours 2:\n";
// 		cout << coefa*sin(my_alpha2)+coefb*cos(my_alpha2)+coefc;
// 		cout << endl;
// 		cout << "my_alpha2 = " << my_alpha2 << endl;
	}
	else
	{
		if (coefb != 0)
		{
			my_alpha = acos(-coefc/coefb);
			my_alpha2 = -my_alpha;
// 			cout << "my_alpha = " << my_alpha << endl;
// 			cout << "my_alpha2 = " << my_alpha2 << endl;
		}
		else
		{
// 				cout << "cannot solve alpha because a==b==0\n";
			return solutions;
		}
	}
	
	vector<double> thetas;
	thetas.push_back(my_alpha);
	thetas.push_back(my_alpha2);
	
	// for each theta(yaw)
	for (int i=0; i < thetas.size(); ++i)
	{
		double alpha = thetas[i];
		double f_d=Ryx(0,0)*cos(alpha)-Ryx(1,0)*sin(alpha);
		double f_e=Ryx(1,0)*cos(alpha)+Ryx(0,0)*sin(alpha);
		double f_f=Ryx(2,0);
		double m_d=Ryx(0,1)*cos(alpha)-Ryx(1,1)*sin(alpha);
		double m_e=Ryx(1,1)*cos(alpha)+Ryx(0,1)*sin(alpha);
		double m_f=Ryx(2,1);
		double n_d=Ryx(0,2)*cos(alpha)-Ryx(1,2)*sin(alpha);
		double n_e=Ryx(1,2)*cos(alpha)+Ryx(0,2)*sin(alpha);
		double n_f=Ryx(2,2);
		
		
		double pT1_d=-TC0_Cp(0,0)*f_d-TC0_Cp(0,1)*m_d-TC0_Cp(0,2)*n_d;
		double pT1_e=-TC0_Cp(0,0)*f_e-TC0_Cp(0,1)*m_e-TC0_Cp(0,2)*n_e;
		double pT1_f=-TC0_Cp(0,0)*f_f-TC0_Cp(0,1)*m_f-TC0_Cp(0,2)*n_f;
		double pT1_g=pT1_a*sin(alpha)+pT1_b*cos(alpha)+pT1_c;
		
		double T2_d=-TC0_Cl(1,0)*f_d-TC0_Cl(1,1)*m_d-TC0_Cl(1,2)*n_d;
		double T2_e=-TC0_Cl(1,0)*f_e-TC0_Cl(1,1)*m_e-TC0_Cl(1,2)*n_e;
		double T2_f=-TC0_Cl(1,0)*f_f-TC0_Cl(1,1)*m_f-TC0_Cl(1,2)*n_f;
		double T2_g=T2_a*sin(alpha)+T2_b*cos(alpha)+T2_c;
		
		double pT2_d=-TC0_Cp(1,0)*f_d-TC0_Cp(1,1)*m_d-TC0_Cp(1,2)*n_d;
		double pT2_e=-TC0_Cp(1,0)*f_e-TC0_Cp(1,1)*m_e-TC0_Cp(1,2)*n_e;
		double pT2_f=-TC0_Cp(1,0)*f_f-TC0_Cp(1,1)*m_f-TC0_Cp(1,2)*n_f;
		double pT2_g=pT2_a*sin(alpha)+pT2_b*cos(alpha)+pT2_c;
		
		double pT3_d=-TC0_Cp(2,0)*f_d-TC0_Cp(2,1)*m_d-TC0_Cp(2,2)*n_d;
		double pT3_e=-TC0_Cp(2,0)*f_e-TC0_Cp(2,1)*m_e-TC0_Cp(2,2)*n_e;
		double pT3_f=-TC0_Cp(2,0)*f_f-TC0_Cp(2,1)*m_f-TC0_Cp(2,2)*n_f;
		double pT3_g=pT3_a*sin(alpha)+pT3_b*cos(alpha)+pT3_c;
		
		
		
		double coefd1=T2_d;
		double coefe1=T2_e;
		double coeff1=T2_f;
		double coefg1=(R21_a*sin(alpha)+R21_b*cos(alpha)+R21_c)*x3+
						(R22_a*sin(alpha)+R22_b*cos(alpha)+R22_c)*y3
						+(R23_a*sin(alpha)+R23_b*cos(alpha)+R23_c)*z3+T2_g;
		double coefd2=a1*pT2_d-b1*pT1_d;double coefe2=a1*pT2_e-b1*pT1_e;
		double coeff2=a1*pT2_f-b1*pT1_f;double coefg2=a1*(pT2_g)-b1*(pT1_g);
		double coefd3=b1*pT3_d-pT2_d;double coefe3=b1*pT3_e-pT2_e;
		double coeff3=b1*pT3_f-pT2_f;double coefg3=b1*pT3_g-pT2_g+b1;
		
		Matrix3d A;
		A << coefd1, coefe1, coeff1, coefd2, coefe2, coeff2, coefd3, coefe3, coeff3;
		Vector3d B;
		B << -coefg1, -coefg2, -coefg3;
		
		
// 		Vector3d loc_t;
// 		loc_t<<1.36075,-0.422468,1.1324;
// 		cout << "Let's check the gt: \n";
// 			cout << A*B << endl;
// 			cout << A.determinant() << endl;
		
		if (abs(A.determinant())<0.00001)
		{
// 			cout << "cout << A is singular can not solve t\n";
// 			cout << A.determinant() << endl;
			return solutions;
		}
		Vector3d loc_tj = A.inverse()*B;
		Eigen::Matrix3d Rz,loc_Rzyx;
		Rz << cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1;
		loc_Rzyx = Rz*Ryx; // camera to world
		Eigen::Matrix4d T_G_C = Eigen::Matrix4d::Identity();
		T_G_C.block(0,0,3,3) = loc_Rzyx;
		T_G_C.block(0,3,3,1) = loc_tj;
		Eigen::Matrix4d T_G_B = T_G_C*T_B_C1.inverse();
		rotation_t rotation = T_G_B.block(0,0,3,3);
		translation_t translation = T_G_B.block(0,3,3,1);
		transformation_t transformation;
		transformation.block<3,3>(0,0) = rotation;
		transformation.col(3) = translation;
		solutions.push_back(transformation);
		
	}
		

	return solutions;
}


opengv::transformations_t 
opengv::absolute_pose::mc2p(
    const AbsoluteAdapterBase & adapter,
    const double & mpitch,
    const double & mroll,
	const double & myaw,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return mc2p( adapter, mpitch, mroll, myaw, indices[0], indices[1] );
}

opengv::transformations_t 
opengv::absolute_pose::mc2p(
    const AbsoluteAdapterBase & adapter,
    const double & myPitch,
    const double & mxRoll,
	const double & mzYaw,
    size_t index0,
    size_t index1)
{
  transformations_t solutions;
  // 3D point and 2D correspondences idx in mvCorrespondencesIdx
    // the aligned pitch and roll is mxroll mypitch
    // the camera intrinsic parameters is mK (Matrix3d)
  
// 	std::cout << "get into mc2p solution\n";
  
  rotation_t R1 = adapter.getCamRotation(index0); // T_B_C camera in world
  translation_t t1 = adapter.getCamOffset(index0);
  Eigen::Matrix4d T_B_C1 = Eigen::Matrix4d::Identity();
  T_B_C1.block(0,0,3,3) = R1;
  T_B_C1.block(0,3,3,1) = t1;
  
  rotation_t R2 = adapter.getCamRotation(index1); // T_B_C
  translation_t t2 = adapter.getCamOffset(index1);
  Eigen::Matrix4d T_B_C2 = Eigen::Matrix4d::Identity();
  T_B_C2.block(0,0,3,3) = R2;
  T_B_C2.block(0,3,3,1) = t2;
  
  Eigen::Matrix4d T_C1_C2 = T_B_C2.inverse()*T_B_C1; // C1 in C2
  Eigen::Matrix3d R_C1_C2 = T_C1_C2.block<3,3>(0,0);
  Eigen::Vector3d t_C1_C2 = T_C1_C2.block<3,1>(0,3);

  
  Eigen::Vector3d e1 = adapter.getBearingVector(index0);
  Eigen::Vector3d e2 = adapter.getBearingVector(index1);
  Eigen::Vector3d D1 = e1/e1(2);
  Eigen::Vector3d D2 = e2/e2(2);
  Eigen::Vector3d p3d1=adapter.getPoint(index0);
  Eigen::Vector3d p3d2=adapter.getPoint(index1);
  
  if((p3d1-p3d2).norm()<1e-3)
  {
    return solutions;
  }
  
  // C0 coordinates
  double a1=D1[0];double b1=D1[1];
  double a2=D2[0];double b2=D2[1];
  
  if(a1==0)
    return solutions;
  
  // W0 coordinates
  Eigen::Vector3d P10=p3d1;
  Eigen::Vector3d P20=p3d2;
  
  //W1 coordinates
  Eigen::Vector3d P1(0,0,0);
  Eigen::Matrix3d RW0_W1=Eigen::Matrix3d::Identity(3,3);
  Eigen::Vector3d tW0_W1=P1-P10;
  Eigen::Vector3d P2=RW0_W1*P20+tW0_W1;
  double x2=P2[0]; double y2=P2[1]; double z2=P2[2];
  
  Eigen::Matrix4d TW0_W1=Eigen::Matrix4d::Identity(4,4);
  TW0_W1.block(0,3,3,1)=tW0_W1;

  Eigen::Matrix4d TW1_W0=TW0_W1.inverse();
  Eigen::Matrix4d TC0_C1=Eigen::Matrix4d::Identity(4,4);
  Eigen::Matrix3d RC1_C0=Eigen::Matrix3d::Identity(3,3);
  
//   // debug
//   Eigen::Matrix3d Rgt;
//   Rgt <<  0.874865 , 0.396999 , 0.277495,
// -0.272976 , 0.877371, -0.394594,
//  -0.40012,  0.269467,  0.875952;
//  Eigen::Vector3d tgt(  1.36075,-0.422468,1.1324 );
//  Eigen::Matrix4d Tgt_G_B=Eigen::Matrix4d::Identity(4,4);
//  Tgt_G_B.block<3,3>(0,0)=Rgt;
//  Tgt_G_B.block<3,1>(0,3)=tgt;
// //  Eigen::Matrix4d T_G_B = T_G_C*T_B_C1.inverse();
//  Eigen::Matrix4d Tgt_W_C; // c1 in w
//  Tgt_W_C = Tgt_G_B*T_B_C1;
 
 
 Eigen::Matrix3d Rx_Roll;
  Rx_Roll <<  1, 0, 0,
        0,cos(mxRoll),-sin(mxRoll),
        0,sin(mxRoll),cos(mxRoll);
  Eigen::Matrix3d Ry_Pitch;
  Ry_Pitch << cos(myPitch), 0, sin(myPitch),
        0, 1, 0,
        -sin(myPitch),0, cos(myPitch);
  Eigen::Matrix3d Rz_Yaw;
    Rz_Yaw << cos(mzYaw),-sin(mzYaw),0,sin(mzYaw),cos(mzYaw),0,0,0,1;
  Eigen::Matrix3d R_W_B=Rz_Yaw*Ry_Pitch*Rx_Roll;
 
//  Eigen::Matrix3d R_W_B; // restore from eular
 Eigen::Matrix3d Rgt_W_C = R_W_B*R1;
//  Eigen::Matrix3d Rgt_W_C = Tgt_W_C.block<3,3>(0,0);
//  rotationMatrixToEulerAngles()
 Eigen::Vector3d Eulerxyz=rotationMatrixToEulerAngles(Rgt_W_C);
	double mxroll = Eulerxyz(0);
	double mypitch = Eulerxyz(1);
	
	// 	// debug
// 	double mzyaw = Eulerxyz(2);
// 	std::cout << "the C1_W pitch roll and yaw is " << mypitch << "," << mxroll << "," <<mzyaw << endl;
// 	myPitch=mypitch;
// 	mxRoll=mxroll;


//   Eigen::Matrix4d TW1_C1 = TC0_C1*Tgt_W_C.inverse()*TW1_W0;
// //   Eigen::Matrix4d TW1_C2 = T_C1_C2*TW1_C1;
//   
//    std::cout << "let us check the frist eq:\n";
// 	std::cout << a1*TW1_C1(1,3)-b1*TW1_C1(0,3)<< ", "
// 	<< a1*TW1_C1(2,3)-TW1_C1(0,3) << std::endl;
	
	
	
  
  Eigen::Matrix3d Rx;
  Rx <<  1, 0, 0,
        0,cos(mxroll),-sin(mxroll),
        0,sin(mxroll),cos(mxroll);
  Eigen::Matrix3d Ry;
  Ry << cos(mypitch), 0, sin(mypitch),
        0, 1, 0,
        -sin(mypitch),0, cos(mypitch);
  Eigen::Matrix3d Ryx=Ry*Rx;
  
  // R11=coefa*sin(mzYaw)+coefb*cos(mzYaw)+coefc
//    Vector3d h_a=  //1*3
  Eigen::MatrixX3d h_a(1,3),h_b(1,3),g_a(1,3),g_b(1,3),k(1,3);
  h_a = -Ryx.row(1)*RC1_C0;
  h_b = Ryx.row(0)*RC1_C0;
  g_a = h_b;
  g_b = -h_a;
  k = Ryx.row(2)*RC1_C0;
//    cout << "h_a = \n";
//    cout << h_a << endl;
  
  double R11_a=TW1_W0(0,0)*h_a(0)+TW1_W0(1,0)*g_a(0);
  double R11_b=TW1_W0(0,0)*h_b(0)+TW1_W0(1,0)*g_b(0);
  double R11_c=TW1_W0(2,0)*k(0);

  double R12_a=TW1_W0(0,1)*h_a(0)+TW1_W0(1,1)*g_a(0); //sin
  double R12_b=TW1_W0(0,1)*h_b(0)+TW1_W0(1,1)*g_b(0);
  double R12_c=TW1_W0(2,1)*k(0);

  double R13_a=TW1_W0(0,2)*h_a(0)+TW1_W0(1,2)*g_a(0); //sin
  double R13_b=TW1_W0(0,2)*h_b(0)+TW1_W0(1,2)*g_b(0);
  double R13_c=TW1_W0(2,2)*k(0);

  double T1_a=TW1_W0(0,3)*h_a(0)+TW1_W0(1,3)*g_a(0); //sin
  double T1_b=TW1_W0(0,3)*h_b(0)+TW1_W0(1,3)*g_b(0);
  double T1_c=TW1_W0(2,3)*k(0)+TC0_C1(0,3);



  double R21_a=TW1_W0(0,0)*h_a(1)+TW1_W0(1,0)*g_a(1); //sin
  double R21_b=TW1_W0(0,0)*h_b(1)+TW1_W0(1,0)*g_b(1);
  double R21_c=TW1_W0(2,0)*k(1);

  double R22_a=TW1_W0(0,1)*h_a(1)+TW1_W0(1,1)*g_a(1); //sin
  double R22_b=TW1_W0(0,1)*h_b(1)+TW1_W0(1,1)*g_b(1);
  double R22_c=TW1_W0(2,1)*k(1);

  double R23_a=TW1_W0(0,2)*h_a(1)+TW1_W0(1,2)*g_a(1); //sin
  double R23_b=TW1_W0(0,2)*h_b(1)+TW1_W0(1,2)*g_b(1);
  double R23_c=TW1_W0(2,2)*k(1);

  double T2_a=TW1_W0(0,3)*h_a(1)+TW1_W0(1,3)*g_a(1); //sin
  double T2_b=TW1_W0(0,3)*h_b(1)+TW1_W0(1,3)*g_b(1);
  double T2_c=TW1_W0(2,3)*k(1)+TC0_C1(1,3);


  double R31_a=TW1_W0(0,0)*h_a(2)+TW1_W0(1,0)*g_a(2); //sin
  double R31_b=TW1_W0(0,0)*h_b(2)+TW1_W0(1,0)*g_b(2);
  double R31_c=TW1_W0(2,0)*k(2);

  double R32_a=TW1_W0(0,1)*h_a(2)+TW1_W0(1,1)*g_a(2); //sin
  double R32_b=TW1_W0(0,1)*h_b(2)+TW1_W0(1,1)*g_b(2);
  double R32_c=TW1_W0(2,1)*k(2);

  double R33_a=TW1_W0(0,2)*h_a(2)+TW1_W0(1,2)*g_a(2); //sin
  double R33_b=TW1_W0(0,2)*h_b(2)+TW1_W0(1,2)*g_b(2);
  double R33_c=TW1_W0(2,2)*k(2);

  double T3_a=TW1_W0(0,3)*h_a(2)+TW1_W0(1,3)*g_a(2); //sin
  double T3_b=TW1_W0(0,3)*h_b(2)+TW1_W0(1,3)*g_b(2);
  double T3_c=TW1_W0(2,3)*k(2)+TC0_C1(2,3);
  
  double ma1=b1/a1;
  double mb1=1/a1;
  Eigen::Matrix3d R_a;
  R_a << R11_a, R12_a, R13_a,
          R21_a, R22_a, R23_a,
          R31_a, R32_a, R33_a;
  Eigen::Matrix3d R_b;
  R_b << R11_b, R12_b, R13_b,
          R21_b, R22_b, R23_b,
          R31_b, R32_b, R33_b;
  Eigen::Matrix3d R_c;
  R_c << R11_c, R12_c, R13_c,
          R21_c, R22_c, R23_c,
          R31_c, R32_c, R33_c;
  Eigen::Matrix3d mR_a = R_C1_C2*R_a;
  Eigen::Matrix3d mR_b = R_C1_C2*R_b;
  Eigen::Matrix3d mR_c = R_C1_C2*R_c;
  double mt1 = R_C1_C2(0,0)+ma1*R_C1_C2(0,1)+mb1*R_C1_C2(0,2);
  double m14 = t_C1_C2(0);
  double mt2 = R_C1_C2(1,0)+ma1*R_C1_C2(1,1)+mb1*R_C1_C2(1,2);
  double m24 = t_C1_C2(1);
  double mt3 = R_C1_C2(2,0)+ma1*R_C1_C2(2,1)+mb1*R_C1_C2(2,2);
  double m34 = t_C1_C2(2);
  double mc = (a2*mt3-mt1)*(a2*m24-b2*m14)+(b2*mt1-a2*mt2)*(a2*m34-m14);
  
  Eigen::Matrix3d mR1,mR2,mR3;
  mR1.block(0,0,1,3)=mR_a.row(1);
  mR1.block(1,0,1,3)=mR_b.row(1);
  mR1.block(2,0,1,3)=mR_c.row(1);
  Eigen::Vector3d m1_abc = mR1*P2;
  
  mR2.block(0,0,1,3)=mR_a.row(0);
  mR2.block(1,0,1,3)=mR_b.row(0);
  mR2.block(2,0,1,3)=mR_c.row(0);
  Eigen::Vector3d n1_abc = mR2*P2;
  
  mR3.block(0,0,1,3)=mR_a.row(2);
  mR3.block(1,0,1,3)=mR_b.row(2);
  mR3.block(2,0,1,3)=mR_c.row(2);
  Eigen::Vector3d m2_abc = mR3*P2;
  
 
  
  
  Eigen::Vector3d coefabc = a2*(a2*mt3-mt1)*m1_abc+a2*(b2*mt1-a2*mt2)*m2_abc+a2*(mt2-b2*mt3)*n1_abc;
  coefabc(2)+=mc;
  
  double coefa=coefabc(0);//sin
  double coefb=coefabc(1);//cos
  double coefc=coefabc(2);
  
//    
  double my_alpha;
  double my_alpha2;
  if (coefa != 0)
  {
    if (coefa*coefa+coefb*coefb-coefc*coefc<0)
    {
//        cout << "can not solve zYaw from this 2 Points\n";
      return solutions;
    }
    double cos_alpha = (-coefb*coefc+abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
    double sin_alpha=-coefc/coefa-coefb/coefa*cos_alpha;
    my_alpha=atan2(sin_alpha,cos_alpha);
    
    double cos_alpha2=(-coefb*coefc-abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
    double sin_alpha2 = -coefc/coefa-coefb/coefa*cos_alpha2;
    my_alpha2=atan2(sin_alpha2,cos_alpha2);
  }
  else
  {
    if (coefb != 0)
    {
      my_alpha = acos(-coefc/coefb);
      my_alpha2 = -my_alpha;
//        cout << "my_alpha = " << my_alpha << endl;
//        cout << "my_alpha2 = " << my_alpha2 << endl;
    }
    else
    {
//        cout << "cannot solve alpha because a==b==0\n";
      return solutions;
    }
  }
  
  std::vector<double> thetas;
  thetas.push_back(my_alpha);
  thetas.push_back(my_alpha2);
//   std::cout << "my_alpha = " << my_alpha << endl;
//   std::cout << "my_alpha2 = " << my_alpha2 << endl;
  
  // for each theta(yaw)
  for (int i=0; i < thetas.size(); ++i)
  {
    double alpha = thetas[i];
    double f_d=Ryx(0,0)*cos(alpha)-Ryx(1,0)*sin(alpha);
    double f_e=Ryx(1,0)*cos(alpha)+Ryx(0,0)*sin(alpha);
    double f_f=Ryx(2,0);
    double m_d=Ryx(0,1)*cos(alpha)-Ryx(1,1)*sin(alpha);
    double m_e=Ryx(1,1)*cos(alpha)+Ryx(0,1)*sin(alpha);
    double m_f=Ryx(2,1);
    double n_d=Ryx(0,2)*cos(alpha)-Ryx(1,2)*sin(alpha);
    double n_e=Ryx(1,2)*cos(alpha)+Ryx(0,2)*sin(alpha);
    double n_f=Ryx(2,2);
    
    
    double T1_d=-TC0_C1(0,0)*f_d-TC0_C1(0,1)*m_d-TC0_C1(0,2)*n_d;
    double T1_e=-TC0_C1(0,0)*f_e-TC0_C1(0,1)*m_e-TC0_C1(0,2)*n_e;
    double T1_f=-TC0_C1(0,0)*f_f-TC0_C1(0,1)*m_f-TC0_C1(0,2)*n_f;
    double T1_g=T1_a*sin(alpha)+T1_b*cos(alpha)+T1_c;
    
    double T2_d=-TC0_C1(1,0)*f_d-TC0_C1(1,1)*m_d-TC0_C1(1,2)*n_d;
    double T2_e=-TC0_C1(1,0)*f_e-TC0_C1(1,1)*m_e-TC0_C1(1,2)*n_e;
    double T2_f=-TC0_C1(1,0)*f_f-TC0_C1(1,1)*m_f-TC0_C1(1,2)*n_f;
    double T2_g=T2_a*sin(alpha)+T2_b*cos(alpha)+T2_c;
    
    double T3_d=-TC0_C1(2,0)*f_d-TC0_C1(2,1)*m_d-TC0_C1(2,2)*n_d;
    double T3_e=-TC0_C1(2,0)*f_e-TC0_C1(2,1)*m_e-TC0_C1(2,2)*n_e;
    double T3_f=-TC0_C1(2,0)*f_f-TC0_C1(2,1)*m_f-TC0_C1(2,2)*n_f;
    double T3_g=T3_a*sin(alpha)+T3_b*cos(alpha)+T3_c;
    
    Eigen::Vector4d T1_defg(T1_d,T1_e,T1_f,T1_g);
    Eigen::Vector4d M14(0,0,0,m14);
    Eigen::Vector4d M24(0,0,0,m24);
    Eigen::Vector4d M34(0,0,0,m34);
    Eigen::Vector4d mT1_defg = mt1*T1_defg+M14;
    Eigen::Vector4d mT2_defg = mt2*T1_defg+M24;
    Eigen::Vector4d mT3_defg = mt3*T1_defg+M34;
    
    
    
    double coefd1=a2*mT2_defg(0)-b2*mT1_defg(0);
	double coefe1=a2*mT2_defg(1)-b2*mT1_defg(1);
	double coeff1=a2*mT2_defg(2)-b2*mT1_defg(2);
    double coefg1=( (mR_b(1,0)*cos(alpha)+mR_a(1,0)*sin(alpha)+mR_c(1,0))*x2
    +(mR_b(1,1)*cos(alpha)+mR_a(1,1)*sin(alpha)+mR_c(1,1))*y2
    +(mR_b(1,2)*cos(alpha)+mR_a(1,2)*sin(alpha)+mR_c(1,2))*z2 )*a2 
    - ( (mR_b(0,0)*cos(alpha)+mR_a(0,0)*sin(alpha)+mR_c(0,0))*x2
    +(mR_b(0,1)*cos(alpha)+mR_a(0,1)*sin(alpha)+mR_c(0,1))*y2
    +(mR_b(0,2)*cos(alpha)+mR_a(0,2)*sin(alpha)+mR_c(0,2))*z2 )*b2 
    + mT2_defg(3)*a2-mT1_defg(3)*b2;
    double coefd2=a1*T2_d-b1*T1_d;
	double coefe2=a1*T2_e-b1*T1_e;
	double coeff2=a1*T2_f-b1*T1_f;
    double coefg2=a1*(T2_g)-b1*(T1_g);
    double coefd3=a1*T3_d-T1_d;double coefe3=a1*T3_e-T1_e;double coeff3=a1*T3_f-T1_f;
    double coefg3=a1*(T3_g)-(T1_g);
    
    Eigen::Matrix3d A;
    A << coefd1, coefe1, coeff1, coefd2, coefe2, coeff2, coefd3, coefe3, coeff3;
    Eigen::Vector3d B;
    B << -coefg1, -coefg2, -coefg3;
//      cout << "Let's check the gt: \n";
//      cout << A*loc_t-B << endl;
    
    if (abs(A.determinant())<0.00001)
    {
//        cout << "cout << A is singular can not solve t\n";
      return solutions;
    }
    Eigen::Vector3d loc_tj = A.inverse()*B; // camera to world
//      cout << "loc_t = " << endl;
//      cout << loc_t << endl;
//      cout << "the estimated loc_tj = \n";
//      cout << loc_tj << endl;
    Eigen::Matrix3d Rz,loc_Rzyx;
    Rz << cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1;
    loc_Rzyx = Rz*Ryx; // camera to world
    Eigen::Matrix4d T_G_C = Eigen::Matrix4d::Identity();
    T_G_C.block(0,0,3,3) = loc_Rzyx;
    T_G_C.block(0,3,3,1) = loc_tj;
    Eigen::Matrix4d T_G_B = T_G_C*T_B_C1.inverse();
    rotation_t rotation = T_G_B.block(0,0,3,3);
    translation_t translation = T_G_B.block(0,3,3,1);
    transformation_t transformation;
    transformation.block<3,3>(0,0) = rotation;
    transformation.col(3) = translation;
    solutions.push_back(transformation);

  }
  
//    cout << "Rs.size()=" << Rs.size() <<  ",ts.size()=" << ts.size() << endl;
  return solutions;
}


opengv::transformations_t 
opengv::absolute_pose::mc2p(
    const AbsoluteLineAdapterBase & adapter,
    const double & mpitch,
    const double & mroll,
	const double & myaw,
    const std::vector<int> & indices )
{
  assert(indices.size()>1);
  return mc2p( adapter, mpitch, mroll, myaw, indices[0], indices[1] );
}

opengv::transformations_t 
opengv::absolute_pose::mc2p(
    const AbsoluteLineAdapterBase & adapter,
    const double & myPitch,
    const double & mxRoll,
	const double & mzYaw,
    size_t index0,
    size_t index1)
{
  transformations_t solutions;
  // 3D point and 2D correspondences idx in mvCorrespondencesIdx
    // the aligned pitch and roll is mxroll mypitch
    // the camera intrinsic parameters is mK (Matrix3d)
  
// 	std::cout << "get into mc2p solution\n";
  
  rotation_t R1 = adapter.getCamRotation(index0); // T_B_C camera in world
  translation_t t1 = adapter.getCamOffset(index0);
  Eigen::Matrix4d T_B_C1 = Eigen::Matrix4d::Identity();
  T_B_C1.block(0,0,3,3) = R1;
  T_B_C1.block(0,3,3,1) = t1;
  
  rotation_t R2 = adapter.getCamRotation(index1); // T_B_C
  translation_t t2 = adapter.getCamOffset(index1);
  Eigen::Matrix4d T_B_C2 = Eigen::Matrix4d::Identity();
  T_B_C2.block(0,0,3,3) = R2;
  T_B_C2.block(0,3,3,1) = t2;
  
  Eigen::Matrix4d T_C1_C2 = T_B_C2.inverse()*T_B_C1; // C1 in C2
  Eigen::Matrix3d R_C1_C2 = T_C1_C2.block<3,3>(0,0);
  Eigen::Vector3d t_C1_C2 = T_C1_C2.block<3,1>(0,3);

  
  Eigen::Vector3d e1 = adapter.getBearingVector(index0);
  Eigen::Vector3d e2 = adapter.getBearingVector(index1);
  Eigen::Vector3d D1 = e1/e1(2);
  Eigen::Vector3d D2 = e2/e2(2);
  Eigen::Vector3d p3d1=adapter.getPoint(index0);
  Eigen::Vector3d p3d2=adapter.getPoint(index1);
  
  if((p3d1-p3d2).norm()<1e-3)
  {
    return solutions;
  }
  
  // C0 coordinates
  double a1=D1[0];double b1=D1[1];
  double a2=D2[0];double b2=D2[1];
  
  if(a1==0)
    return solutions;
  
  // W0 coordinates
  Eigen::Vector3d P10=p3d1;
  Eigen::Vector3d P20=p3d2;
  
  //W1 coordinates
  Eigen::Vector3d P1(0,0,0);
  Eigen::Matrix3d RW0_W1=Eigen::Matrix3d::Identity(3,3);
  Eigen::Vector3d tW0_W1=P1-P10;
  Eigen::Vector3d P2=RW0_W1*P20+tW0_W1;
  double x2=P2[0]; double y2=P2[1]; double z2=P2[2];
  
  Eigen::Matrix4d TW0_W1=Eigen::Matrix4d::Identity(4,4);
  TW0_W1.block(0,3,3,1)=tW0_W1;

  Eigen::Matrix4d TW1_W0=TW0_W1.inverse();
  Eigen::Matrix4d TC0_C1=Eigen::Matrix4d::Identity(4,4);
  Eigen::Matrix3d RC1_C0=Eigen::Matrix3d::Identity(3,3);
  
//   // debug
//   Eigen::Matrix3d Rgt;
//   Rgt <<  0.874865 , 0.396999 , 0.277495,
// -0.272976 , 0.877371, -0.394594,
//  -0.40012,  0.269467,  0.875952;
//  Eigen::Vector3d tgt(  1.36075,-0.422468,1.1324 );
//  Eigen::Matrix4d Tgt_G_B=Eigen::Matrix4d::Identity(4,4);
//  Tgt_G_B.block<3,3>(0,0)=Rgt;
//  Tgt_G_B.block<3,1>(0,3)=tgt;
// //  Eigen::Matrix4d T_G_B = T_G_C*T_B_C1.inverse();
//  Eigen::Matrix4d Tgt_W_C; // c1 in w
//  Tgt_W_C = Tgt_G_B*T_B_C1;
 
 
 Eigen::Matrix3d Rx_Roll;
  Rx_Roll <<  1, 0, 0,
        0,cos(mxRoll),-sin(mxRoll),
        0,sin(mxRoll),cos(mxRoll);
  Eigen::Matrix3d Ry_Pitch;
  Ry_Pitch << cos(myPitch), 0, sin(myPitch),
        0, 1, 0,
        -sin(myPitch),0, cos(myPitch);
  Eigen::Matrix3d Rz_Yaw;
    Rz_Yaw << cos(mzYaw),-sin(mzYaw),0,sin(mzYaw),cos(mzYaw),0,0,0,1;
  Eigen::Matrix3d R_W_B=Rz_Yaw*Ry_Pitch*Rx_Roll;
 
//  Eigen::Matrix3d R_W_B; // restore from eular
 Eigen::Matrix3d Rgt_W_C = R_W_B*R1;
//  Eigen::Matrix3d Rgt_W_C = Tgt_W_C.block<3,3>(0,0);
//  rotationMatrixToEulerAngles()
 Eigen::Vector3d Eulerxyz=rotationMatrixToEulerAngles(Rgt_W_C);
	double mxroll = Eulerxyz(0);
	double mypitch = Eulerxyz(1);
	
	// 	// debug
// 	double mzyaw = Eulerxyz(2);
// 	std::cout << "the C1_W pitch roll and yaw is " << mypitch << "," << mxroll << "," <<mzyaw << endl;
// 	myPitch=mypitch;
// 	mxRoll=mxroll;


//   Eigen::Matrix4d TW1_C1 = TC0_C1*Tgt_W_C.inverse()*TW1_W0;
// //   Eigen::Matrix4d TW1_C2 = T_C1_C2*TW1_C1;
//   
//    std::cout << "let us check the frist eq:\n";
// 	std::cout << a1*TW1_C1(1,3)-b1*TW1_C1(0,3)<< ", "
// 	<< a1*TW1_C1(2,3)-TW1_C1(0,3) << std::endl;
	
	
	
  
  Eigen::Matrix3d Rx;
  Rx <<  1, 0, 0,
        0,cos(mxroll),-sin(mxroll),
        0,sin(mxroll),cos(mxroll);
  Eigen::Matrix3d Ry;
  Ry << cos(mypitch), 0, sin(mypitch),
        0, 1, 0,
        -sin(mypitch),0, cos(mypitch);
  Eigen::Matrix3d Ryx=Ry*Rx;
  
  // R11=coefa*sin(mzYaw)+coefb*cos(mzYaw)+coefc
//    Vector3d h_a=  //1*3
  Eigen::MatrixX3d h_a(1,3),h_b(1,3),g_a(1,3),g_b(1,3),k(1,3);
  h_a = -Ryx.row(1)*RC1_C0;
  h_b = Ryx.row(0)*RC1_C0;
  g_a = h_b;
  g_b = -h_a;
  k = Ryx.row(2)*RC1_C0;
//    cout << "h_a = \n";
//    cout << h_a << endl;
  
  double R11_a=TW1_W0(0,0)*h_a(0)+TW1_W0(1,0)*g_a(0);
  double R11_b=TW1_W0(0,0)*h_b(0)+TW1_W0(1,0)*g_b(0);
  double R11_c=TW1_W0(2,0)*k(0);

  double R12_a=TW1_W0(0,1)*h_a(0)+TW1_W0(1,1)*g_a(0); //sin
  double R12_b=TW1_W0(0,1)*h_b(0)+TW1_W0(1,1)*g_b(0);
  double R12_c=TW1_W0(2,1)*k(0);

  double R13_a=TW1_W0(0,2)*h_a(0)+TW1_W0(1,2)*g_a(0); //sin
  double R13_b=TW1_W0(0,2)*h_b(0)+TW1_W0(1,2)*g_b(0);
  double R13_c=TW1_W0(2,2)*k(0);

  double T1_a=TW1_W0(0,3)*h_a(0)+TW1_W0(1,3)*g_a(0); //sin
  double T1_b=TW1_W0(0,3)*h_b(0)+TW1_W0(1,3)*g_b(0);
  double T1_c=TW1_W0(2,3)*k(0)+TC0_C1(0,3);



  double R21_a=TW1_W0(0,0)*h_a(1)+TW1_W0(1,0)*g_a(1); //sin
  double R21_b=TW1_W0(0,0)*h_b(1)+TW1_W0(1,0)*g_b(1);
  double R21_c=TW1_W0(2,0)*k(1);

  double R22_a=TW1_W0(0,1)*h_a(1)+TW1_W0(1,1)*g_a(1); //sin
  double R22_b=TW1_W0(0,1)*h_b(1)+TW1_W0(1,1)*g_b(1);
  double R22_c=TW1_W0(2,1)*k(1);

  double R23_a=TW1_W0(0,2)*h_a(1)+TW1_W0(1,2)*g_a(1); //sin
  double R23_b=TW1_W0(0,2)*h_b(1)+TW1_W0(1,2)*g_b(1);
  double R23_c=TW1_W0(2,2)*k(1);

  double T2_a=TW1_W0(0,3)*h_a(1)+TW1_W0(1,3)*g_a(1); //sin
  double T2_b=TW1_W0(0,3)*h_b(1)+TW1_W0(1,3)*g_b(1);
  double T2_c=TW1_W0(2,3)*k(1)+TC0_C1(1,3);


  double R31_a=TW1_W0(0,0)*h_a(2)+TW1_W0(1,0)*g_a(2); //sin
  double R31_b=TW1_W0(0,0)*h_b(2)+TW1_W0(1,0)*g_b(2);
  double R31_c=TW1_W0(2,0)*k(2);

  double R32_a=TW1_W0(0,1)*h_a(2)+TW1_W0(1,1)*g_a(2); //sin
  double R32_b=TW1_W0(0,1)*h_b(2)+TW1_W0(1,1)*g_b(2);
  double R32_c=TW1_W0(2,1)*k(2);

  double R33_a=TW1_W0(0,2)*h_a(2)+TW1_W0(1,2)*g_a(2); //sin
  double R33_b=TW1_W0(0,2)*h_b(2)+TW1_W0(1,2)*g_b(2);
  double R33_c=TW1_W0(2,2)*k(2);

  double T3_a=TW1_W0(0,3)*h_a(2)+TW1_W0(1,3)*g_a(2); //sin
  double T3_b=TW1_W0(0,3)*h_b(2)+TW1_W0(1,3)*g_b(2);
  double T3_c=TW1_W0(2,3)*k(2)+TC0_C1(2,3);
  
  double ma1=b1/a1;
  double mb1=1/a1;
  Eigen::Matrix3d R_a;
  R_a << R11_a, R12_a, R13_a,
          R21_a, R22_a, R23_a,
          R31_a, R32_a, R33_a;
  Eigen::Matrix3d R_b;
  R_b << R11_b, R12_b, R13_b,
          R21_b, R22_b, R23_b,
          R31_b, R32_b, R33_b;
  Eigen::Matrix3d R_c;
  R_c << R11_c, R12_c, R13_c,
          R21_c, R22_c, R23_c,
          R31_c, R32_c, R33_c;
  Eigen::Matrix3d mR_a = R_C1_C2*R_a;
  Eigen::Matrix3d mR_b = R_C1_C2*R_b;
  Eigen::Matrix3d mR_c = R_C1_C2*R_c;
  double mt1 = R_C1_C2(0,0)+ma1*R_C1_C2(0,1)+mb1*R_C1_C2(0,2);
  double m14 = t_C1_C2(0);
  double mt2 = R_C1_C2(1,0)+ma1*R_C1_C2(1,1)+mb1*R_C1_C2(1,2);
  double m24 = t_C1_C2(1);
  double mt3 = R_C1_C2(2,0)+ma1*R_C1_C2(2,1)+mb1*R_C1_C2(2,2);
  double m34 = t_C1_C2(2);
  double mc = (a2*mt3-mt1)*(a2*m24-b2*m14)+(b2*mt1-a2*mt2)*(a2*m34-m14);
  
  Eigen::Matrix3d mR1,mR2,mR3;
  mR1.block(0,0,1,3)=mR_a.row(1);
  mR1.block(1,0,1,3)=mR_b.row(1);
  mR1.block(2,0,1,3)=mR_c.row(1);
  Eigen::Vector3d m1_abc = mR1*P2;
  
  mR2.block(0,0,1,3)=mR_a.row(0);
  mR2.block(1,0,1,3)=mR_b.row(0);
  mR2.block(2,0,1,3)=mR_c.row(0);
  Eigen::Vector3d n1_abc = mR2*P2;
  
  mR3.block(0,0,1,3)=mR_a.row(2);
  mR3.block(1,0,1,3)=mR_b.row(2);
  mR3.block(2,0,1,3)=mR_c.row(2);
  Eigen::Vector3d m2_abc = mR3*P2;
  
 
  
  
  Eigen::Vector3d coefabc = a2*(a2*mt3-mt1)*m1_abc+a2*(b2*mt1-a2*mt2)*m2_abc+a2*(mt2-b2*mt3)*n1_abc;
  coefabc(2)+=mc;
  
  double coefa=coefabc(0);//sin
  double coefb=coefabc(1);//cos
  double coefc=coefabc(2);
  
//    
  double my_alpha;
  double my_alpha2;
  if (coefa != 0)
  {
    if (coefa*coefa+coefb*coefb-coefc*coefc<0)
    {
//        cout << "can not solve zYaw from this 2 Points\n";
      return solutions;
    }
    double cos_alpha = (-coefb*coefc+abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
    double sin_alpha=-coefc/coefa-coefb/coefa*cos_alpha;
    my_alpha=atan2(sin_alpha,cos_alpha);
    
    double cos_alpha2=(-coefb*coefc-abs(coefa)*sqrt(coefa*coefa+coefb*coefb-coefc*coefc))/(coefa*coefa+coefb*coefb);
    double sin_alpha2 = -coefc/coefa-coefb/coefa*cos_alpha2;
    my_alpha2=atan2(sin_alpha2,cos_alpha2);
  }
  else
  {
    if (coefb != 0)
    {
      my_alpha = acos(-coefc/coefb);
      my_alpha2 = -my_alpha;
//        cout << "my_alpha = " << my_alpha << endl;
//        cout << "my_alpha2 = " << my_alpha2 << endl;
    }
    else
    {
//        cout << "cannot solve alpha because a==b==0\n";
      return solutions;
    }
  }
  
  std::vector<double> thetas;
  thetas.push_back(my_alpha);
  thetas.push_back(my_alpha2);
//   std::cout << "my_alpha = " << my_alpha << endl;
//   std::cout << "my_alpha2 = " << my_alpha2 << endl;
  
  // for each theta(yaw)
  for (int i=0; i < thetas.size(); ++i)
  {
    double alpha = thetas[i];
    double f_d=Ryx(0,0)*cos(alpha)-Ryx(1,0)*sin(alpha);
    double f_e=Ryx(1,0)*cos(alpha)+Ryx(0,0)*sin(alpha);
    double f_f=Ryx(2,0);
    double m_d=Ryx(0,1)*cos(alpha)-Ryx(1,1)*sin(alpha);
    double m_e=Ryx(1,1)*cos(alpha)+Ryx(0,1)*sin(alpha);
    double m_f=Ryx(2,1);
    double n_d=Ryx(0,2)*cos(alpha)-Ryx(1,2)*sin(alpha);
    double n_e=Ryx(1,2)*cos(alpha)+Ryx(0,2)*sin(alpha);
    double n_f=Ryx(2,2);
    
    
    double T1_d=-TC0_C1(0,0)*f_d-TC0_C1(0,1)*m_d-TC0_C1(0,2)*n_d;
    double T1_e=-TC0_C1(0,0)*f_e-TC0_C1(0,1)*m_e-TC0_C1(0,2)*n_e;
    double T1_f=-TC0_C1(0,0)*f_f-TC0_C1(0,1)*m_f-TC0_C1(0,2)*n_f;
    double T1_g=T1_a*sin(alpha)+T1_b*cos(alpha)+T1_c;
    
    double T2_d=-TC0_C1(1,0)*f_d-TC0_C1(1,1)*m_d-TC0_C1(1,2)*n_d;
    double T2_e=-TC0_C1(1,0)*f_e-TC0_C1(1,1)*m_e-TC0_C1(1,2)*n_e;
    double T2_f=-TC0_C1(1,0)*f_f-TC0_C1(1,1)*m_f-TC0_C1(1,2)*n_f;
    double T2_g=T2_a*sin(alpha)+T2_b*cos(alpha)+T2_c;
    
    double T3_d=-TC0_C1(2,0)*f_d-TC0_C1(2,1)*m_d-TC0_C1(2,2)*n_d;
    double T3_e=-TC0_C1(2,0)*f_e-TC0_C1(2,1)*m_e-TC0_C1(2,2)*n_e;
    double T3_f=-TC0_C1(2,0)*f_f-TC0_C1(2,1)*m_f-TC0_C1(2,2)*n_f;
    double T3_g=T3_a*sin(alpha)+T3_b*cos(alpha)+T3_c;
    
    Eigen::Vector4d T1_defg(T1_d,T1_e,T1_f,T1_g);
    Eigen::Vector4d M14(0,0,0,m14);
    Eigen::Vector4d M24(0,0,0,m24);
    Eigen::Vector4d M34(0,0,0,m34);
    Eigen::Vector4d mT1_defg = mt1*T1_defg+M14;
    Eigen::Vector4d mT2_defg = mt2*T1_defg+M24;
    Eigen::Vector4d mT3_defg = mt3*T1_defg+M34;
    
    
    
    double coefd1=a2*mT2_defg(0)-b2*mT1_defg(0);
	double coefe1=a2*mT2_defg(1)-b2*mT1_defg(1);
	double coeff1=a2*mT2_defg(2)-b2*mT1_defg(2);
    double coefg1=( (mR_b(1,0)*cos(alpha)+mR_a(1,0)*sin(alpha)+mR_c(1,0))*x2
    +(mR_b(1,1)*cos(alpha)+mR_a(1,1)*sin(alpha)+mR_c(1,1))*y2
    +(mR_b(1,2)*cos(alpha)+mR_a(1,2)*sin(alpha)+mR_c(1,2))*z2 )*a2 
    - ( (mR_b(0,0)*cos(alpha)+mR_a(0,0)*sin(alpha)+mR_c(0,0))*x2
    +(mR_b(0,1)*cos(alpha)+mR_a(0,1)*sin(alpha)+mR_c(0,1))*y2
    +(mR_b(0,2)*cos(alpha)+mR_a(0,2)*sin(alpha)+mR_c(0,2))*z2 )*b2 
    + mT2_defg(3)*a2-mT1_defg(3)*b2;
    double coefd2=a1*T2_d-b1*T1_d;
	double coefe2=a1*T2_e-b1*T1_e;
	double coeff2=a1*T2_f-b1*T1_f;
    double coefg2=a1*(T2_g)-b1*(T1_g);
    double coefd3=a1*T3_d-T1_d;double coefe3=a1*T3_e-T1_e;double coeff3=a1*T3_f-T1_f;
    double coefg3=a1*(T3_g)-(T1_g);
    
    Eigen::Matrix3d A;
    A << coefd1, coefe1, coeff1, coefd2, coefe2, coeff2, coefd3, coefe3, coeff3;
    Eigen::Vector3d B;
    B << -coefg1, -coefg2, -coefg3;
//      cout << "Let's check the gt: \n";
//      cout << A*loc_t-B << endl;
    
    if (abs(A.determinant())<0.00001)
    {
//        cout << "cout << A is singular can not solve t\n";
      return solutions;
    }
    Eigen::Vector3d loc_tj = A.inverse()*B; // camera to world
//      cout << "loc_t = " << endl;
//      cout << loc_t << endl;
//      cout << "the estimated loc_tj = \n";
//      cout << loc_tj << endl;
    Eigen::Matrix3d Rz,loc_Rzyx;
    Rz << cos(alpha),-sin(alpha),0,sin(alpha),cos(alpha),0,0,0,1;
    loc_Rzyx = Rz*Ryx; // camera to world
    Eigen::Matrix4d T_G_C = Eigen::Matrix4d::Identity();
    T_G_C.block(0,0,3,3) = loc_Rzyx;
    T_G_C.block(0,3,3,1) = loc_tj;
    Eigen::Matrix4d T_G_B = T_G_C*T_B_C1.inverse();
    rotation_t rotation = T_G_B.block(0,0,3,3);
    translation_t translation = T_G_B.block(0,3,3,1);
    transformation_t transformation;
    transformation.block<3,3>(0,0) = rotation;
    transformation.col(3) = translation;
    solutions.push_back(transformation);

  }
  
//    cout << "Rs.size()=" << Rs.size() <<  ",ts.size()=" << ts.size() << endl;
  return solutions;
}



opengv::transformations_t
opengv::absolute_pose::p3p_kneip(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  assert(indices.size()>2);
  return p3p_kneip( adapter, indices[0], indices[1], indices[2] );
}

opengv::transformations_t
opengv::absolute_pose::p3p_kneip(
    const AbsoluteAdapterBase & adapter,
    size_t index0,
    size_t index1,
    size_t index2)
{
  bearingVectors_t f;
  f.push_back(adapter.getBearingVector(index0));
  f.push_back(adapter.getBearingVector(index1));
  f.push_back(adapter.getBearingVector(index2));
  points_t p;
  p.push_back(adapter.getPoint(index0));
  p.push_back(adapter.getPoint(index1));
  p.push_back(adapter.getPoint(index2));
  transformations_t solutions;
  modules::p3p_kneip_main( f, p, solutions );
  return solutions;
}

opengv::transformations_t
opengv::absolute_pose::p3p_gao(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  assert(indices.size()>2);
  return p3p_gao( adapter, indices[0], indices[1], indices[2] );
}

opengv::transformations_t
opengv::absolute_pose::p3p_gao(
    const AbsoluteAdapterBase & adapter,
    size_t index0,
    size_t index1,
    size_t index2)
{
  bearingVectors_t f;
  f.push_back(adapter.getBearingVector(index0));
  f.push_back(adapter.getBearingVector(index1));
  f.push_back(adapter.getBearingVector(index2));
  points_t p;
  p.push_back(adapter.getPoint(index0));
  p.push_back(adapter.getPoint(index1));
  p.push_back(adapter.getPoint(index2));
  transformations_t solutions;
  modules::p3p_gao_main( f, p, solutions );
  return solutions;
}

opengv::transformations_t
opengv::absolute_pose::gp3p(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  assert(indices.size()>2);

  Eigen::Matrix3d f;
  Eigen::Matrix3d v;
  Eigen::Matrix3d p;

  for(size_t i = 0; i < 3; i++)
  {
    f.col(i) = adapter.getBearingVector(indices[i]);
    rotation_t R = adapter.getCamRotation(indices[i]);
// 	std::cout<<"========================"<<endl;
// 	std::cout<<R<<endl;
    
    //unrotate the bearingVectors already so the camera rotation doesn't appear
    //in the problem
    f.col(i) = R * f.col(i);
    v.col(i) = adapter.getCamOffset(indices[i]);
    p.col(i) = adapter.getPoint(indices[i]);
  }

  transformations_t solutions;
  modules::gp3p_main(f,v,p,solutions);

  return solutions;
}

opengv::transformations_t
opengv::absolute_pose::gp3p(
    const AbsoluteAdapterBase & adapter,
    size_t index0,
    size_t index1,
    size_t index2)
{
  std::vector<int> indices;
  indices.push_back(index0);
  indices.push_back(index1);
  indices.push_back(index2);

  return gp3p(adapter,indices);
}

namespace opengv
{
namespace absolute_pose
{

transformation_t epnp(
    const AbsoluteAdapterBase & adapter,
    const Indices & indices )
{
  //starting from 4 points, we have a unique solution
  assert(indices.size() > 3);

  modules::Epnp PnP;
  PnP.set_maximum_number_of_correspondences(indices.size());
  PnP.reset_correspondences();

  for( size_t i = 0; i < indices.size(); i++ )
  {
    point_t p = adapter.getPoint(indices[i]);
    bearingVector_t f = adapter.getBearingVector(indices[i]);
    PnP.add_correspondence(p[0], p[1], p[2], f[0], f[1], f[2]);
  }

  double R_epnp[3][3], t_epnp[3];
  PnP.compute_pose(R_epnp, t_epnp);

  rotation_t rotation;
  translation_t translation;

  for(int r = 0; r < 3; r++)
  {
    for(int c = 0; c < 3; c++)
      rotation(r,c) = R_epnp[r][c];
  }

  translation[0] = t_epnp[0];
  translation[1] = t_epnp[1];
  translation[2] = t_epnp[2];

  //take inverse transformation
  rotation.transposeInPlace();
  translation = -rotation * translation;

  transformation_t transformation;
  transformation.col(3) = translation;
  transformation.block<3,3>(0,0) = rotation;
  return transformation;
}

}
}

opengv::transformation_t
opengv::absolute_pose::epnp( const AbsoluteAdapterBase & adapter )
{
  Indices idx(adapter.getNumberCorrespondences());
  return epnp(adapter,idx);
}

opengv::transformation_t
opengv::absolute_pose::epnp(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  Indices idx(indices);
  return epnp(adapter,idx);
}

namespace opengv
{
namespace absolute_pose
{

transformation_t gpnp(
    const AbsoluteAdapterBase & adapter,
    const Indices & indices )
{
  assert( indices.size() > 5 );

  //compute the centroid
  point_t c0 = Eigen::Vector3d::Zero();
  for( size_t i = 0; i < indices.size(); i++ )
    c0 = c0 + adapter.getPoint(indices[i]);
  c0 = c0 / indices.size();

  //compute the point-cloud
  Eigen::MatrixXd p(3,indices.size());
  for( size_t i = 0; i < indices.size(); i++ )
    p.col(i) = adapter.getPoint(indices[i]) - c0;

  //compute the moment
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(
      p,
      Eigen::ComputeThinU | Eigen::ComputeThinV );

  //define the control points
  points_t c;
  c.push_back(c0);
  //c.push_back(c0 + SVD.singularValues()[0] * SVD.matrixU().col(0));
  //c.push_back(c0 + SVD.singularValues()[1] * SVD.matrixU().col(1));
  //c.push_back(c0 + SVD.singularValues()[2] * SVD.matrixU().col(2));
  c.push_back(c0 + 15.0 * SVD.matrixU().col(0));
  c.push_back(c0 + 15.0 * SVD.matrixU().col(1));
  c.push_back(c0 + 15.0 * SVD.matrixU().col(2));

  //derive the barycentric frame
  Eigen::Vector3d e1 = c[1]-c0;
  double e1dote1 = e1.dot(e1);
  Eigen::Vector3d e2 = c[2]-c0;
  double e2dote2 = e2.dot(e2);
  Eigen::Vector3d e3 = c[3]-c0;
  double e3dote3 = e3.dot(e3);

  //derive the weighting factors
  Eigen::MatrixXd weights(4,indices.size());
  for( size_t i = 0; i < indices.size(); i++ )
  {
    Eigen::Vector3d temp = p.col(i);
    weights(1,i) = temp.dot(e1)/e1dote1;
    weights(2,i) = temp.dot(e2)/e2dote2;
    weights(3,i) = temp.dot(e3)/e3dote3;
    weights(0,i) = 1.0-(weights(1,i)+weights(2,i)+weights(3,i));
  }

  //setup matrix A and vector b
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*indices.size(),12);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(2*indices.size(),1);
  for( size_t i = 0; i < indices.size(); i++ )
  {
    translation_t camOffset = adapter.getCamOffset(indices[i]);
    rotation_t camRotation = adapter.getCamRotation(indices[i]);
    //respect the rotation
    bearingVector_t f = camRotation * adapter.getBearingVector(indices[i]);

    A(2*i,0)  =  weights(0,i)*f[2];
    A(2*i,2)  = -weights(0,i)*f[0];
    A(2*i,3)  =  weights(1,i)*f[2];
    A(2*i,5)  = -weights(1,i)*f[0];
    A(2*i,6)  =  weights(2,i)*f[2];
    A(2*i,8)  = -weights(2,i)*f[0];
    A(2*i,9)  =  weights(3,i)*f[2];
    A(2*i,11) = -weights(3,i)*f[0];

    A(2*i+1,1)  =  weights(0,i)*f[2];
    A(2*i+1,2)  = -weights(0,i)*f[1];
    A(2*i+1,4)  =  weights(1,i)*f[2];
    A(2*i+1,5)  = -weights(1,i)*f[1];
    A(2*i+1,7)  =  weights(2,i)*f[2];
    A(2*i+1,8)  = -weights(2,i)*f[1];
    A(2*i+1,10) =  weights(3,i)*f[2];
    A(2*i+1,11) = -weights(3,i)*f[1];

    b(2*i,0)   = f[2]*camOffset[0]-f[0]*camOffset[2];
    b(2*i+1,0) = f[2]*camOffset[1]-f[1]*camOffset[2];
  }

  //computing the SVD
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD2(
      A,
      Eigen::ComputeThinV | Eigen::ComputeThinU );

  //computing the pseudoinverse
  Eigen::MatrixXd invD = Eigen::MatrixXd::Zero(12,12);
  Eigen::MatrixXd D = SVD2.singularValues();
  for( size_t i = 0; i < 12; i++ )
  {
    if( D(i,0) > 1.e-6 )
      invD(i,i) = 1.0/D(i,0);
    else
      invD(i,i) = 0.0;
  }

  //Extract the nullsapce vectors;
  Eigen::MatrixXd V = SVD2.matrixV();

  //computing the nullspace intercept
  Eigen::MatrixXd pinvA = V * invD * SVD2.matrixU().transpose();

  //compute the intercept
  Eigen::Matrix<double,12,1> a = pinvA * b;

  //compute the solution
  transformation_t transformation;
  modules::gpnp_main( a, V, c, transformation );
  return transformation;
}

}
}

opengv::transformation_t
opengv::absolute_pose::gpnp( const AbsoluteAdapterBase & adapter )
{
  Indices idx(adapter.getNumberCorrespondences());
  return gpnp(adapter,idx);
}

opengv::transformation_t
opengv::absolute_pose::gpnp(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  Indices idx(indices);
  return gpnp(adapter,idx);
}

namespace opengv
{
namespace absolute_pose
{

void fill3x10( const Eigen::Vector3d & x, Eigen::Matrix<double,3,10> & Phi )
{
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];
  
  Phi << x1,  x1, -x1, -x1,     0.0,  2.0*x3, -2.0*x2, 2.0*x2, 2.0*x3,    0.0,
         x2, -x2,  x2, -x2, -2.0*x3,     0.0,  2.0*x1, 2.0*x1,    0.0, 2.0*x3,
         x3, -x3, -x3,  x3,  2.0*x2, -2.0*x1,     0.0,    0.0, 2.0*x1, 2.0*x2;
}

void f(
    const Eigen::Matrix<double,10,10> & M,
    const Eigen::Matrix<double,1,10> & C,
    double gamma,
    Eigen::Vector3d & f )
{
  f[0] = (2*M(0,4)+2*C(0,4));
  f[1] = (2*M(0,5)+2*C(0,5));
  f[2] = (2*M(0,6)+2*C(0,6));
}

void Jac(
    const Eigen::Matrix<double,10,10> & M,
    const Eigen::Matrix<double,1,10> & C,
    double gamma,
    Eigen::Matrix3d & Jac )
{
  Jac(0,0) = (2*M(4,4)+4*M(0,1)-4*M(0,0)+4*C(0,1)-4*C(0,0));
  Jac(0,1) = (2*M(5,4)+2*M(0,7)+2*C(0,7));
  Jac(0,2) = (2*M(6,4)+2*M(0,8)+2*C(0,8));
  Jac(1,0) = (2*M(4,5)+2*M(0,7)+2*C(0,7));
  Jac(1,1) = (2*M(5,5)+4*M(0,2)-4*M(0,0)+4*C(0,2)-4*C(0,0));
  Jac(1,2) = (2*M(6,5)+2*M(0,9)+2*C(0,9));
  Jac(2,0) = (2*M(4,6)+2*M(0,8)+2*C(0,8));
  Jac(2,1) = (2*M(5,6)+2*M(0,9)+2*C(0,9));
  Jac(2,2) = (2*M(6,6)+4*M(0,3)-4*M(0,0)+4*C(0,3)-4*C(0,0));
}



transformations_t upnp(
    const AbsoluteAdapterBase & adapter,
    const Indices & indices )
{
  assert( indices.size() > 2 );
    
  Eigen::Matrix<double,3,3> F = Eigen::Matrix3d::Zero();
  for( int i = 0; i < (int) indices.size(); i++ )
  {
    Eigen::Matrix<double,3,1> f = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    F += f * f.transpose();
  }
  
  Eigen::Matrix<double,3,3> H_inv = (indices.size() * Eigen::Matrix<double,3,3>::Identity()) - F;
  Eigen::Matrix<double,3,3> H = H_inv.inverse();
  
  Eigen::Matrix<double,3,10> I = Eigen::Matrix<double,3,10>::Zero();
  Eigen::Matrix<double,3,1> J = Eigen::Matrix<double,3,1>::Zero();
  Eigen::Matrix<double,3,10> Phi;
  
  for( int i = 0; i < (int) indices.size(); i++ )
  {
    Eigen::Matrix<double,3,1> f = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    Eigen::Matrix<double,3,3> Vk = H * ( f * f.transpose() - Eigen::Matrix<double,3,3>::Identity() );
    Eigen::Matrix<double,3,1> p = adapter.getPoint(indices[i]);
    Eigen::Matrix<double,3,1> v = adapter.getCamOffset(indices[i]);
    
    fill3x10(p,Phi);
    I += Vk * Phi;
    J += Vk * v;
  }
  
  Eigen::Matrix<double,10,10> M = Eigen::Matrix<double,10,10>::Zero();
  Eigen::Matrix<double,1,10>  C = Eigen::Matrix<double,1,10>::Zero();
  double gamma = 0.0;
  
  for(int i = 0; i < (int) indices.size(); i++ )
  {    
    Eigen::Matrix<double,3,1> f = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
    Eigen::Matrix<double,3,1> v = adapter.getCamOffset(indices[i]);
    Eigen::Matrix<double,3,1> p = adapter.getPoint(indices[i]);
    
    fill3x10(p,Phi);
    Eigen::Matrix<double,3,3> temp = f*f.transpose() - Eigen::Matrix<double,3,3>::Identity();
    Eigen::Matrix<double,3,10> Ai =  temp * (Phi + I);
    Eigen::Matrix<double,3, 1> bi = -temp * (  v + J);
    
    M     += (Ai.transpose() * Ai);
    C     += (bi.transpose() * Ai);
    gamma += (bi.transpose() * bi);
  }
  
  //now do the main computation
  std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > > quaternions1;
  if( indices.size() > 4 )
    modules::upnp_main_sym( M, C, gamma, quaternions1 );
  else
    modules::upnp_main( M, C, gamma, quaternions1 );
  
  //prepare the output vector
  transformations_t transformations;
  
  //Round 1: chirality check
  std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > > quaternions2;
  for( size_t i = 0; i < quaternions1.size(); i++ )
  {
    rotation_t Rinv = math::quaternion2rot(quaternions1[i].second);
    
    Eigen::Matrix<double,10,1> s;
    modules::upnp_fill_s( quaternions1[i].second, s );
    translation_t tinv = I*s - J;
    
    if( transformations.size() == 0 )
    {
      transformation_t newTransformation;
      newTransformation.block<3,3>(0,0) = Rinv.transpose();
      newTransformation.block<3,1>(0,3) = -newTransformation.block<3,3>(0,0) * tinv;
      transformations.push_back(newTransformation);
    }
    
    int count_negative = 0;
    
    for( int j = 0; j < (int) indices.size(); j++ )
    {
      Eigen::Matrix<double,3,1> f = adapter.getCamRotation(indices[j]) * adapter.getBearingVector(indices[j]);
      Eigen::Matrix<double,3,1> p = adapter.getPoint(indices[j]);
      Eigen::Matrix<double,3,1> v = adapter.getCamOffset(indices[j]);
      
      Eigen::Vector3d p_est = Rinv*p + tinv - v;
      
      if( p_est.transpose()*f < 0.0 )
        count_negative++;
    }
    
    if( count_negative < floor(0.2 * indices.size() + 0.5) )
      quaternions2.push_back(quaternions1[i]);
  }
  
  if( quaternions2.size() == 0 )
    return transformations;
  else
    transformations.clear();
  
  //Round 2: Second order optimality (plus polishing)
  Eigen::Matrix<double,3,10> I_cay;
  Eigen::Matrix<double,10,10> M_cay;
  Eigen::Matrix<double,1,10>  C_cay;
  double gamma_cay;
  
  for( size_t q = 0; q < quaternions2.size(); q++ )
  {    
    I_cay = Eigen::Matrix<double,3,10>::Zero();
    rotation_t Rinv = math::quaternion2rot(quaternions2[q].second);
    
    for( int i = 0; i < (int) indices.size(); i++ )
    {
      Eigen::Matrix<double,3,1> f = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
      Eigen::Matrix<double,3,3> Vk = H * ( f * f.transpose() - Eigen::Matrix<double,3,3>::Identity() );
      Eigen::Matrix<double,3,1> p = Rinv * adapter.getPoint(indices[i]);
      
      fill3x10(p,Phi);
      I_cay += Vk * Phi;
    }
    
    M_cay = Eigen::Matrix<double,10,10>::Zero();
    C_cay = Eigen::Matrix<double,1,10>::Zero();
    gamma_cay = 0.0;
    
    for(int i = 0; i < (int) indices.size(); i++ )
    {    
      Eigen::Matrix<double,3,1> f = adapter.getCamRotation(indices[i]) * adapter.getBearingVector(indices[i]);
      Eigen::Matrix<double,3,1> v = adapter.getCamOffset(indices[i]);
      Eigen::Matrix<double,3,1> p = Rinv * adapter.getPoint(indices[i]);
      
      fill3x10(p,Phi);
      Eigen::Matrix<double,3,3> temp = f*f.transpose() - Eigen::Matrix<double,3,3>::Identity();
      Eigen::Matrix<double,3,10> Ai =  temp * (Phi + I_cay);
      Eigen::Matrix<double,3,1> bi = -temp * (  v + J);
      
      M_cay     += (Ai.transpose() * Ai);
      C_cay     += (bi.transpose() * Ai);
      gamma_cay += (bi.transpose() * bi);
    }
    
    //now analyze the eigenvalues of the "Hessian"
    Eigen::Vector3d val;
    Eigen::Matrix3d Jacobian;
    f( M_cay, C_cay, gamma_cay, val );
    Jac( M_cay, C_cay, gamma_cay, Jacobian );
    std::vector<double> characteristicPolynomial;
    characteristicPolynomial.push_back(-1.0);
    characteristicPolynomial.push_back(Jacobian(2,2)+Jacobian(1,1)+Jacobian(0,0));
    characteristicPolynomial.push_back(-Jacobian(2,2)*Jacobian(1,1)-Jacobian(2,2)*Jacobian(0,0)-Jacobian(1,1)*Jacobian(0,0)+pow(Jacobian(1,2),2)+pow(Jacobian(0,2),2)+pow(Jacobian(0,1),2));
    characteristicPolynomial.push_back(Jacobian(2,2)*Jacobian(1,1)*Jacobian(0,0)+2*Jacobian(1,2)*Jacobian(0,2)*Jacobian(0,1)-Jacobian(2,2)*pow(Jacobian(0,1),2)-pow(Jacobian(1,2),2)*Jacobian(0,0)-Jacobian(1,1)*pow(Jacobian(0,2),2));

    /* This is commented in the upstream repository.
    std::vector<double> roots = opengv::math::o3_roots( characteristicPolynomial );

    bool allPositive = true;
    for( size_t i = 0; i < roots.size(); i++ )
    {
      if( roots[i] < 0.0 )
      {
        allPositive = false;
        break;
      }
    }

    if( allPositive )
    */
    {
      //perform the polishing step
      Eigen::Vector3d cay = - Jacobian.inverse() * val;
      rotation_t Rinv2 = math::cayley2rot(cay) * Rinv;
      quaternion_t q = math::rot2quaternion(Rinv2);
      
      Eigen::Matrix<double,10,1> s;
      modules::upnp_fill_s(q,s);
      translation_t tinv = I*s - J;
      
      transformation_t newTransformation;
      newTransformation.block<3,3>(0,0) = Rinv2.transpose();
      newTransformation.block<3,1>(0,3) = -newTransformation.block<3,3>(0,0) * tinv;
      transformations.push_back(newTransformation);
    }
  }
  
  //if there are no results, simply add the one with lowest score
  if( transformations.size() == 0 )
  {
    Eigen::Vector4d q = quaternions2[0].second;
    Eigen::Matrix<double,10,1> s;
    modules::upnp_fill_s(q,s);
    translation_t tinv = I*s - J;
    rotation_t Rinv = math::quaternion2rot(q);
    
    transformation_t newTransformation;
    newTransformation.block<3,3>(0,0) = Rinv.transpose();
    newTransformation.block<3,1>(0,3) = -newTransformation.block<3,3>(0,0) * tinv;
    transformations.push_back(newTransformation);
  }
  
  return transformations;
}

}
}

opengv::transformations_t
opengv::absolute_pose::upnp( const AbsoluteAdapterBase & adapter )
{
  Indices idx(adapter.getNumberCorrespondences());
  return upnp(adapter,idx);
}

opengv::transformations_t
opengv::absolute_pose::upnp(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  Indices idx(indices);
  return upnp(adapter,idx);
}

opengv::transformations_t
opengv::absolute_pose::upnp_ransac(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
	std::vector<int>  new_indices;
// 	for(int i=0;i < indices.size();++i)
// 		new_indices.push_back(indices[i]);
	new_indices.push_back(indices[0]);
	new_indices.push_back(indices[1]);
	new_indices.push_back(indices[2]);
  Indices idx(new_indices);
  return upnp(adapter,idx);
}

namespace opengv
{
namespace absolute_pose
{

struct OptimizeNonlinearFunctor1 : OptimizationFunctor<double>
{
  const AbsoluteAdapterBase & _adapter;
  const Indices & _indices;

  OptimizeNonlinearFunctor1(
      const AbsoluteAdapterBase & adapter,
      const Indices & indices ) :
      OptimizationFunctor<double>(6,indices.size()),
      _adapter(adapter),
      _indices(indices) {}

  int operator()(const VectorXd &x, VectorXd &fvec) const
  {
    assert( x.size() == 6 );
    assert( (unsigned int) fvec.size() == _indices.size());

    //compute the current position
    translation_t translation = x.block<3,1>(0,0);
    cayley_t cayley = x.block<3,1>(3,0);
    rotation_t rotation = math::cayley2rot(cayley);

    //compute inverse transformation
    transformation_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = rotation.transpose();
    inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*translation;

    Eigen::Matrix<double,4,1> p_hom;
    p_hom[3] = 1.0;

    for(size_t i = 0; i < _indices.size(); i++)
    {
      //get point in homogeneous form
      p_hom.block<3,1>(0,0) = _adapter.getPoint(_indices[i]);

      //compute the reprojection (this is working for both central and
      //non-central case)
      point_t bodyReprojection = inverseSolution * p_hom;
      point_t reprojection = _adapter.getCamRotation(_indices[i]).transpose() *
          (bodyReprojection - _adapter.getCamOffset(_indices[i]));
      reprojection = reprojection / reprojection.norm();

      //compute the score
      double factor = 1.0;
      fvec[i] = factor *
          (1.0 -
          (reprojection.transpose() * _adapter.getBearingVector(_indices[i])));
    }

    return 0;
  }
};

transformation_t optimize_nonlinear(
    const AbsoluteAdapterBase & adapter,
    const Indices & indices )
{
  const int n=6;
  VectorXd x(n);

  x.block<3,1>(0,0) = adapter.gett();
  x.block<3,1>(3,0) = math::rot2cayley(adapter.getR());

  OptimizeNonlinearFunctor1 functor( adapter, indices );
  NumericalDiff<OptimizeNonlinearFunctor1> numDiff(functor);
  LevenbergMarquardt< NumericalDiff<OptimizeNonlinearFunctor1> > lm(numDiff);

  lm.resetParameters();
  lm.parameters.ftol = 1.E1*NumTraits<double>::epsilon();
  lm.parameters.xtol = 1.E1*NumTraits<double>::epsilon();
  lm.parameters.maxfev = 1000;
  lm.minimize(x);

  transformation_t transformation;
  transformation.col(3) = x.block<3,1>(0,0);
  transformation.block<3,3>(0,0) = math::cayley2rot(x.block<3,1>(3,0));
  return transformation;
}

}
}

opengv::transformation_t
opengv::absolute_pose::optimize_nonlinear( const AbsoluteAdapterBase & adapter )
{
  Indices idx(adapter.getNumberCorrespondences());
  return optimize_nonlinear(adapter,idx);
}

opengv::transformation_t
opengv::absolute_pose::optimize_nonlinear(
    const AbsoluteAdapterBase & adapter,
    const std::vector<int> & indices )
{
  Indices idx(indices);
  return optimize_nonlinear(adapter,idx);
}

// optimize with point or/and with line

namespace opengv
{
namespace absolute_pose
{

struct OptimizeNonlinearFunctor2 : OptimizationFunctor<double>
{
  const AbsoluteLineAdapterBase & _adapter;
  const Indices & _indices;

  OptimizeNonlinearFunctor2(
      const AbsoluteLineAdapterBase & adapter,
      const Indices & indices ) :
      OptimizationFunctor<double>(6,indices.size()),
      _adapter(adapter),
      _indices(indices) {}

  int operator()(const VectorXd &x, VectorXd &fvec) const
  {
    assert( x.size() == 6 );
    assert( (unsigned int) fvec.size() == _indices.size());

    //compute the current position
    translation_t translation = x.block<3,1>(0,0);
    cayley_t cayley = x.block<3,1>(3,0);
    rotation_t rotation = math::cayley2rot(cayley);

    //compute inverse transformation
    transformation_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = rotation.transpose();
    inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*translation;

    Eigen::Matrix<double,4,1> p_hom;
    p_hom[3] = 1.0;

    for(size_t i = 0; i < _indices.size(); i++)
    {
      //get point in homogeneous form
      p_hom.block<3,1>(0,0) = _adapter.getPoint(_indices[i]);

      //compute the reprojection (this is working for both central and
      //non-central case)
      point_t bodyReprojection = inverseSolution * p_hom;
      point_t reprojection = _adapter.getCamRotation(_indices[i]).transpose() *
          (bodyReprojection - _adapter.getCamOffset(_indices[i]));
      reprojection = reprojection / reprojection.norm();

      //compute the score
      double factor = 1.0;
      fvec[i] = factor *
          (1.0 -
          (reprojection.transpose() * _adapter.getBearingVector(_indices[i])));
    }

    return 0;
  }
};

transformation_t optimize_nonlinear(
    const AbsoluteLineAdapterBase & adapter,
    const Indices & indices )
{
  const int n=6;
  VectorXd x(n);

  x.block<3,1>(0,0) = adapter.gett();
  x.block<3,1>(3,0) = math::rot2cayley(adapter.getR());

  OptimizeNonlinearFunctor2 functor( adapter, indices );
  NumericalDiff<OptimizeNonlinearFunctor2> numDiff(functor);
  LevenbergMarquardt< NumericalDiff<OptimizeNonlinearFunctor2> > lm(numDiff);

  lm.resetParameters();
  lm.parameters.ftol = 1.E1*NumTraits<double>::epsilon();
  lm.parameters.xtol = 1.E1*NumTraits<double>::epsilon();
  lm.parameters.maxfev = 1000;
  lm.minimize(x);

  transformation_t transformation;
  transformation.col(3) = x.block<3,1>(0,0);
  transformation.block<3,3>(0,0) = math::cayley2rot(x.block<3,1>(3,0));
  return transformation;
}

}
}

opengv::transformation_t
opengv::absolute_pose::optimize_nonlinear( const AbsoluteLineAdapterBase & adapter )
{
  Indices idx(adapter.getNumberCorrespondences());
  return optimize_nonlinear(adapter,idx);
}

opengv::transformation_t
opengv::absolute_pose::optimize_nonlinear(
    const AbsoluteLineAdapterBase & adapter,
    const std::vector<int> & indices )
{
  Indices idx(indices);
  return optimize_nonlinear(adapter,idx);
}



// // optimize with point and with line
// 
// namespace opengv
// {
// namespace absolute_pose
// {
// 
// struct OptimizeNonlinearFunctor3 : OptimizationFunctor<double>
// {
//   const AbsoluteLineAdapterBase & _adapter;
//   const Indices & _indices;
//   const Indices & _line_indices;
// 
//   OptimizeNonlinearFunctor3(
//       const AbsoluteLineAdapterBase & adapter,
//       const Indices & indices ,const Indices & line_indices ) :
//       OptimizationFunctor<double>(6,indices.size()),
//       _adapter(adapter),
//       _indices(indices),
//       _line_indices(line_indices){}
// 
//   int operator()(const VectorXd &x, VectorXd &fvec) const
//   {
//     assert( x.size() == 6 );
//     assert( (unsigned int) fvec.size() == _indices.size());
// 
//     //compute the current position
//     translation_t translation = x.block<3,1>(0,0);
//     cayley_t cayley = x.block<3,1>(3,0);
//     rotation_t rotation = math::cayley2rot(cayley);
// 
//     //compute inverse transformation
//     transformation_t inverseSolution;
//     inverseSolution.block<3,3>(0,0) = rotation.transpose();
//     inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*translation;
// 
//     Eigen::Matrix<double,4,1> p_hom;
//     p_hom[3] = 1.0;
// 
//     for(size_t i = 0; i < _indices.size(); i++)
//     {
//       //get point in homogeneous form
//       p_hom.block<3,1>(0,0) = _adapter.getPoint(_indices[i]);
// 
//       //compute the reprojection (this is working for both central and
//       //non-central case)
//       point_t bodyReprojection = inverseSolution * p_hom;
//       point_t reprojection = _adapter.getCamRotation(_indices[i]).transpose() *
//           (bodyReprojection - _adapter.getCamOffset(_indices[i]));
//       reprojection = reprojection / reprojection.norm();
// 
//       //compute the score
//       double factor = 1.0;
//       fvec[i] = factor *
//           (1.0 -
//           (reprojection.transpose() * _adapter.getBearingVector(_indices[i])));
//     }
// 
//     return 0;
//   }
// };
// 
// transformation_t optimize_nonlinear(
//     const AbsoluteLineAdapterBase & adapter,
//     const Indices & indices,
// 	const Indices & line_indices )
// {
//   const int n=6;
//   VectorXd x(n);
// 
//   x.block<3,1>(0,0) = adapter.gett();
//   x.block<3,1>(3,0) = math::rot2cayley(adapter.getR());
// 
//   OptimizeNonlinearFunctor3 functor( adapter, indices, line_indices );
//   NumericalDiff<OptimizeNonlinearFunctor3> numDiff(functor);
//   LevenbergMarquardt< NumericalDiff<OptimizeNonlinearFunctor3> > lm(numDiff);
// 
//   lm.resetParameters();
//   lm.parameters.ftol = 1.E1*NumTraits<double>::epsilon();
//   lm.parameters.xtol = 1.E1*NumTraits<double>::epsilon();
//   lm.parameters.maxfev = 1000;
//   lm.minimize(x);
// 
//   transformation_t transformation;
//   transformation.col(3) = x.block<3,1>(0,0);
//   transformation.block<3,3>(0,0) = math::cayley2rot(x.block<3,1>(3,0));
//   return transformation;
// }
// 
// }
// }
// 
// 
// opengv::transformation_t
// opengv::absolute_pose::optimize_nonlinear(
//     const AbsoluteLineAdapterBase & adapter,
//     const std::vector<int> & indices,
// 	const std::vector<int> & line_indices)
// {
//   Indices idx(indices);
//   indices lidx(line_indices);
//   return optimize_nonlinear(adapter,idx,lidx);
// }

