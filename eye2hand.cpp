// eye2hand.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fstream>
#include <vector>
#include<fstream>
#include<cmath>
#include <Eigen/SVD>
using namespace std;
using namespace Eigen;
typedef vector<MatrixXd>poses;
void UseEigen()
{
	Eigen::VectorXd a(4);
	a << 1, 2, 3, 1;

	Eigen::MatrixXd T(4, 4);
	T << 1, 0, 0, 1,
		0, 1, 0, 0.2,
		0, 0, 1, 0.1,
		0, 0, 0, 1;

	//cout << endl << (T * a) << endl;
}

void EigenQuat()
{
	Eigen::Quaterniond q(0.5, 0.5, 0, 0); // Eigen quaternion is in the form of [ow, ox, oy, oz]
	cout << endl << q.coeffs() << endl;

	Eigen::VectorXd vec(3);
	vec << 1, 2, 3;
	cout << endl << q._transformVector(vec) << endl;

	Eigen::MatrixXd T = q.matrix();
	cout << endl << (T * vec) << endl;
}
int CsvRead(MatrixXd& outputMatrix, const string& fileName, const streamsize dPrec) {
	ifstream inputData;
	inputData.open(fileName);
	cout.precision(dPrec);
	if (!inputData)
		return -1;
	string fileline, filecell;
	unsigned int  noOfRows = 0, noOfCols = 0;
	//确定输入矩阵行列
	while (getline(inputData, fileline)) {
		noOfCols = 0;
		stringstream linestream(fileline);
		while (getline(linestream, filecell, ',')) {
			try {
				stod(filecell);
			}
			catch (...) {
				return -1;
			}
			noOfCols++;
		}
		noOfRows++;
	}
	inputData.close();
	outputMatrix.resize(noOfRows, noOfCols);
	inputData.open(fileName);
	noOfRows = 0;
	while (getline(inputData, fileline)) {
		noOfCols = 0;
		stringstream linestream(fileline);
		while (getline(linestream, filecell, ',')) {
			outputMatrix(noOfRows, noOfCols++) = stod(filecell);
		}
		noOfRows++;
	}
	return 0;
}
//将矢量转化为反对称矩阵
 Matrix3d skew(Vector3d u)
{
	Matrix3d u_hat =MatrixXd::Zero(3, 3);
	u_hat(0, 1) = u(2);
	u_hat(1, 0) = -u(2);
	u_hat(0, 2) = -u(1);
	u_hat(2, 0) = u(1);
	u_hat(1, 2) = u(0);
	u_hat(2, 1) = -u(0);

	return u_hat;
}
//将<平移，四元数>转化为齐次矩阵形式
void Vector2Homogeneous(MatrixXd source, poses&target)
{
	for (int i = 0; i < source.rows(); i++) {
		Eigen::Vector3d translation_base = source.block<1, 3>(i, 0);
		Eigen::VectorXd rotation_base = source.block<1, 4>(i, 3);
		Eigen::Quaterniond q_rotaion_base(rotation_base[3], rotation_base[0], rotation_base[1], rotation_base[2]);
		Eigen::Matrix3d Rotation_base = q_rotaion_base.matrix();
		Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
		T.block<3,3>(0,0)=Rotation_base;
		T.block<3,1>(0,3)=translation_base;
		target.push_back(T);
	}
}
MatrixXd SvdInverse(MatrixXd A)
{
	JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);//M=USV*
	double  pinvtoler = 1.e-6; //tolerance
	int row = A.rows();
	int col = A.cols();
	int k = min(row, col);
	MatrixXd X = MatrixXd::Zero(col, row);
	MatrixXd singularValues_inv = svd.singularValues();//奇异值
	MatrixXd singularValues_inv_mat = MatrixXd::Zero(col, row);
	for (long i = 0; i < k; ++i) {
		if (singularValues_inv(i) > pinvtoler)
			singularValues_inv(i) = 1.0 / singularValues_inv(i);
		else singularValues_inv(i) = 0;
	}
	for (long i = 0; i < k; ++i)
	{
		singularValues_inv_mat(i, i) = singularValues_inv(i);
	}
	X = (svd.matrixV())*(singularValues_inv_mat)*(svd.matrixU().transpose());//X=VS+U*

	return X;
}
MatrixXd Solve_Axxb(const poses AA,const poses BB) {
	Matrix4d Hx = Matrix4d::Identity(); //待求解标定矩阵
	int num = AA.size();//观察数据组个数
	cout << "数据对个数:" << num << endl;
    //MatrixXd output=MatrixXd::Identity();
	MatrixXd A = MatrixXd::Zero(num * 3, 3); // solve Ax=B for rotation
	VectorXd B = MatrixXd::Zero(num * 3, 1);
	//solve rotation
	for (size_t i = 0; i < AA.size(); i++) {
		//将旋转矩阵转化为旋转轴形式
		AngleAxisd rgij, rcij;
		rgij.fromRotationMatrix(AA[i].block<3,3>(0,0));
		rcij.fromRotationMatrix(BB[i].block<3,3>(0,0));
		double theta_gij = rgij.angle();
		double theta_cij = rcij.angle();
		Vector3d rngij = rgij.axis();
		Vector3d rncij = rcij.axis();
		AngleAxisd Pgij(2*sin(theta_gij / 2),rngij);
		AngleAxisd Pcij (2 * sin(theta_cij / 2),rncij);
		
		Vector3d pgij_vector = Pgij.axis()*Pgij.angle();
		Vector3d pcij_vector = Pcij.axis()*Pcij.angle();
		Vector3d vec_sum_axis = pgij_vector + pcij_vector;
		Matrix3d S = skew(vec_sum_axis);
		Vector3d vec_sub_axis = pcij_vector - pgij_vector;
				
		
		//cout << S << endl;
		//AngleAxisd sub_axis(Pcij.matrix()*Pgij.matrix().inverse());
		
		//Vector3d vec_sub_axis =sub_axis.axis()*sub_axis.angle();
		
		//A(3 * i - 2 : 3 * i, 1 : 3) = S;
		//B.block<>(3 * i - 2 : 3 * i, 1) = vec_sub_axis;
		A.block<3, 3>(3 * i, 0) = S;
		//cout << A.block<3, 3>(3 * i, 0) << endl;
		B.block<3, 1>(3 * i, 0) = vec_sub_axis;
	}
	//求解旋转
	VectorXd Pce1 = SvdInverse(A)*B;
	VectorXd Pce = 2 * Pce1 / sqrt(1 + Pce1.norm()*Pce1.norm());
	Matrix3d I=Matrix3d::Identity();
	double np2 = Pce.norm()*Pce.norm();
	MatrixXd Rx = (1 - np2 * 0.5)*I + 0.5*(Pce*Pce.transpose() + sqrt(4 - np2)*skew(Pce));
	cout<<"旋转矩阵是否正交:" <<Rx * Rx.transpose()<<endl;
	//求解平移
	A.setZero();
	B.setZero();
	for (size_t i = 0; i < AA.size(); i++) {
		Vector3d T_A, T_B;
		T_A = AA[i].block<3, 1>(0, 3);
		T_B = BB[i].block<3, 1>(0, 3);
		MatrixXd R_A = AA[i].block<3, 3>(0, 0);
		A.block<3, 3>(3 * i, 0) = R_A - I;
	    B.block<3, 1>(3 * i, 0) = Rx*T_B-T_A;
	}
	VectorXd Tx = SvdInverse(A)*B;
	Hx.block<3, 3>(0, 0) = Rx;
	Hx.block<3, 1>(0, 3) = Tx;
	return Hx;
	}
MatrixXd CalibrateHand2Eye(const poses HA, const poses HB) {
	const int n = HA.size();
	poses VA, VB;
	//求解相邻两次运动末端-末端，相机-相机之间变换
	for (int i = 0; i < n-1; i++)
	{
		
			MatrixXd A = HA[i+1] * HA[i].inverse();
			MatrixXd B = HB[i+1] * HB[i].inverse();
			VA.push_back(A);
			VB.push_back(B);
		
	}
	MatrixXd H = Solve_Axxb(VA, VB);
	return H;
}

int main()
{
	//UseEigen();
	//EigenQuat();
	string tool2base = "Data/T_b_t.txt";
	string cal2cam = "Data/T_c_cal.txt";
	MatrixXd Tool2Base, Camera2Cal;
	int errcode1 = CsvRead(Tool2Base, tool2base, 8); //读取csv矩阵,精度为小数点后8位
	int errcode2 = CsvRead(Camera2Cal,cal2cam,8);
	if (errcode1 == 0 && errcode2 == 0) {
		if ((Tool2Base.cols() != Camera2Cal.cols()) && (Tool2Base.rows() != Camera2Cal.rows()))
			cout << "source data's size should be equal" << endl;
		if (Tool2Base.rows() < 3)
			cout << "at least three motions are needed" << endl;
		//cout <<"Data/T_b_t.txt:"<<Tool2Base<<endl<<endl;
		//cout << "Data/T_c_cal.txt:" << Camera2Cal << endl;
		poses A, B;
		Vector2Homogeneous(Tool2Base, A);
		Vector2Homogeneous(Camera2Cal, B);
		MatrixXd sol = CalibrateHand2Eye(A, B);
		cout << "标定矩阵:" << sol;
	}
	else
		cout << "read raw csv data wrong!";
	return 0;

}

