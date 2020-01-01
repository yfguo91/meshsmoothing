/**
* notice:时间仓促，仅仅实现功能，方便使用，代码质量不可参考！！！
*/
#include<iostream>
#include"equation.h"
#include<vector>
#include<iostream>
#include"OpenNL_psm.h"

using namespace std;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;

Vec3d calc1(double matrix1[3][4]){
	double matrix[4][5];
	Vec3d point;
	for (int i = 1; i<4; i++){
		for (int j = 1; j<5; j++){
			matrix[i][j] = matrix1[i-1][j-1];
		}
	}
	double base_D = matrix[1][1] * matrix[2][2] * matrix[3][3] + matrix[2][1] * matrix[3][2] * matrix[1][3] + matrix[3][1] * matrix[1][2] * matrix[2][3];//计算行列式 
	base_D = base_D - (matrix[1][3] * matrix[2][2] * matrix[3][1] + matrix[1][1] * matrix[2][3] * matrix[3][2] + matrix[1][2] * matrix[2][1] * matrix[3][3]);

	if (abs(base_D) >= 0.000001){
		double     x_D = matrix[1][4] * matrix[2][2] * matrix[3][3] + matrix[2][4] * matrix[3][2] * matrix[1][3] + matrix[3][4] * matrix[1][2] * matrix[2][3];
		x_D = x_D - (matrix[1][3] * matrix[2][2] * matrix[3][4] + matrix[1][4] * matrix[2][3] * matrix[3][2] + matrix[1][2] * matrix[2][4] * matrix[3][3]);
		double     y_D = matrix[1][1] * matrix[2][4] * matrix[3][3] + matrix[2][1] * matrix[3][4] * matrix[1][3] + matrix[3][1] * matrix[1][4] * matrix[2][3];
		y_D = y_D - (matrix[1][3] * matrix[2][4] * matrix[3][1] + matrix[1][1] * matrix[2][3] * matrix[3][4] + matrix[1][4] * matrix[2][1] * matrix[3][3]);
		double     z_D = matrix[1][1] * matrix[2][2] * matrix[3][4] + matrix[2][1] * matrix[3][2] * matrix[1][4] + matrix[3][1] * matrix[1][2] * matrix[2][4];
		z_D = z_D - (matrix[1][4] * matrix[2][2] * matrix[3][1] + matrix[1][1] * matrix[2][4] * matrix[3][2] + matrix[1][2] * matrix[2][1] * matrix[3][4]);

		double x = x_D / base_D;
		double y = y_D / base_D;
		double z = z_D / base_D;
		point[0] = x;
		point[1] = y;
		point[2] = z;
		//cout << "[ x:" << x << "; y:" << y << "; z:" << z << " ]" << endl;
	}
	else{
		//cout << "【无解】";
		//        return DBL_MIN;
	}
	return point;
}
Vec3d calc(double matrix[3][4]){
	Vec3d point;
	nlNewContext();
	nlSolverParameteri(NL_NB_VARIABLES, 3);
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);

	nlBegin(NL_SYSTEM);
	nlBegin(NL_MATRIX);
	for (int i = 0; i < 3; i++) {
		nlBegin(NL_ROW);
		for (int j = 0; j < 3; j++) {
			nlCoefficient(j, matrix[i][j]);
		}
		nlRightHandSide(matrix[i][2]);
		nlEnd(NL_ROW);
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	nlSolve();
	for (int i = 0; i < 3; i++) {
		point[i] = nlGetVariable(i);
	}
	nlDeleteContext(nlGetCurrent());
	return point;
}
Vec3d calc(double matrix[2][3]){

	Vec3d point;
	point[1] = (matrix[0][0] * matrix[1][2] - matrix[1][0] * matrix[0][2]) / (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]);
	if (matrix[0][0]!=0)
		point[0] = (matrix[0][2]-matrix[0][1]*point[1]) / (matrix[0][0]);
	else if (matrix[1][0] != 0)
		point[0] = (matrix[1][2] - matrix[1][1] * point[1]) / (matrix[1][0]);
	else{
		system("pause");
		exit(1);
	}
	point[2] = 0;
	//cout << (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) << endl;
	//cout << point[0] << " " << point[1] << endl;
	return point;
}

Vec3d Solve_OpenNL(double matrix[2][3]) {
	Vec3d point;
	nlNewContext();
	nlSolverParameteri(NL_NB_VARIABLES, 2);
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);

	nlBegin(NL_SYSTEM);
	nlBegin(NL_MATRIX);
	for (int i = 0; i < 2; i++) {
		nlBegin(NL_ROW);
		for (int j = 0; j < 2; j++) {
			nlCoefficient(j, matrix[i][j]);
		}
		nlRightHandSide(matrix[i][2]);
		nlEnd(NL_ROW);
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	nlSolve();
	for (int i = 0; i < 2; i++) {
		point[i] = nlGetVariable(i);
	}
	nlDeleteContext(nlGetCurrent());
	return point;
}
/*
demo

2x-y+z=10;
3x+2y-z=16;
x+6y-z=28;
2 -1 1 10
3 2 -1 16
1 6 -1 28

output:input ok[ x:4.18182; y:5.09091; z:6.72727 ]
*/