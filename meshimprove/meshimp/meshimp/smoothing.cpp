#include"smoothing.h"
#include<vector>
#include<iostream>
#include"equation.h"
#include <graphics.h>      // ����ͼ�ο�ͷ�ļ�
#include <conio.h>
#include "meshquality.h"
#define PI 3.14159265
using namespace std;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;
#include <Dense>
using namespace Eigen;
void angleBasedSmoothing(MyMesh &mesh, int iternum)
{
	for (int i = 0; i < iternum; i++)
	{
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			MyMesh::Point cog;
			MyMesh::Scalar valence;
			if (mesh.is_boundary(*v_it))
				continue;
			cog[0] = cog[1] = cog[2] = valence = 0.0;

			for (auto vv_it = mesh.vv_begin(*v_it); vv_it.is_valid(); ++vv_it)
			{
				auto vvp_it = vv_it, vvm_it = vv_it, vvn_it = vv_it;
				vvm_it++;
				if (vvm_it == mesh.vv_end(*v_it)){
					vvm_it = mesh.vv_begin(*v_it);
					vvn_it = vvm_it; vvn_it++;
				}
				else {
					vvn_it = vvm_it; vvn_it++;
					if (vvn_it == mesh.vv_end(*v_it)){
						vvn_it = mesh.vv_begin(*v_it);
					}
				}
				MyMesh::Normal v, v1, v2;
				v = mesh.point(*v_it) - mesh.point(*vvm_it);
				v1 = mesh.point(*vvp_it) - mesh.point(*vvm_it);
				v2 = mesh.point(*vvn_it) - mesh.point(*vvm_it);

				MyMesh::Point ptest = mesh.point(*v_it);

				v.normalize();
				v1.normalize();
				v2.normalize();
				double beta = -(0.5*(acos(v1 | v) - acos(v2 | v)));

				cog[0] += (mesh.point(*vvm_it)[0] + (mesh.point(*v_it)[0] - mesh.point(*vvm_it)[0])*cos(beta) - (mesh.point(*v_it)[1] - mesh.point(*vvm_it)[1])*sin(beta));
				cog[1] += (mesh.point(*vvm_it)[1] + (mesh.point(*v_it)[0] - mesh.point(*vvm_it)[0])*sin(beta) + (mesh.point(*v_it)[1] - mesh.point(*vvm_it)[1])*cos(beta));
				++valence;
			}
			cog = cog / valence;
			mesh.set_point(*v_it, cog);
		}
	}
}

void smartangleBasedSmoothing(MyMesh &mesh, int iternum)
{
	for (int i = 0; i < iternum; i++)
	{
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			MyMesh::Point cog;
			MyMesh::Scalar valence;
			if (mesh.is_boundary(*v_it))
				continue;
			cog[0] = cog[1] = cog[2] = valence = 0.0;

			for (auto vv_it = mesh.vv_begin(*v_it); vv_it.is_valid(); ++vv_it)
			{
				auto vvp_it = vv_it, vvm_it = vv_it, vvn_it = vv_it;
				vvm_it++;
				if (vvm_it == mesh.vv_end(*v_it)){
					vvm_it = mesh.vv_begin(*v_it);
					vvn_it = vvm_it; vvn_it++;
				}
				else {
					vvn_it = vvm_it; vvn_it++;
					if (vvn_it == mesh.vv_end(*v_it)){
						vvn_it = mesh.vv_begin(*v_it);
					}
				}
				MyMesh::Normal v, v1, v2;
				v = mesh.point(*v_it) - mesh.point(*vvm_it);
				v1 = mesh.point(*vvp_it) - mesh.point(*vvm_it);
				v2 = mesh.point(*vvn_it) - mesh.point(*vvm_it);

				MyMesh::Point ptest = mesh.point(*v_it);

				v.normalize();
				v1.normalize();
				v2.normalize();
				double beta = -(0.5*(acos(v1 | v) - acos(v2 | v)));

				cog[0] += (mesh.point(*vvm_it)[0] + (mesh.point(*v_it)[0] - mesh.point(*vvm_it)[0])*cos(beta) - (mesh.point(*v_it)[1] - mesh.point(*vvm_it)[1])*sin(beta));
				cog[1] += (mesh.point(*vvm_it)[1] + (mesh.point(*v_it)[0] - mesh.point(*vvm_it)[0])*sin(beta) + (mesh.point(*v_it)[1] - mesh.point(*vvm_it)[1])*cos(beta));
				++valence;
			}
			cog = cog / valence;
			if (isImprovedLocally(mesh, *v_it, cog))
			mesh.set_point(*v_it, cog);
		}
	}
}
double calq(MyMesh::Point p0, MyMesh::Point p1, MyMesh::Point p2)
{
	MyMesh::Normal v0 = p1 - p0;
	MyMesh::Normal v1 = p2 - p0;
	MyMesh::Normal v2 = p1 - p2;
	MyMesh::Scalar area = (v0%v1).length();
	MyMesh::Scalar quality = 2 * sqrt(3) * area / (v0.length()*v0.length() + v1.length()*v1.length() + v2.length()*v2.length());
	return quality;
}
double cals(MyMesh::Point p0, MyMesh::Point p1, MyMesh::Point p2)
{
	MyMesh::Normal v0 = p1 - p0;
	MyMesh::Normal v1 = p2 - p0;
	MyMesh::Normal v2 = p1 - p2;
	MyMesh::Scalar area = (v0%v1).length()/2.;
	return area;
}
void GetMe(MyMesh &mesh, int iternum)
{
	for (int i = 0; i < iternum; i++)
	{
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			MyMesh::Point cog;
			vector<MyMesh::Point>  oldpoints;
			vector<MyMesh::Point>  newpoints;
			vector<MyMesh::Scalar>  newpointweight;
			MyMesh::Scalar weights = 0.;

			MyMesh::Scalar valence;
			MyMesh::Point me = mesh.point(*v_it);
			//����Ǳ߽��Ļ�ֱ������
			if (mesh.is_boundary(*v_it))
				continue;
			weights = 0.0;
			//�ҳ�ring�����е���
			for (auto vv_it = mesh.vv_ccwbegin(*v_it); vv_it != mesh.vv_ccwend(*v_it); ++vv_it)
			{
				oldpoints.push_back(mesh.point(*vv_it));
			}
			//�����ǰ����Լ��������
			oldpoints.push_back(oldpoints[0]);
			//�����µĽڵ㼯
			for (int i = 0; i < oldpoints.size()-1; i++)
			{
				//ÿ����������д���
				MyMesh::Point p1 = me;
				MyMesh::Point p2 = oldpoints[i];
				MyMesh::Point p3 = oldpoints[i+1];
				MyMesh::Point y1, y2, y3;
				MyMesh::Point z1,z2,z3;
				//���㵥Ԫ����
				MyMesh::Scalar qbefore = calq(p1, p2, p3);
				MyMesh::Scalar sbefore = cals(p1, p2, p3);
				//���е�Ԫ��ת��
				y1[0] = (p1[0] + p3[0]) / 2. + sqrt(3) / 2. *(p1[1] - p3[1]);
				y1[1] = (p1[1] + p3[1]) / 2. + sqrt(3) / 2. * (p3[0] - p1[0]);
				y1[2] = 0;
				y2[0] = (p2[0] + p1[0]) / 2. + sqrt(3) / 2. *(p2[1] - p1[1]);
				y2[1] = (p2[1] + p1[1]) / 2. + sqrt(3) / 2. * (p1[0] - p2[0]);
				y2[2] = 0;
				y3[0] = (p3[0] + p2[0]) / 2. + sqrt(3) / 2. *(p3[1] - p2[1]);
				y3[1] = (p3[1] + p2[1]) / 2. + sqrt(3) / 2. * (p2[0] - p3[0]);
				y3[2] = 0;
				z1[0] = (y2[0] + y1[0]) / 2. + sqrt(3) / 2. *(y2[1] - y1[1]);
				z1[1] = (y2[1] + y1[1]) / 2. + sqrt(3) / 2. * (y1[0] - y2[0]);
				z1[2] = 0;
				z2[0] = (y3[0] + y2[0]) / 2. + sqrt(3) / 2. *(y3[1] - y2[1]);
				z2[1] = (y3[1] + y2[1]) / 2. + sqrt(3) / 2. * (y2[0] - y3[0]);
				z2[2] = 0;
				z3[0] = (y1[0] + y3[0]) / 2. + sqrt(3) / 2. *(y1[1] - y3[1]);
				z3[1] = (y1[1] + y3[1]) / 2. + sqrt(3) / 2. * (y3[0] - y1[0]);
				z3[2] = 0;
				MyMesh::Scalar qafter = calq(z1, z2, z3);
				MyMesh::Scalar safter = cals(z1, z2, z3);
				//�����µĽڵ�
				MyMesh::Point newpoint = (p1 + p2 + p3) / 3. + sqrt(sbefore / safter)*(z1 - (p1 + p2 + p3) / 3.);
				newpoints.push_back(newpoint);
				//newpointweight.push_back(1. / ((pi - pj).length()*(pi - pj).length()));
				//newpointweight.push_back(0.5+1-qbefore);
				newpointweight.push_back(qafter / qbefore);
			}
			//����Ȩ��
			for (int i = 0; i < newpointweight.size(); i++)
			{
				weights += newpointweight[i];
			}
			//�������սڵ�λ��
			cog[0] = cog[1] = cog[2] = 0.0;
			for (int i = 0; i < newpoints.size(); i++)
			{
				cog += newpoints[i] * newpointweight[i] / weights;
			}
			mesh.set_point(*v_it, cog);
		}
	}
}

void negprove(MyMesh &mesh, int iternum)
{
	for (int i = 0; i < iternum; i++)
	{
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			MyMesh::Point cog;
			if (mesh.is_boundary(*v_it))
				continue;
			int kk = 0;
			cog = 2.*mesh.point(*v_it);
			for (auto vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
			{
				cog += mesh.point(*vv_it);
				kk++;
				if (kk == 2)
					break;
			}
			cog = cog / 4.;
			mesh.set_point(*v_it, cog);
		}
	}
}

//������������Ż�
void NNSmoothing(MyMesh &mesh, int iternum)
{
	for (int i = 0; i < iternum; i++)
	{
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			MyMesh::Point cog;
			MyMesh::Scalar valence;
			MyMesh::Point position;
			vector<MyMesh::Point>  oldpoints;
			if (mesh.is_boundary(*v_it))
				continue;
			position[2] = 0.0;
			cog[0] = cog[1] = cog[2] = valence = 0.0;
			//�ҳ�ring�����е���
			MyMesh::Scalar ringxmin = mesh.point(*v_it)[0];
			MyMesh::Scalar ringxmax = mesh.point(*v_it)[0];
			MyMesh::Scalar ringymin = mesh.point(*v_it)[1];
			MyMesh::Scalar ringymax = mesh.point(*v_it)[1];
			for (auto vv_it = mesh.vv_ccwbegin(*v_it); vv_it != mesh.vv_ccwend(*v_it); ++vv_it)
			{
				oldpoints.push_back(mesh.point(*vv_it));
				cog += mesh.point(*vv_it);
				++valence;
				if (ringxmin > mesh.point(*vv_it)[0]) ringxmin = mesh.point(*vv_it)[0];
				if (ringxmax < mesh.point(*vv_it)[0]) ringxmax = mesh.point(*vv_it)[0];
				if (ringymin > mesh.point(*vv_it)[1]) ringymin = mesh.point(*vv_it)[1];
				if (ringymax < mesh.point(*vv_it)[1]) ringymax = mesh.point(*vv_it)[1];
			}
			if (oldpoints.size()>9)
			{
				cog = cog / valence;
				mesh.set_point(*v_it, cog);
				continue;
			}
			//��һ����0-1֮��
			VectorXd vput(oldpoints.size());
			MyMesh::Scalar len;
			if (ringxmax - ringxmin > ringymax - ringymin)
				len = ringxmax - ringxmin;
			else
				len = ringymax - ringymin;
			for (int i = 0; i < oldpoints.size(); i++)
			{
				oldpoints[i][0] = (oldpoints[i][0] - ringxmin) / len;
				//vput[2 * i] = oldpoints[i][0];
				oldpoints[i][1] = (oldpoints[i][1] - ringymin) / len;
				//vput[2 * i + 1] = oldpoints[i][1];
			}
			//�Ż�
			Vector2d vout;
			if (oldpoints.size() == 3) vout = nnopt3(oldpoints);
			if (oldpoints.size() == 4) vout = nnopt4(oldpoints);
			if (oldpoints.size() == 5) vout = nnopt5(oldpoints);
			if (oldpoints.size() == 6) vout = nnopt6(oldpoints);
			if (oldpoints.size() == 7) vout = nnopt7(oldpoints);
			if (oldpoints.size() == 8) vout = nnopt8(oldpoints);
			if (oldpoints.size() == 9) vout = nnopt9(oldpoints);
			//ӳ���ȥ
			position[0] = vout[0];
			position[1] = vout[1];
			position[0] = position[0] * len + ringxmin;
			position[1] = position[1] * len + ringymin;

			mesh.set_point(*v_it, position);
		}
	}
}


//����ÿ�����ͻ����Ż�
//3
MatrixXd opt3h1(20, 6);
VectorXd opt3b1(20);
MatrixXd opt3h2(6, 20);
VectorXd opt3b2(6);
MatrixXd opt3pr(2, 6);
VectorXd opt3br(2);
Vector2d nnopt3(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(6);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt3h1*vput + opt3b1;
	for (int i = 0; i < 20; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt3h2*vput1 + opt3b2;
	for (int i = 0; i < 6; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt3pr*vput2 + opt3br;
	return vput3;
}
void readopt3(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 20; i++)
			{
				for (int j = 0; j < 6; j++)
				{
					ist >> opt3h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 20; i++)
			{
				ist >> opt3b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 6; i++)
			{
				for (int j = 0; j < 20; j++)
				{
					ist >> opt3h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 6; i++)
			{
				ist >> opt3b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 6; j++)
				{
					ist >> opt3pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt3br(i);
			}
		}
	}
}
//4
MatrixXd opt4h1(20, 8);
VectorXd opt4b1(20);
MatrixXd opt4h2(8, 20);
VectorXd opt4b2(8);
MatrixXd opt4pr(2, 8);
VectorXd opt4br(2);
Vector2d nnopt4(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(8);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt4h1*vput + opt4b1;
	for (int i = 0; i < 20; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt4h2*vput1 + opt4b2;
	for (int i = 0; i < 8; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt4pr*vput2 + opt4br;
	return vput3;
}
void readopt4(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 20; i++)
			{
				for (int j = 0; j < 8; j++)
				{
					ist >> opt4h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 20; i++)
			{
				ist >> opt4b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 20; j++)
				{
					ist >> opt4h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 8; i++)
			{
				ist >> opt4b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 8; j++)
				{
					ist >> opt4pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt4br(i);
			}
		}
	}
}
//5
MatrixXd opt5h1(20, 10);
VectorXd opt5b1(20);
MatrixXd opt5h2(10, 20);
VectorXd opt5b2(10);
MatrixXd opt5pr(2, 10);
VectorXd opt5br(2);
Vector2d nnopt5(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(10);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt5h1*vput + opt5b1;
	for (int i = 0; i < 20; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt5h2*vput1 + opt5b2;
	for (int i = 0; i < 10; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt5pr*vput2 + opt5br;
	return vput3;
}
void readopt5(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 20; i++)
			{
				for (int j = 0; j < 10; j++)
				{
					ist >> opt5h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 20; i++)
			{
				ist >> opt5b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 10; i++)
			{
				for (int j = 0; j < 20; j++)
				{
					ist >> opt5h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 10; i++)
			{
				ist >> opt5b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 10; j++)
				{
					ist >> opt5pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt5br(i);
			}
		}
	}
}
//6
MatrixXd opt6h1(24, 12);
VectorXd opt6b1(24);
MatrixXd opt6h2(12, 24);
VectorXd opt6b2(12);
MatrixXd opt6pr(2, 12);
VectorXd opt6br(2);
Vector2d nnopt6(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(12);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt6h1*vput + opt6b1;
	for (int i = 0; i < 24; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt6h2*vput1 + opt6b2;
	for (int i = 0; i < 12; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt6pr*vput2 + opt6br;
	return vput3;
}
void readopt6(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 24; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					ist >> opt6h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 24; i++)
			{
				ist >> opt6b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 24; j++)
				{
					ist >> opt6h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 12; i++)
			{
				ist >> opt6b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					ist >> opt6pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt6br(i);
			}
		}
	}
}
//7
MatrixXd opt7h1(28, 14);
VectorXd opt7b1(28);
MatrixXd opt7h2(14, 28);
VectorXd opt7b2(14);
MatrixXd opt7pr(2, 14);
VectorXd opt7br(2);
Vector2d nnopt7(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(14);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt7h1*vput + opt7b1;
	for (int i = 0; i < 28; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt7h2*vput1 + opt7b2;
	for (int i = 0; i < 14; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt7pr*vput2 + opt7br;
	return vput3;
}
void readopt7(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 28; i++)
			{
				for (int j = 0; j < 14; j++)
				{
					ist >> opt7h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 28; i++)
			{
				ist >> opt7b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 14; i++)
			{
				for (int j = 0; j < 28; j++)
				{
					ist >> opt7h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 14; i++)
			{
				ist >> opt7b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 14; j++)
				{
					ist >> opt7pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt7br(i);
			}
		}
	}
}
//8
MatrixXd opt8h1(32, 16);
VectorXd opt8b1(32);
MatrixXd opt8h2(16, 32);
VectorXd opt8b2(16);
MatrixXd opt8pr(2, 16);
VectorXd opt8br(2);
Vector2d nnopt8(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(16);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt8h1*vput + opt8b1;
	for (int i = 0; i < 32; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt8h2*vput1 + opt8b2;
	for (int i = 0; i < 16; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt8pr*vput2 + opt8br;
	return vput3;
}
void readopt8(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 32; i++)
			{
				for (int j = 0; j < 16; j++)
				{
					ist >> opt8h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 32; i++)
			{
				ist >> opt8b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 16; i++)
			{
				for (int j = 0; j < 32; j++)
				{
					ist >> opt8h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 16; i++)
			{
				ist >> opt8b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 16; j++)
				{
					ist >> opt8pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt8br(i);
			}
		}
	}
}
//9
MatrixXd opt9h1(30, 18);
VectorXd opt9b1(30);
MatrixXd opt9h2(18, 30);
VectorXd opt9b2(18);
MatrixXd opt9pr(2, 18);
VectorXd opt9br(2);
Vector2d nnopt9(vector<MyMesh::Point> ppoint)
{
	VectorXd vput(18);
	for (int i = 0; i < ppoint.size(); i++)
	{
		vput[2 * i] = ppoint[i][0];
		vput[2 * i + 1] = ppoint[i][1];
	}
	VectorXd vput1 = opt9h1*vput + opt9b1;
	for (int i = 0; i < 30; i++)
	{
		vput1[i] = 1. / (1 + exp(-vput1[i]));
	}
	VectorXd vput2 = opt9h2*vput1 + opt9b2;
	for (int i = 0; i < 18; i++)
	{
		vput2[i] = 1. / (1 + exp(-vput2[i]));
	}
	Vector2d vput3 = opt9pr*vput2 + opt9br;
	return vput3;
}
void readopt9(char* filename)
{
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		if (strcmp(buf, "hidden1.weight") == 0)
		{
			for (int i = 0; i < 30; i++)
			{
				for (int j = 0; j < 18; j++)
				{
					ist >> opt9h1(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden1.bias") == 0)
		{
			for (int i = 0; i < 30; i++)
			{
				ist >> opt9b1(i);
			}
		}
		if (strcmp(buf, "hidden2.weight") == 0)
		{
			for (int i = 0; i < 18; i++)
			{
				for (int j = 0; j < 30; j++)
				{
					ist >> opt9h2(i, j);
				}
			}
		}
		if (strcmp(buf, "hidden2.bias") == 0)
		{
			for (int i = 0; i < 18; i++)
			{
				ist >> opt9b2(i);
			}
		}
		if (strcmp(buf, "predict.weight") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 18; j++)
				{
					ist >> opt9pr(i, j);
				}
			}
		}
		if (strcmp(buf, "predict.bias") == 0)
		{
			for (int i = 0; i < 2; i++)
			{
				ist >> opt9br(i);
			}
		}
	}
}


