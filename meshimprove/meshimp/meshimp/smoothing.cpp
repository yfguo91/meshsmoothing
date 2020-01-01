#include"smoothing.h"
#include<vector>
#include<iostream>
#include"equation.h"
#include <graphics.h>      // 引用图形库头文件
#include <conio.h>
#include "meshquality.h"
#define PI 3.14159265
using namespace std;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;

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
			//如果是边界点的话直接跳过
			if (mesh.is_boundary(*v_it))
				continue;
			weights = 0.0;
			//找出ring的所有点来
			for (auto vv_it = mesh.vv_ccwbegin(*v_it); vv_it != mesh.vv_ccwend(*v_it); ++vv_it)
			{
				oldpoints.push_back(mesh.point(*vv_it));
			}
			//后面减前面可以计算出边来
			oldpoints.push_back(oldpoints[0]);
			//计算新的节点集
			for (int i = 0; i < oldpoints.size()-1; i++)
			{
				//每次两个点进行处理
				MyMesh::Point p1 = me;
				MyMesh::Point p2 = oldpoints[i];
				MyMesh::Point p3 = oldpoints[i+1];
				MyMesh::Point y1, y2, y3;
				MyMesh::Point z1,z2,z3;
				//计算单元质量
				MyMesh::Scalar qbefore = calq(p1, p2, p3);
				MyMesh::Scalar sbefore = cals(p1, p2, p3);
				//进行单元的转化
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
				//计算新的节点
				MyMesh::Point newpoint = (p1 + p2 + p3) / 3. + sqrt(sbefore / safter)*(z1 - (p1 + p2 + p3) / 3.);
				newpoints.push_back(newpoint);
				//newpointweight.push_back(1. / ((pi - pj).length()*(pi - pj).length()));
				//newpointweight.push_back(0.5+1-qbefore);
				newpointweight.push_back(qafter / qbefore);
			}
			//计算权重
			for (int i = 0; i < newpointweight.size(); i++)
			{
				weights += newpointweight[i];
			}
			//计算最终节点位置
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

