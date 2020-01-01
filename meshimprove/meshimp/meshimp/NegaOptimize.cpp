#include "NegaOptimize.h"
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

void NegaOptimize(MyTetMesh &mesh, int iternum)
	{
		//循环次数
		for (int i = 0; i < iternum; i++){
			cout << i << endl;
			double num = 0;
			//循环所有节点
			for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++){
				int flag = 0;
				Vec3d cog;
				if (mesh.is_boundary(*v_it))
					continue;
				//节点周围所有单元
				for (auto vc_it = mesh.vc_iter(*v_it); vc_it.valid(); ++vc_it){
					vector<OpenVolumeMesh::VertexHandle> vhs;
					Vec3d target;
					for (auto cv_it = mesh.cv_iter(*vc_it); cv_it.valid(); ++cv_it){
						if (*cv_it == *v_it){
							vhs.push_back(*cv_it);
							target += mesh.vertex(*cv_it);
						}
					}
					target /= vhs.size();
					Vec3d v = target - mesh.vertex(*v_it);
					cog = mesh.vertex(*v_it) + 0.8*v;
				}
				mesh.set_vertex(*v_it, cog);
			}
		}
	}