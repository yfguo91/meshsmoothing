#ifndef _SMOOTHINGFILE_
#define _SMOOTHINGFILE_

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include<OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
typedef OpenMesh::PolyMesh_ArrayKernelT<>  PolyMesh;


#include"smoothing.h"
#include<vector>
#include<iostream>
#include"equation.h"
#include <graphics.h>      // 引用图形库头文件
#include <conio.h>
#include "meshquality.h"

typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d MyTetMesh;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;

void angleBasedSmoothing(MyMesh &mesh, int iternum);
void smartangleBasedSmoothing(MyMesh &mesh, int iternum);

void GetMe(MyMesh &mesh, int iternum);

void negprove(MyMesh &mesh, int iternum);

#include <Dense>
using namespace Eigen;
Vector2d nnopt3(vector<MyMesh::Point> ppoint);
Vector2d nnopt4(vector<MyMesh::Point> ppoint);
Vector2d nnopt5(vector<MyMesh::Point> ppoint);
Vector2d nnopt6(vector<MyMesh::Point> ppoint);
Vector2d nnopt7(vector<MyMesh::Point> ppoint);
Vector2d nnopt8(vector<MyMesh::Point> ppoint);
Vector2d nnopt9(vector<MyMesh::Point> ppoint);
void readopt3(char* filename);
void readopt4(char* filename);
void readopt5(char* filename);
void readopt6(char* filename);
void readopt7(char* filename);
void readopt8(char* filename);
void readopt9(char* filename);


void NNSmoothing(MyMesh &mesh, int iternum);





#endif // !_SMOOTHINGFILE_

