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
#endif // !_SMOOTHINGFILE_

