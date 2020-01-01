#ifndef _WRITEFILE_
#define _WRITEFILE_
#include <iostream>
#include <vector>
using namespace std;
// -------------------- OpenMesh
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

// -------------------- OpenVolumeMesh
#include<OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d MyTetMesh;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;
// ----------------------------------------------------------------------------

//用于四边形
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
typedef OpenMesh::PolyMesh_ArrayKernelT<>  PolyMesh;

void writeVTK(MyMesh &mesh, char* filename);
void writeVTKWithFixedBoundery(MyMesh &mesh, char* filename);

void readEleFile(MyTetMesh &mytet, string filename);
//输入四面体vtk格式
bool readTetVTK(MyTetMesh &mesh, string fname);

//输入四边形网格
void readQuadVTK(PolyMesh &mesh, char* filename);

//输出四边形网格
void writeQuadVTK(PolyMesh &mesh, char* filename);


//输入三角形网格
void readTriVTK(MyMesh &mesh, char* filename);

//输出三角形网格
void writeTriVTK(MyMesh &mesh, char* filename);

//输入三角形网格
void readTritxt(MyMesh &mesh, char* filename1, char* filename2);

void readliu(MyMesh &mesh, char* filename);

#endif