#ifndef _MESHQUALITYFILE_
#define _MESHQUALITYFILE_

#include <iostream>
#include <vector>
#include "smoothing.h"
#include"meshquality.h"
#include"writefile.h"
using namespace std;
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
// ---------------------openvolumemesh
// Include vector classes
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
typedef OpenMesh::PolyMesh_ArrayKernelT<>  PolyMesh;

// Make some typedefs to facilitate your life
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;
typedef OpenVolumeMesh::Geometry::Vec3f         Vec3f;
typedef OpenVolumeMesh::GeometryKernel<Vec3f>   PolyhedralMeshV3f;


#include<OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d MyTetMesh;
using namespace std;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d MyTetMesh;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;



vector<double> computesizelenQuality(MyMesh& mesh);


vector<double> computeMeshIdealElementQuality(MyMesh& mesh);
vector<double> computeMeshIdealElementQuality(MyTetMesh& mesh);

vector<double> computeMeshAngleQuality(MyMesh& mesh);
vector<double> computeMeshAngleQuality(MyTetMesh& mesh);
vector<double> computeMeshAngleQuality(PolyMesh& mesh);

double computeMeshAveang(MyMesh& mesh);

//输出网格质量
void outMeshIdealElementQuality(vector<double> &qualities);
void outMeshAngleQuality(vector<double> &qualities);
void outFileMeshAngleQuality(string fname, vector<double> &qualities2);
void outFileMeshIdealElementQuality(string fname, vector<double> &qualities2);


double computeLocalAngleQuality(MyMesh::Point p, vector<MyMesh::Point> ring);

bool isImprovedLocally(MyMesh &mesh, MyMesh::VertexHandle vh, MyMesh::Point &cog);
bool isImprovedLocally(MyTetMesh &mesh, OpenVolumeMesh::VertexHandle vh, Vec3d &cog);

#endif