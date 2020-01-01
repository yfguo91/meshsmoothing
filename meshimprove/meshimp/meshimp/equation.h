#include<iostream>
#include<vector>
#include<iostream>
using namespace std;
// Include vector classes
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;
using namespace std;
Vec3d calc(double matrix1[3][4]);
Vec3d calc(double matrix[2][3]);
Vec3d Solve_OpenNL(double matrix[2][3]);