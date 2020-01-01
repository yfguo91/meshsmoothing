#include<OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d MyTetMesh;
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;
void NegaOptimize(MyTetMesh &mesh, int iternum);