#include <iostream>
#include <vector>
#include "smoothing.h"
#include"meshquality.h"
#include"writefile.h"
#include "equation.h"
#include "NegaOptimize.h"
//#include <graphics.h>      // 引用图形库头文件
#include <conio.h>
using namespace std;
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <OpenMesh/Tools/Smoother/LaplaceSmootherT.hh>
// ---------------------openvolumemesh
// Include vector classes
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

// Make some typedefs to facilitate your life
typedef OpenVolumeMesh::Geometry::Vec3d         Vec3d;
typedef OpenVolumeMesh::Geometry::Vec3f         Vec3f;
typedef OpenVolumeMesh::GeometryKernel<Vec3f>   PolyhedralMeshV3f;


#include<OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d MyTetMesh;


int main2(int _argc, char** _argv) {

	// Create mesh object
	PolyhedralMeshV3f myMesh;

	// Add eight vertices
	OpenVolumeMesh::VertexHandle v0 = myMesh.add_vertex(Vec3f(-1.0, 0.0, 0.0));
	OpenVolumeMesh::VertexHandle v1 = myMesh.add_vertex(Vec3f(0.0, 0.0, 1.0));
	OpenVolumeMesh::VertexHandle v2 = myMesh.add_vertex(Vec3f(1.0, 0.0, 0.0));
	OpenVolumeMesh::VertexHandle v3 = myMesh.add_vertex(Vec3f(0.0, 0.0, -1.0));
	OpenVolumeMesh::VertexHandle v4 = myMesh.add_vertex(Vec3f(0.0, 1.0, 0.0));

	std::vector<OpenVolumeMesh::VertexHandle> vertices;

	// Add faces
	vertices.push_back(v0); vertices.push_back(v1); vertices.push_back(v4);
	OpenVolumeMesh::FaceHandle f0 = myMesh.add_face(vertices);

	vertices.clear();
	vertices.push_back(v1); vertices.push_back(v2); vertices.push_back(v4);
	OpenVolumeMesh::FaceHandle f1 = myMesh.add_face(vertices);

	vertices.clear();
	vertices.push_back(v0); vertices.push_back(v1); vertices.push_back(v2);
	OpenVolumeMesh::FaceHandle f2 = myMesh.add_face(vertices);

	vertices.clear();
	vertices.push_back(v0); vertices.push_back(v4); vertices.push_back(v2);
	OpenVolumeMesh::FaceHandle f3 = myMesh.add_face(vertices);

	vertices.clear();
	vertices.push_back(v0); vertices.push_back(v4); vertices.push_back(v3);
	OpenVolumeMesh::FaceHandle f4 = myMesh.add_face(vertices);

	vertices.clear();
	vertices.push_back(v2); vertices.push_back(v3); vertices.push_back(v4);
	OpenVolumeMesh::FaceHandle f5 = myMesh.add_face(vertices);

	vertices.clear();
	vertices.push_back(v0); vertices.push_back(v2); vertices.push_back(v3);
	OpenVolumeMesh::FaceHandle f6 = myMesh.add_face(vertices);

	std::vector<OpenVolumeMesh::HalfFaceHandle> halffaces;

	// Add first tetrahedron
	halffaces.push_back(myMesh.halfface_handle(f0, 1));
	halffaces.push_back(myMesh.halfface_handle(f1, 1));
	halffaces.push_back(myMesh.halfface_handle(f2, 0));
	halffaces.push_back(myMesh.halfface_handle(f3, 1));
	myMesh.add_cell(halffaces);

	// Add second tetrahedron
	halffaces.clear();
	halffaces.push_back(myMesh.halfface_handle(f4, 1));
	halffaces.push_back(myMesh.halfface_handle(f5, 1));
	halffaces.push_back(myMesh.halfface_handle(f3, 0));
	halffaces.push_back(myMesh.halfface_handle(f6, 0));
	myMesh.add_cell(halffaces);

	// Print positions of vertices to std out
	for (OpenVolumeMesh::VertexIter v_it = myMesh.vertices_begin();
		v_it != myMesh.vertices_end(); ++v_it) {

		std::cout << "Position of vertex " << v_it->idx() << ": " <<
			myMesh.vertex(*v_it) << std::endl;
	}

	return 0;
}


void main(){
	MyMesh mymesh1;
	MyMesh mymesh2;
	MyMesh mymesh3;
	MyMesh mymesh4;
	//MyMesh mymesh3;
	string filename = "sresult2.stl";
	OpenMesh::IO::read_mesh(mymesh2, filename);
	OpenMesh::IO::read_mesh(mymesh3, filename);

	for (int i = -1; i < 10; i++)
	{
		if (i >= 0){
			angleBasedSmoothing(mymesh2, 1);
			GetMe(mymesh3, 1);
		}
		vector<double> qualities2 = computeMeshAngleQuality(mymesh2);
		vector<double> qualities3 = computeMeshAngleQuality(mymesh3);
		vector<double> qualities5 = computeMeshIdealElementQuality(mymesh2);
		vector<double> qualities6 = computeMeshIdealElementQuality(mymesh3);
		//vector<double> qualities4 = computeMeshAngleQuality(mymesh4);
		//输出网格质量
		cout << "angleBasedSmoothing    第" << i+1 << "次：";
		outMeshAngleQuality(qualities2);
		outMeshIdealElementQuality(qualities5);
		cout << "GetMe      第" << i+1 << "次：";
		outMeshAngleQuality( qualities3);
		outMeshIdealElementQuality(qualities6);
		//cout << "SmartSplitAngle 第" << i+1 << "次：";
		//outMeshAngleQuality(qualities4);
		cout << endl;
		//bool out1 = OpenMesh::IO::write_mesh(mymesh1, "1.stl");
	}

	bool out2 = OpenMesh::IO::write_mesh(mymesh2, "ange2000.stl");
	bool out3 = OpenMesh::IO::write_mesh(mymesh3, "get2000.stl");
}


void main3(){
	// 测试比较优化算法
	MyMesh mesh0;
	MyMesh mesh1;
	MyMesh mesh2;
	MyMesh mesh3;
	bool result0 = OpenMesh::IO::read_mesh(mesh0, "getsresult3.stl");
	bool result1 = OpenMesh::IO::read_mesh(mesh1, "myoptsresult3.stl");
	//bool result2 = OpenMesh::IO::read_mesh(mesh2, "foot.stl");
	//bool result3 = OpenMesh::IO::read_mesh(mesh3, "foot_eas.stl");
	if (result1 == false || result1 == false)
	{
		return;
	}
	vector<double> qualities0 = computeMeshAngleQuality(mesh0);
	vector<double> qualities1 = computeMeshAngleQuality(mesh1);
	vector<double> qualities2 = computeMeshIdealElementQuality(mesh0);
	vector<double> qualities3 = computeMeshIdealElementQuality(mesh1);
	//vector<double> qualities3 = computeMeshIdealElementQuality(mesh3);
	outMeshAngleQuality(qualities0);
	outMeshIdealElementQuality(qualities2);
	outMeshAngleQuality(qualities1);
	outMeshIdealElementQuality(qualities3);
}

void main12()
{
	MyMesh mesh;
	bool result0 = OpenMesh::IO::read_mesh(mesh, "sresult3.stl");
	double aveang = computeMeshAveang(mesh);
	cout << aveang << endl;
	MyMesh mesh2;
	bool result2 = OpenMesh::IO::read_mesh(mesh2, "angsresult3.stl");
	double aveang2 = computeMeshAveang(mesh2);
	cout << aveang2 << endl;
	MyMesh mesh3;
	bool result3 = OpenMesh::IO::read_mesh(mesh3, "getsresult3.stl");
	double aveang3 = computeMeshAveang(mesh3);
	cout << aveang3 << endl;
	MyMesh mesh4;
	bool result4 = OpenMesh::IO::read_mesh(mesh4, "optsresult3.stl");
	double aveang4 = computeMeshAveang(mesh4);
	cout << aveang4 << endl;
	MyMesh mesh5;
	bool result5 = OpenMesh::IO::read_mesh(mesh5, "myoptsresult3.stl");
	double aveang5 = computeMeshAveang(mesh5);
	cout << aveang5 << endl;
}

void main13()
{
	MyMesh mesh;
	bool result0 = OpenMesh::IO::read_mesh(mesh, "sss3.stl");
	negprove(mesh,1);
	bool out1 = OpenMesh::IO::write_mesh(mesh, "sresult3.stl");

}

void main14()
{
	readopt3("opt3.txt");
	readopt4("opt4.txt");
	readopt5("opt5.txt");
	readopt6("opt6.txt");
	readopt7("opt7.txt");
	readopt8("opt8.txt");
	readopt9("opt9.txt");
	MyMesh mesh1, mesh2, mesh3;
	bool result0 = OpenMesh::IO::read_mesh(mesh1, "e2000.stl");
	bool result1 = OpenMesh::IO::read_mesh(mesh2, "e2000.stl");
	bool result2 = OpenMesh::IO::read_mesh(mesh3, "e2000.stl");
	clock_t ts = clock();
	GetMe(mesh1, 10);
	clock_t te = clock();
	cout << 0 << " tests " << (te - ts) / 1000. << " sec." << endl;
	OpenMesh::IO::write_mesh(mesh1, "1111.stl");
	ts = clock();
	laplacianSmoothing(mesh2, 10);
	te = clock();
	cout << 1 << " tests " << (te - ts) / 1000. << " sec." << endl;
	OpenMesh::IO::write_mesh(mesh2, "11111.stl");
	ts = clock();
	NNSmoothing(mesh3, 10);
	te = clock();
	cout << 2 << " tests " << (te - ts) / 1000. << " sec." << endl;
	OpenMesh::IO::write_mesh(mesh3, "nne2000.stl");
	//bool out1 = OpenMesh::IO::write_mesh(mesh3, "shenjingwangluo.stl");
	cout << "hello" << endl;

}